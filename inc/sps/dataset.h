#pragma once

#include "sps/desc.h"
#include "sps/overlay.h"
#include "sps/point.h"
#include "sps/sparse_coordinate.h"
#include "sps/thread_pool.h"
#include "sps/type_defs.h"
#include "sps/util.h"
#include <cassert>
#include <functional>
#include <mutex>
#include <string>


namespace sps
{


template <typename type_defs> class Dataset
{
    EXTRACT_TYPE_DEFS; // macro call

    using sparse_coord_t = SparseCoord<type_defs>;
    using overlay_t = AlignedPower2<Overlay<type_defs>>;
    using overlay_grid_t = typename Overlay<type_defs>::overlay_grid_t;
    using prefix_sum_grid_t = typename Overlay<type_defs>::prefix_sum_grid_t;
    using point_t = AlignedPower2<Point<type_defs>>;
    using points_t = Points<type_defs>;
    using desc_t = Desc<type_defs>;

    std::array<typename sparse_coord_t::Entry, D> vSparseCoords;

    typename sparse_coord_t::EntryArray xSparseCoordsDependantDimension;

    typename overlay_grid_t::template Entry<D> xOverlays;

    class DimIterator
    {
        typename points_t::EntryIterator xIt, xItEnd;
        size_t uiDimension;
        size_t uiBlockSize;

      public:
        DimIterator( typename points_t::EntryIterator xIt, typename points_t::EntryIterator xItEnd, size_t uiDimension,
                     size_t uiBlockSize )
            : xIt( xIt ), xItEnd( xItEnd ), uiDimension( uiDimension ), uiBlockSize( uiBlockSize )
        {}

        void operator++( )
        {
            for( size_t uiI = 0; uiI < uiBlockSize && xIt != xItEnd; uiI++ )
            {
                coordinate_t uiCurr = **this;
                while( xIt != xItEnd && uiCurr == **this )
                    ++xIt;
            }
        }

        const coordinate_t operator*( ) const
        {
            return xIt->vPos[ uiDimension ];
        }

        bool operator!=( const DimIterator& rOther ) const
        {
            return xIt != rOther.xIt;
        }

        friend std::ostream& operator<<( std::ostream& os, const DimIterator& rIt )
        {
            os << rIt.xIt;

            os << " d";
            os << rIt.uiDimension;

            return os;
        }
    };

    using points_it_t = typename points_t::points_vec_t::iterator;
    struct PointsBinComperator
    {
        const Dataset& rDataset;
        overlay_grid_t& rOverlays;
        sparse_coord_t& rSparseCoords;
        std::mutex& xLock;

        PointsBinComperator( const Dataset& rDataset, overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords,
                             std::mutex& xLock )
            : rDataset( rDataset ), rOverlays( rOverlays ), rSparseCoords( rSparseCoords ), xLock( xLock )
        {}

        bool comp( const point_t& a, const point_t& b ) const
        {
            coordinate_t uiA = rDataset.overlayIndex( rOverlays, rSparseCoords, a.vPos );
            coordinate_t uiB = rDataset.overlayIndex( rOverlays, rSparseCoords, b.vPos );
            if( uiA == uiB )
                return a.uiDescOffset < b.uiDescOffset;
            return uiA < uiB;
        }

        bool operator( )( const point_t& a, const point_t& b ) const
        {
            if constexpr( sparse_coord_t::THREADSAVE )
                return comp( a, b );
            else
            {
                std::lock_guard<std::mutex> xGuard( xLock );
                return comp( a, b );
            }
        }

        point_t min_value( ) const
        {
            return point_t(
                rSparseCoords.invSparse( rOverlays.posOf( rDataset.xOverlays.uiStartIndex, rDataset.xOverlays ),
                                         rDataset.vSparseCoords ),
                0 );
        };

        point_t max_value( ) const
        {
            coordinate_t uiIndex = rDataset.xOverlays.uiStartIndex + rOverlays.sizeOf( rDataset.xOverlays ) - 1;

            pos_t vPos = rOverlays.posOf( uiIndex, rDataset.xOverlays );
            pos_t vActPos = rDataset.actualFromGridPos( rSparseCoords, vPos );

            return point_t( vActPos, std::numeric_limits<size_t>::max( ) );
        };
    };
    points_sort_func_t<points_it_t, PointsBinComperator> sort_points_bin =
        points_sort_func_t<points_it_t, PointsBinComperator>( );

    struct PointsComperator
    {
        const size_t uiDim;

        PointsComperator( size_t uiDim ) : uiDim( uiDim )
        {}

        bool operator( )( const point_t& a, const point_t& b ) const
        {
            return a.vPos[ uiDim ] < b.vPos[ uiDim ];
        }

        point_t min_value( ) const
        {
            point_t xRet{ };
            return xRet;
        };

        point_t max_value( ) const
        {
            point_t xRet{ };
            xRet.vPos[ uiDim ] = std::numeric_limits<coordinate_t>::max( );
            return xRet;
        };
    };
    points_sort_func_t<points_it_t, PointsComperator> sort_points =
        points_sort_func_t<points_it_t, PointsComperator>( );

  public:
    Dataset( ) : vSparseCoords( ), xSparseCoordsDependantDimension( )
    {}

    typename sparse_coord_t::Entry
    makeSparseCoords( sparse_coord_t& rSparseCoords, points_t& vPoints, typename points_t::Entry xPoints, size_t uiDim,
                      progress_stream_t xProg, coordinate_t uiFixedStart = std::numeric_limits<coordinate_t>::max( ),
                      coordinate_t uiFixedEnd = std::numeric_limits<coordinate_t>::max( ) ) const
    {
        size_t uiNumCoords = 0;
        coordinate_t uiLast = std::numeric_limits<coordinate_t>::max( );
        sort_points( vPoints.vData.begin( ) + xPoints.uiStartIndex, vPoints.vData.begin( ) + xPoints.uiEndIndex,
                     PointsComperator( uiDim ) );
        vPoints.iterate(
            [ & ]( const point_t& xPoint ) {
                // @todo consider all coordinates of point
                coordinate_t uiCurr = xPoint.vPos[ uiDim ];
                if( uiCurr != uiLast )
                {
                    uiLast = uiCurr;
                    ++uiNumCoords;
                }
            },
            xPoints );
        // size_t uiNumPerDimension = (size_t)std::pow( uiNumCoords, 1.0 / (float)( D ) ); //
        size_t uiNumPerDimension = (size_t)std::pow( uiNumCoords, 1.0 / 2.0 ); // @todo why is this better?
        xProg << "generating " << uiNumPerDimension << " overlays in dimension " << uiDim << " for " << uiNumCoords
              << " different coordinates.\n";
        size_t uiBlockSize = std::max( 1ul, uiNumCoords / uiNumPerDimension );
        auto xStart = DimIterator( vPoints.cbegin( xPoints ), vPoints.cend( xPoints ), uiDim, uiBlockSize );
        auto xEnd = DimIterator( vPoints.cend( xPoints ), vPoints.cend( xPoints ), uiDim, uiBlockSize );
        if( uiFixedStart == std::numeric_limits<coordinate_t>::max( ) )
            return rSparseCoords.add( xStart, xEnd );
        return rSparseCoords.addStartEnd( xStart, xEnd, uiFixedStart, uiFixedEnd );
    }

#pragma GCC diagnostic push
// vPosTopRightActual not used with DEPENDANT_DIMENSION == false
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
    std::array<std::vector<coordinate_t>, D> getPredecessor( const overlay_grid_t& rOverlays,
                                                             const sparse_coord_t& rSparseCoords, pos_t vGridPos,
                                                             pos_t vPosBottomLeftActual, pos_t vPosTopRightActual,
                                                             progress_stream_t xProg ) const
    {

        auto vbHasPredecessor = hasPredecessor( rSparseCoords, vGridPos );
        std::array<std::vector<coordinate_t>, D> vPredecessors{ };
        for( size_t uiD = 0; uiD < D; uiD++ )
        {
            if( vbHasPredecessor[ uiD ] )
            {
                if constexpr( DEPENDANT_DIMENSION )
                {
#ifndef NDEBUG
                    coordinate_t uiIdx = rOverlays.indexOf( vGridPos, xOverlays );
#endif
                    if( uiD == 0 )
                    {
                        pos_t vItrPos = vPosBottomLeftActual;
                        assert( vItrPos[ uiD ] > 0 );
                        // move outside the current overlay
                        --vItrPos[ uiD ];
                        pos_t vItrTopRight;
                        do
                        {
                            coordinate_t uiItrIndex = overlayIndex( rOverlays, rSparseCoords, vItrPos );
                            vPredecessors[ uiD ].push_back( uiItrIndex );
                            pos_t vItrGridPos = rOverlays.posOf( uiItrIndex, xOverlays );
                            vItrTopRight = actualTopRightFromGridPos( rSparseCoords, vItrGridPos );

                            xProg << Verbosity( 2 ) << "uiD " << uiD << " vItrPos " << vItrPos << " uiItrIndex "
                                  << uiItrIndex << " vItrGridPos " << vItrGridPos << " vItrTopRight " << vItrTopRight
                                  << " vPosBottomLeftActual " << vPosBottomLeftActual << " vPosTopRightActual "
                                  << vPosTopRightActual << "\n";
                            assert( uiItrIndex != uiIdx );
                            assert( exists( rSparseCoords, vItrGridPos ) );

                            // move upwards
                            assert( vItrTopRight[ 1 ] > vItrPos[ 1 ] ||
                                    vItrTopRight[ 1 ] == std::numeric_limits<coordinate_t>::max( ) );
                            vItrPos[ 1 ] = vItrTopRight[ 1 ];
                        } while( vItrTopRight[ 1 ] < vPosTopRightActual[ 1 ] &&
                                 vItrTopRight[ 1 ] != std::numeric_limits<coordinate_t>::max( ) );

                        continue;
                    }
                }

                --vPosBottomLeftActual[ uiD ];

                vPredecessors[ uiD ].push_back( overlayIndex( rOverlays, rSparseCoords, vPosBottomLeftActual ) );

                xProg << Verbosity( 3 ) << "predecessor dim " << uiD << " pos " << vPosBottomLeftActual << " index "
                      << vPredecessors[ uiD ].back( ) << "\n";
                assert( exists( rSparseCoords, rOverlays.posOf( vPredecessors[ uiD ].back( ), xOverlays ) ) );

                ++vPosBottomLeftActual[ uiD ];
            }
            else
                xProg << Verbosity( 2 ) << "predecessor dim " << uiD << " nonexistant\n";
        }
        return vPredecessors;
    }
#pragma GCC diagnostic pop

    Dataset( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords, prefix_sum_grid_t& rPrefixSums,
             points_t& vPoints, typename points_t::Entry xPoints, progress_stream_t xProg )
        : vSparseCoords( ), xSparseCoordsDependantDimension( )
    {
        if( xPoints.uiEndIndex == xPoints.uiStartIndex )
            return;
        // generate the overall sparse coordinates
        Profiler xProfiler( "setting up overlay sparse coords" );
        ThreadPool xPool;
        xProg << Verbosity( 0 ) << "generating overlay grid.\n";
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            xProg << Verbosity( 0 );
            if constexpr( DEPENDANT_DIMENSION )
            {
                if( uiI == 1 )
                {
                    coordinate_t uiAbsoluteStart = std::numeric_limits<coordinate_t>::max( );
                    coordinate_t uiAbsoluteEnd = 0;
                    vPoints.iterate(
                        [ & ]( const point_t& xPoint ) {
                            uiAbsoluteStart = std::min( uiAbsoluteStart, xPoint.vPos[ 1 ] );
                            uiAbsoluteEnd = std::max( uiAbsoluteEnd, xPoint.vPos[ 1 ] );
                        },
                        xPoints );

                    coordinate_t uiNumBins = rSparseCoords.axisSize( vSparseCoords[ 0 ] );
                    typename points_t::Entry xCurrPoints{ };
                    xCurrPoints.uiEndIndex = xPoints.uiStartIndex;
                    vSparseCoords[ uiI ].uiStartIndex = std::numeric_limits<coordinate_t>::max( );
                    for( size_t uiBin = 0; uiBin < uiNumBins; uiBin++ )
                    {
                        // collect points for overlay uiI (points still sorted by dim 0)
                        xCurrPoints.uiStartIndex = xCurrPoints.uiEndIndex;
                        while( xCurrPoints.uiEndIndex < xPoints.uiEndIndex &&
                               rSparseCoords.replace( vPoints.vData[ xCurrPoints.uiEndIndex ].vPos[ 0 ],
                                                      vSparseCoords[ 0 ] ) == uiBin )
                            ++xCurrPoints.uiEndIndex;

                        auto xCurr = makeSparseCoords( rSparseCoords, vPoints, xCurrPoints, 1, xProg, uiAbsoluteStart,
                                                       uiAbsoluteEnd );
                        assert( rSparseCoords.axisSize( xCurr ) > 0 );
                        if( vSparseCoords[ uiI ].uiStartIndex == std::numeric_limits<coordinate_t>::max( ) ||
                            rSparseCoords.axisSize( vSparseCoords[ uiI ] ) < rSparseCoords.axisSize( xCurr ) )
                            vSparseCoords[ uiI ] = xCurr;
                        sparse_coord_t::append( xSparseCoordsDependantDimension, xCurr );
                        xProg << Verbosity( 1 );
                    }
                    assert( xCurrPoints.uiEndIndex == xPoints.uiEndIndex );

                    continue;
                }
            }
            vSparseCoords[ uiI ] = makeSparseCoords( rSparseCoords, vPoints, xPoints, uiI, xProg );
        }
        xProg << Verbosity( 2 );
        if( xProg.active( ) )
        {
            vSparseCoords[ 0 ].stream( std::cout << "sparse coords 0: ", rSparseCoords ) << std::endl;
            if constexpr( DEPENDANT_DIMENSION )
                xSparseCoordsDependantDimension.stream( std::cout << "sparse coords 1: ", rSparseCoords ) << std::endl;
            else
                vSparseCoords[ 1 ].stream( std::cout << "sparse coords 1: ", rSparseCoords ) << std::endl;
            for( size_t uiI = 2; uiI < D; uiI++ )
                vSparseCoords[ uiI ].stream( std::cout << "sparse coords " << uiI << ": ", rSparseCoords ) << std::endl;
        }


        // generate overlay grid
        xOverlays = rOverlays.add( rSparseCoords.axisSizes( vSparseCoords ) );

        // sort points so that they match the overlay grid order
        xProg << Verbosity( 0 ) << "sorting points into overlays.\n";
        /* Weirdly stxxl::sort seems to be parallel eventhough the access to stxxl::vector is required to be sequential.
         * This does not seem to be a problem for the vector that is currently sorted.
         * However, my PointsBinComperator accesses another stxxl::vector and therefore has to be guarded with a mutex.
         */
        std::mutex xLock;
        xProfiler.step( "sorting points into overlays" );
        sort_points_bin( vPoints.vData.begin( ) + xPoints.uiStartIndex, vPoints.vData.begin( ) + xPoints.uiEndIndex,
                         PointsBinComperator( *this, rOverlays, rSparseCoords, xLock ) );

        // generate all overlays
        coordinate_t uiNumTotal = rOverlays.sizeOf( xOverlays ); // @todo count considering non existant entries

        std::vector<typename points_t::Entry> vSplitPoints( uiNumTotal );
        for( coordinate_t uiI = 0; uiI < uiNumTotal; uiI++ )
        {
            vSplitPoints[ uiI ].uiStartIndex = uiI > 0 ? vSplitPoints[ uiI - 1 ].uiEndIndex : xPoints.uiStartIndex;
            vSplitPoints[ uiI ].uiEndIndex = vSplitPoints[ uiI ].uiStartIndex;
            // collect points for overlay uiI
            while( vSplitPoints[ uiI ].uiEndIndex < xPoints.uiEndIndex &&
                   overlayIndex( rOverlays, rSparseCoords, vPoints.vData[ vSplitPoints[ uiI ].uiEndIndex ].vPos ) ==
                       uiI + xOverlays.uiStartIndex )
                ++vSplitPoints[ uiI ].uiEndIndex;
        }
        assert( vSplitPoints[ uiNumTotal - 1 ].uiEndIndex == xPoints.uiEndIndex );

        xProg << Verbosity( 0 ) << "generating overlays.\n";
        xProfiler.step( "overlay gen loop" );

        auto xIterator = rOverlays.template genIterator<D>(
            xOverlays, [ & ]( size_t uiIdx, size_t uiDim,
                              // return vector of indices for successor overlays
                              const typename overlay_grid_t::template Entry<D>& xEntry ) {
                std::vector<size_t> vRet;
                auto vPos = overlay_grid_t::posOf( uiIdx, xEntry );

                bool bFits = true;
                for( size_t uiI = 0; uiI < D && bFits; uiI++ )
                {
                    ++vPos[ uiI ];
                    auto xPartialLock1 = rSparseCoords.xLockable.partialLock( );
                    if( DEPENDANT_DIMENSION && uiI == 1 )
                    {
                        if( vPos[ 1 ] < xEntry.vAxisSizes[ 1 ] && exists( rSparseCoords, vPos ) )
                            bFits = false;
                    }
                    else if( vPos[ uiI ] < xEntry.vAxisSizes[ uiI ] )
                        bFits = false;
                    --vPos[ uiI ];
                }
                if( bFits )
                {
                    vRet.push_back( std::numeric_limits<size_t>::max( ) ); // poison
                    return vRet;
                }

                ++vPos[ uiDim ];
                if( vPos[ uiDim ] < xEntry.vAxisSizes[ uiDim ] )
                {
                    auto xPartialLock1 = rSparseCoords.xLockable.partialLock( );
                    if( DEPENDANT_DIMENSION && uiDim == 0 )
                    {
                        --vPos[ uiDim ];

                        pos_t vBottomLeft = actualFromGridPos( rSparseCoords, vPos );
                        pos_t vTopRight = actualTopRightFromGridPos( rSparseCoords, vPos );

                        pos_t vItrPos = vBottomLeft;
                        vItrPos[ 0 ] = vTopRight[ 0 ];
                        pos_t vItrTopRight;
                        do
                        {
                            coordinate_t uiItrIndex = overlayIndex( rOverlays, rSparseCoords, vItrPos );
                            pos_t vItrGridPos = rOverlays.posOf( uiItrIndex, xOverlays );
                            vItrTopRight = actualTopRightFromGridPos( rSparseCoords, vItrGridPos );

                            if( vItrTopRight[ 1 ] <= vTopRight[ 1 ] && exists( rSparseCoords, vItrGridPos ) )
                                vRet.push_back( uiItrIndex );

                            // move upwards
                            vItrPos[ 1 ] = vItrTopRight[ 1 ];
                        } while( vItrTopRight[ 1 ] < vTopRight[ 1 ] &&
                                 vItrTopRight[ 1 ] != std::numeric_limits<coordinate_t>::max( ) );
                    }
                    else if( exists( rSparseCoords, vPos ) )
                        vRet.push_back( overlay_grid_t::indexOf( vPos, xEntry ) );
                }
                return vRet;
            } );
        size_t uiNumDone = 0;
        // actually process the overlays
        xIterator.process(
            rSparseCoords.THREADSAVE && rOverlays.THREADSAVE && vPoints.THREADSAVE
                ? std::thread::hardware_concurrency( )
                : 0,
            [ & ]( size_t uiI ) {
                typename points_t::Entry& xCurrPoints = vSplitPoints[ uiI - xOverlays.uiStartIndex ];

                // get bottom left position (compressed)
                pos_t vGridPos = rOverlays.posOf( uiI, xOverlays );
                if constexpr( DEPENDANT_DIMENSION )
                {
                    auto xPartialLock1 = rSparseCoords.xLockable.partialLock( );

                    assert( exists( rSparseCoords, vGridPos ) );
                    if( !exists( rSparseCoords, vGridPos ) )
                    {
#ifndef NDEBUG
                        xProg << Verbosity( 2 ) << "skipping nonexistant overlap\n";
#endif
                        // continue;
                        return;
                    }
                }

                pos_t vPosBottomLeftActual, vPosTopRightActual;
                std::array<std::vector<coordinate_t>, D> vPredecessors;
                {
                    auto xPartialLock1 = rSparseCoords.xLockable.partialLock( );
                    vPosBottomLeftActual = actualFromGridPos( rSparseCoords, vGridPos );
                    vPosTopRightActual = actualTopRightFromGridPos( rSparseCoords, vGridPos );

#ifndef NDEBUG
                    xProg << Verbosity( 2 ) << "overlay anchor is " << vGridPos << " actual pos is "
                          << vPosBottomLeftActual << "\n";
#endif

                    // collect direct predecessor overlays for each dimension
                    vPredecessors = getPredecessor( rOverlays, rSparseCoords, vGridPos, vPosBottomLeftActual,
                                                    vPosTopRightActual, xProg );

#ifndef NDEBUG
                    xProg << Verbosity( 2 );
                    if( xProg.active( ) )
                        for( size_t uiD = 0; uiD < D; uiD++ )
                            for( size_t uiJ = 0; uiJ < vPredecessors[ uiD ].size( ); uiJ++ )
                            {
                                auto xPartialLock2 = rPrefixSums.xLockable.partialLock( );
                                xProg << "predecessor " << uiJ << " dim " << uiD << " is "
                                      << vPredecessors[ uiD ][ uiJ ] << "\n";
                                assert( vPredecessors[ uiD ][ uiJ ] < uiI );

                                pos_t vGridPosPred = rOverlays.posOf( vPredecessors[ uiD ][ uiJ ], xOverlays );
                                rOverlays.vData[ vPredecessors[ uiD ][ uiJ ] ].stream(
                                    std::cout, vGridPosPred, rSparseCoords, rPrefixSums, *this, vPoints )
                                    << std::endl;
                            }
#endif
                }
#ifndef NDEBUG
                xProg << Verbosity( 1 ) << "generating overlay ouf of "
                      << xCurrPoints.uiEndIndex - xCurrPoints.uiStartIndex << " points now...\n";
#endif
                // generate the overlay
                overlay_t xOverlay;
                // generate, then copy over to array, so that we do not mess up any cache in the vector.
                xOverlay.generate( rOverlays, rSparseCoords, rPrefixSums, vPoints, xCurrPoints, vPredecessors,
                                   vPosBottomLeftActual, vPosTopRightActual, this, uiNumDone, uiNumTotal, xProg,
                                   xProfiler, xPool );

                rOverlays.vData[ uiI ] = xOverlay;
#ifndef NDEBUG
                xProfiler.step( "overlay gen loop" );
#endif

                std::unique_lock xGuard( xLock );
                xProg << Verbosity( 2 );
                if( xProg.active( ) )
                {
                    auto xPartialLock1 = rSparseCoords.xLockable.partialLock( );
                    auto xPartialLock2 = rPrefixSums.xLockable.partialLock( );
                    std::cout << rOverlays.vData[ uiI ] << std::endl;
                    rOverlays.vData[ uiI ].stream( std::cout, vGridPos, rSparseCoords, rPrefixSums, *this, vPoints )
                        << std::endl;
                }

                ++uiNumDone;
                if( xProg.printAgain( ) )
                    xProg << Verbosity( 0 ) << "computed " << uiNumDone << " out of " << uiNumTotal
                          << " overlays, thats " << 100 * uiNumDone / (double)uiNumTotal << "%.\n";
            } );
        xProg << Verbosity( 0 ) << "done.\n";
    }

    std::array<bool, D> hasPredecessor( const sparse_coord_t& /*rSparseCoords*/, pos_t vPos ) const
    {
        std::array<bool, D> vRet;

        for( size_t uiI = 0; uiI < D; uiI++ )
            vRet[ uiI ] = vPos[ uiI ] > 0;

        return vRet;
    }

#pragma GCC diagnostic push
// vPos not used with DEPENDANT_DIMENSION == false
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
    bool exists( const sparse_coord_t& rSparseCoords, pos_t vPos ) const
    {
        if constexpr( DEPENDANT_DIMENSION )
            return rSparseCoords.axisSize( sparse_coord_t::at( xSparseCoordsDependantDimension, vPos[ 0 ] ) ) >
                   vPos[ 1 ];
        return true;
    }
#pragma GCC diagnostic pop

    pos_t actualFromGridPos( const sparse_coord_t& rSparseCoords, pos_t vPos ) const
    {
        pos_t vRet = rSparseCoords.invSparse( vPos, vSparseCoords );

        // fix dependant dimension
        if constexpr( DEPENDANT_DIMENSION )
            vRet[ 1 ] = rSparseCoords.invReplace(
                vPos[ 1 ],
                sparse_coord_t::at( xSparseCoordsDependantDimension,
                                    std::min( vPos[ 0 ], rSparseCoords.axisSize( vSparseCoords[ 0 ] ) - 1 ) ) );

        return vRet;
    }

    pos_t actualTopRightFromGridPos( const sparse_coord_t& rSparseCoords, pos_t vPos ) const
    {
        for( size_t uiI = 0; uiI < D; uiI++ )
            ++vPos[ uiI ];

        pos_t vRet = actualFromGridPos( rSparseCoords, vPos );

        // fix dependant dimension
        if constexpr( DEPENDANT_DIMENSION )
            vRet[ 1 ] = rSparseCoords.invReplace(
                vPos[ 1 ], sparse_coord_t::at( xSparseCoordsDependantDimension, vPos[ 0 ] - 1 ) );

        return vRet;
    }

    pos_t overlayCoord( const sparse_coord_t& rSparseCoords, const pos_t& vPos ) const
    {
        pos_t vRet = rSparseCoords.sparse( vPos, vSparseCoords );

        // fix dependant dimension
        if constexpr( DEPENDANT_DIMENSION )
        {
            if( vRet[ 0 ] != std::numeric_limits<coordinate_t>::max( ) )
                vRet[ 1 ] = rSparseCoords.replace( vPos[ 1 ],
                                                   sparse_coord_t::at( xSparseCoordsDependantDimension, vRet[ 0 ] ) );
            else
                vRet[ 1 ] =
                    rSparseCoords.replace( vPos[ 1 ],
                                           sparse_coord_t::at( xSparseCoordsDependantDimension,
                                                               rSparseCoords.axisSize( vSparseCoords[ 0 ] ) - 1 ) );
        }

        return vRet;
    }

    struct OverlayInfo
    {
        pos_t vBottomLeft, vTopRight;
        pos_t vGridPos;
        coordinate_t uiIdx;
        std::array<std::vector<coordinate_t>, D> vPredIds;
        std::vector<pos_t> vvPoints{ };
    };

    std::vector<OverlayInfo> getOverlayInfo( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
                                             const points_t& vPoints ) const
    {
        std::vector<OverlayInfo> vRet;

        for( coordinate_t uiI = 0; uiI < rOverlays.sizeOf( xOverlays ); uiI++ )
        {
            vRet.emplace_back( );
            vRet.back( ).uiIdx = uiI;
            vRet.back( ).vGridPos = rOverlays.posOf( uiI, xOverlays );
            vRet.back( ).vBottomLeft = actualFromGridPos( rSparseCoords, vRet.back( ).vGridPos );
            vRet.back( ).vTopRight = actualTopRightFromGridPos( rSparseCoords, vRet.back( ).vGridPos );
            vRet.back( ).vPredIds =
                getPredecessor( rOverlays, rSparseCoords, vRet.back( ).vGridPos, vRet.back( ).vBottomLeft,
                                vRet.back( ).vTopRight, typename type_defs::progress_stream_t( 0 ) );
            vPoints.iterate( [ & ]( const point_t& xP ) { vRet.back( ).vvPoints.push_back( xP.vPos ); },
                             rOverlays.vData[ uiI ].xPoints );
        }

        return vRet;
    }

    std::vector<std::array<pos_t, 3>> getOverlayGrid( const overlay_grid_t& rOverlays,
                                                      const sparse_coord_t& rSparseCoords ) const
    {
        std::vector<std::array<pos_t, 3>> vRet;

        for( coordinate_t uiI = 0; uiI < rOverlays.sizeOf( xOverlays ); uiI++ )
        {
            vRet.emplace_back( );
            pos_t vGridPos = rOverlays.posOf( uiI, xOverlays );
            vRet.back( )[ 0 ] = vGridPos;
            vRet.back( )[ 1 ] = actualFromGridPos( rSparseCoords, vGridPos );
            vRet.back( )[ 2 ] = actualTopRightFromGridPos( rSparseCoords, vGridPos );
        }
        return vRet;
    }

    coordinate_t overlayIndex( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
                               const pos_t& vPos ) const
    {
        pos_t vPosOverlay = overlayCoord( rSparseCoords, vPos );
        return rOverlays.indexOf( vPosOverlay, xOverlays );
    }

    sps_t get( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
               const prefix_sum_grid_t& rPrefixSums, const pos_t& vPos, progress_stream_t& xProg
#if PROFILE_GET
               ,
               std::shared_ptr<Profiler> pProfiler = std::make_shared<Profiler>( )
#endif
    ) const
    {
#if PROFILE_GET
        pProfiler->step( "dataset_get" );
#endif
        auto vSparsePos = overlayCoord( rSparseCoords, vPos );
#if GET_PROG_PRINTS
        xProg << Verbosity( 2 );
        if( xProg.active( ) )
            xProg << "\t" << vPos << " -> " << vSparsePos << "; that's overlay "
                  << rOverlays.indexOf( vSparsePos, xOverlays ) << "\n";
#endif
        for( size_t uiI = 0; uiI < D; uiI++ )
            if( vSparsePos[ uiI ] == std::numeric_limits<coordinate_t>::max( ) )
                return sps_t{ };
        return rOverlays.get( vSparsePos, xOverlays )
            .get( rSparseCoords, rPrefixSums, vPos, actualFromGridPos( rSparseCoords, vSparsePos ), xProg
#if PROFILE_GET
                  ,
                  pProfiler
#endif
            );
    }

    friend std::ostream& operator<<( std::ostream& os, const Dataset& rDataset )
    {
        os << "<" << std::endl;
        os << "\tvSparseCoords: ";
        os << rDataset.vSparseCoords << std::endl;

        os << "\txOverlayGrid: ";
        os << rDataset.xOverlays << std::endl;
        os << ">";

        return os;
    }

    std::ostream& stream( std::ostream& os, const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
                          const prefix_sum_grid_t& rPrefixSums, const points_t& vPoints, const desc_t& vDesc ) const
    {
        os << "<" << std::endl;
        os << "\tvSparseCoords: ";
        for( size_t uiI = 0; uiI < D; uiI++ )
            vSparseCoords[ uiI ].stream( os, rSparseCoords ) << " ";
        os << std::endl;
        if constexpr( DEPENDANT_DIMENSION )
        {
            os << "\txSparseCoordsDependantDimension: ";
            xSparseCoordsDependantDimension.stream( os, rSparseCoords ) << " ";
            os << std::endl;
        }

        os << "\txOverlayGrid: ";
        xOverlays.stream( os, rOverlays, rSparseCoords, rPrefixSums, *this, vPoints, vDesc ) << std::endl;

        os << ">";

        return os;
    }
};


} // namespace sps