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
#include <random>
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

    pos_t uiSizeOverlays; // only relevant if UNIFORM_OVERLAY_GRID is true
    pos_t uiMaxCoords; // only relevant if UNIFORM_OVERLAY_GRID is true
    pos_t uiMinCoords; // only relevant if UNIFORM_OVERLAY_GRID is true

    static std::mt19937 xGen;

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
    typename points_t::template points_sort_func_t<points_it_t, PointsBinComperator> sort_points_bin =
        typename points_t::template points_sort_func_t<points_it_t, PointsBinComperator>( );

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
    typename points_t::template points_sort_func_t<points_it_t, PointsComperator> sort_points =
        typename points_t::template points_sort_func_t<points_it_t, PointsComperator>( );

    struct PointsComperatorDistinct
    {
        bool operator( )( const point_t& a, const point_t& b ) const
        {
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
                if( a.vPos[ uiI ] > b.vPos[ uiI ] )
                    return false;
                else if( a.vPos[ uiI ] < b.vPos[ uiI ] )
                    return true;
                // else if positions are equal continue
            }
            return false;
        }

        point_t min_value( ) const
        {
            point_t xRet{ };
            return xRet;
        };

        point_t max_value( ) const
        {
            point_t xRet{ };
            for( size_t uiI = 0; uiI < D; uiI++ )
                xRet.vPos[ uiI ] = std::numeric_limits<coordinate_t>::max( );
            return xRet;
        };
    };

    using points_sort_func_t = typename points_t::template points_sort_func_t<points_it_t, PointsComperatorDistinct>;

    static points_sort_func_t sort_points_distinct;

  public:
    Dataset( ) : vSparseCoords( ), xSparseCoordsDependantDimension( )
    {}

    typename sparse_coord_t::Entry
    makeSparseCoords( sparse_coord_t& rSparseCoords, points_t& vPoints, typename points_t::Entry xPoints, size_t uiDim,
                      size_t uiNumBlocks, progress_stream_t xProg,
                      coordinate_t uiFixedStart = std::numeric_limits<coordinate_t>::max( ),
                      coordinate_t uiFixedEnd = std::numeric_limits<coordinate_t>::max( ) ) const
    {
        size_t uiNumCoords = 0;
        coordinate_t uiLast = std::numeric_limits<coordinate_t>::max( );
        coordinate_t uiFirst = std::numeric_limits<coordinate_t>::max( );
        sort_points( vPoints.vData.begin( ) + xPoints.uiStartIndex, vPoints.vData.begin( ) + xPoints.uiEndIndex,
                     PointsComperator( uiDim ) );
        vPoints.iterate(
            [ & ]( const point_t& xPoint ) {
                coordinate_t uiCurr = xPoint.vPos[ uiDim ];
                if( uiFirst == std::numeric_limits<coordinate_t>::max( ) )
                    uiFirst = uiCurr;
                if( uiCurr != uiLast )
                {
                    uiLast = uiCurr;
                    ++uiNumCoords;
                }
            },
            xPoints );
        size_t uiBlockSize = 1 + ( uiNumCoords - 1 ) / uiNumBlocks;
        xProg << "generating " << uiNumBlocks << " overlays in dimension " << uiDim << " for " << uiNumCoords
              << " different coordinates.\n";
        auto xStart = DimIterator( vPoints.cbegin( xPoints ), vPoints.cend( xPoints ), uiDim, uiBlockSize );
        auto xEnd = DimIterator( vPoints.cend( xPoints ), vPoints.cend( xPoints ), uiDim, uiBlockSize );
        if( uiFixedStart == std::numeric_limits<coordinate_t>::max( ) )
        {
            uiFixedEnd = uiLast;
            uiFixedStart = uiFirst;
        }
        return rSparseCoords.template addStartEnd<true>( xStart, xEnd, uiFixedStart, uiFixedEnd );
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


    coordinate_t generateSparseCoords( overlay_grid_t& rOverlays,
                                       sparse_coord_t& rSparseCoords,
                                       points_t& vPoints,
                                       progress_stream_t xProg,
                                       std::function<coordinate_t( size_t )>
                                           fDo )
    {
        rSparseCoords.vData.reserve( rSparseCoords.vData.size( ) + ( 2 << 8 ) );

        std::vector<coordinate_t> vPrefixSumSize;
        size_t uiFrom = xOverlays.uiStartIndex;
        size_t uiTo = uiFrom + overlay_grid_t::sizeOf( xOverlays );
        {
            ThreadPool xPool( rOverlays.THREADSAVE && rSparseCoords.THREADSAVE && vPoints.THREADSAVE
                                  ? std::thread::hardware_concurrency( )
                                  : 0 );
            vPrefixSumSize.resize( std::max( xPool.numThreads( ), (size_t)1 ) );

            std::mutex xPrintMutex;


            for( size_t uiOverlayId = uiFrom; uiOverlayId < uiTo; uiOverlayId++ )
                xPool.enqueue(
                    [ & ]( size_t uiTid, size_t uiOverlayId ) {
                        auto xGuard = rSparseCoords.getCapacityGuard( );
                        vPrefixSumSize[ uiTid ] += fDo( uiOverlayId );

                        std::lock_guard<std::mutex> xPrintGuard( xPrintMutex );
                        if( xProg.printAgain( ) )
                            xProg << Verbosity( 0 ) << "processed " << uiOverlayId - uiFrom << " out of "
                                  << uiTo - uiFrom << " overlays, thats "
                                  << 100.0 * ( (double)( uiOverlayId - uiFrom ) / (double)( uiTo - uiFrom ) ) << "%.\n";
                    },
                    uiOverlayId );

        } // end of scope for xPool


        rSparseCoords.shrink_to_fit( );

        coordinate_t uiRet = 0;
        for( coordinate_t uiX : vPrefixSumSize )
            uiRet += uiX;
        return uiRet;
    }

    coordinate_t generateInternalSparseCoords( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords,
                                               points_t& vPoints, progress_stream_t xProg )
    {
        return generateSparseCoords( rOverlays, rSparseCoords, vPoints, xProg, [ & ]( size_t uiOverlayId ) {
            // get bottom left position (compressed)
            if constexpr( DEPENDANT_DIMENSION )
                if( !exists( rSparseCoords, rOverlays.posOf( uiOverlayId, xOverlays ) ) )
                {
#ifndef NDEBUG
                    xProg << Verbosity( 2 ) << "skipping nonexistant overlap\n";
#endif
                    // continue;
                    return (coordinate_t)0;
                }

            return rOverlays.vData[ uiOverlayId ].generateInternalSparseCoords( rSparseCoords, vPoints, xProg );
        } );
    }

    static coordinate_t maxSpaceForOverlayCoords( pos_t uiNumCoordsPerDim, pos_t uiNumOverlaysPerDim )
    {
        return maxSpaceForInternalCoords( uiNumCoordsPerDim, uiNumOverlaysPerDim );
        return ( 1 + maxSpaceForInternalCoords( uiNumCoordsPerDim, uiNumOverlaysPerDim ) ) * ( D - 1 );
    }


    coordinate_t generateOverlaySparseCoords( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords,
                                              points_t& vPoints, progress_stream_t xProgIn )
    {
        rSparseCoords.vData.reserve( rSparseCoords.vData.size( ) + ( 2 << 8 ) );

        auto xIterator =
            rOverlays.template genIterator<D>( xOverlays, successors, *this, rOverlays, rSparseCoords, xOverlays );


        size_t uiNumThreads = rSparseCoords.THREADSAVE && rOverlays.THREADSAVE && vPoints.THREADSAVE
                                  ? std::thread::hardware_concurrency( )
                                  : 0;
        std::vector<coordinate_t> vPrefixSumSize( std::max( (size_t)1, uiNumThreads ) );
        // actually process the overlays
        xIterator.process( uiNumThreads, xProgIn, [ & ]( size_t uiTid, size_t uiOverlayId ) {
            auto xGuard = rSparseCoords.getCapacityGuard( );
            progress_stream_t xProg = xProgIn;
            // get bottom left position (compressed)
            pos_t vGridPos = rOverlays.posOf( uiOverlayId, xOverlays );
            if constexpr( DEPENDANT_DIMENSION )
                if( !exists( rSparseCoords, vGridPos ) )
                {
#ifndef NDEBUG
                    xProg << Verbosity( 2 ) << "skipping nonexistant overlap\n";
#endif
                    // continue;
                    return;
                }

            pos_t vPosBottomLeftActual = actualFromGridPos( rSparseCoords, vGridPos );
            pos_t vPosTopRightActual = actualTopRightFromGridPos( rSparseCoords, vGridPos );

#ifndef NDEBUG
            xProg << Verbosity( 2 ) << "overlay anchor is " << vGridPos << " actual pos is " << vPosBottomLeftActual
                  << "\n";
#endif

            // collect direct predecessor overlays for each dimension
            std::array<std::vector<coordinate_t>, D> vPredecessors =
                getPredecessor( rOverlays, rSparseCoords, vGridPos, vPosBottomLeftActual, vPosTopRightActual, xProg );


            vPrefixSumSize[ uiTid ] += rOverlays.vData[ uiOverlayId ].generateOverlaySparseCoords(
                rOverlays, rSparseCoords, vPoints, vPredecessors, vPosBottomLeftActual, vPosTopRightActual, xProg );
        } );
        rSparseCoords.shrink_to_fit( );

        coordinate_t uiRet = 0;
        for( coordinate_t uiX : vPrefixSumSize )
            uiRet += uiX;
        return uiRet;
    }

    void generateInternalPrefixSums( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords,
                                     prefix_sum_grid_t& rPrefixSums, points_t& vPoints,
                                     coordinate_t uiSizeInternalPrefixSums, progress_stream_t xProg )
    {
        rPrefixSums.vData.reserve( rPrefixSums.vData.size( ) + uiSizeInternalPrefixSums );

        size_t uiFrom = xOverlays.uiStartIndex;
        size_t uiTo = uiFrom + overlay_grid_t::sizeOf( xOverlays );
        {
            ThreadPool xPool( rSparseCoords.THREADSAVE && rOverlays.THREADSAVE && vPoints.THREADSAVE ? std::thread::hardware_concurrency( ) : 0 );

            for( size_t uiOverlayId = uiFrom; uiOverlayId < uiTo; uiOverlayId++ )
                xPool.enqueue(
                    [ & ]( size_t uiTid, size_t uiOverlayId, progress_stream_t xProgCpy ) {
                        rOverlays.vData[ uiOverlayId ].generateInternalPrefixSums( rOverlays, rSparseCoords,
                                                                                   rPrefixSums, vPoints, xProgCpy );

                        if( uiTid == 0 && xProg.printAgain( ) )
                            xProg << Verbosity( 0 ) << "processed " << uiOverlayId - uiFrom << " out of "
                                  << uiTo - uiFrom << " overlays, thats "
                                  << 100.0 * ( (double)( uiOverlayId - uiFrom ) / (double)( uiTo - uiFrom ) ) << "%.\n";
                    },
                    uiOverlayId, xProg );
        } // end of scope for xPool
    }

    // return vector of indices for successor overlays
    static std::vector<size_t> successors( size_t uiIdx, size_t uiDim,
                                           const typename overlay_grid_t::template Entry<D>& xEntry, Dataset& rThis,
                                           overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords,
                                           typename overlay_grid_t::template Entry<D>& xOverlays )
    {
        if constexpr( UNIFORM_OVERLAY_GRID )
        {
            auto vPos = overlay_grid_t::posOf( uiIdx, xEntry );
            ++vPos[ uiDim ];
            std::vector<size_t> vRet;
            if( vPos[ uiDim ] < xEntry.vAxisSizes[ uiDim ] )
                vRet.push_back( overlay_grid_t::indexOf( vPos, xEntry ) );
            else if( uiIdx + 1 == overlay_grid_t::sizeOf( xEntry ) + xEntry.uiStartIndex )
                vRet.push_back( std::numeric_limits<size_t>::max( ) ); // poison
            return vRet;
        }
        else
        {
            auto xGuard = rSparseCoords.getCapacityGuard( );
            std::vector<size_t> vRet;
            auto vPos = overlay_grid_t::posOf( uiIdx, xEntry );

            bool bFits = true;
            for( size_t uiI = 0; uiI < D && bFits; uiI++ )
            {
                ++vPos[ uiI ];
                if( DEPENDANT_DIMENSION && uiI == 1 )
                {
                    if( vPos[ 1 ] < xEntry.vAxisSizes[ 1 ] && rThis.exists( rSparseCoords, vPos ) )
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
                if( DEPENDANT_DIMENSION && uiDim == 0 )
                {
                    --vPos[ uiDim ];

                    pos_t vBottomLeft = rThis.actualFromGridPos( rSparseCoords, vPos );
                    pos_t vTopRight = rThis.actualTopRightFromGridPos( rSparseCoords, vPos );

                    pos_t vItrPos = vBottomLeft;
                    vItrPos[ 0 ] = vTopRight[ 0 ];
                    pos_t vItrTopRight;
                    do
                    {
                        coordinate_t uiItrIndex = rThis.overlayIndex( rOverlays, rSparseCoords, vItrPos );
                        pos_t vItrGridPos = rOverlays.posOf( uiItrIndex, xOverlays );
                        vItrTopRight = rThis.actualTopRightFromGridPos( rSparseCoords, vItrGridPos );

                        if( vItrTopRight[ 1 ] <= vTopRight[ 1 ] && rThis.exists( rSparseCoords, vItrGridPos ) )
                            vRet.push_back( uiItrIndex );

                        // move upwards
                        vItrPos[ 1 ] = vItrTopRight[ 1 ];
                    } while( vItrTopRight[ 1 ] < vTopRight[ 1 ] &&
                             vItrTopRight[ 1 ] != std::numeric_limits<coordinate_t>::max( ) );
                }
                else if( rThis.exists( rSparseCoords, vPos ) )
                    vRet.push_back( overlay_grid_t::indexOf( vPos, xEntry ) );
            }
            return vRet;
        }
    }

    void generateOverlayPrefixSums( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords,
                                    prefix_sum_grid_t& rPrefixSums, points_t& vPoints,
                                    coordinate_t uiSizeOverlayPrefixSums, progress_stream_t xProgIn )
    {
        rPrefixSums.vData.reserve( rPrefixSums.vData.size( ) + uiSizeOverlayPrefixSums );

        auto xIterator =
            rOverlays.template genIterator<D>( xOverlays, successors, *this, rOverlays, rSparseCoords, xOverlays );
        // actually process the overlays
        xIterator.process(
            rSparseCoords.THREADSAVE && rOverlays.THREADSAVE && rPrefixSums.THREADSAVE
                ? std::thread::hardware_concurrency( )
                : 0,
            xProgIn,
            [ & ]( size_t /* uiTid */, size_t uiOverlayId ) {
                progress_stream_t xProg = xProgIn;
                pos_t vGridPos = rOverlays.posOf( uiOverlayId, xOverlays );
                if constexpr( DEPENDANT_DIMENSION )
                    assert( exists( rSparseCoords, vGridPos ) );

                pos_t vPosBottomLeftActual = actualFromGridPos( rSparseCoords, vGridPos );
                pos_t vPosTopRightActual = actualTopRightFromGridPos( rSparseCoords, vGridPos );

                // collect direct predecessor overlays for each dimension
                std::array<std::vector<coordinate_t>, D> vPredecessors = getPredecessor(
                    rOverlays, rSparseCoords, vGridPos, vPosBottomLeftActual, vPosTopRightActual, xProg );

                rOverlays.vData[ uiOverlayId ].generateOverlayPrefixSums(
                    rOverlays, rSparseCoords, rPrefixSums, vPoints, vPredecessors, vPosBottomLeftActual, this, xProg );
            } );
    }

    static void pickRandomDistinct( size_t uiN, size_t uiMax, std::function<void( size_t )> fDo )
    {
        std::map<size_t, size_t> xPointers;
        std::uniform_int_distribution<size_t> xDist( 0, uiMax - 1 );
        for( size_t uiFrom = 0; uiFrom < uiN; uiFrom++ )
        {
            if( xPointers.count( uiFrom ) == 0 )
                xPointers[ uiFrom ] = uiFrom;

            size_t uiTo = xDist( xGen );

            if( xPointers.count( uiTo ) == 0 )
                xPointers[ uiTo ] = uiTo;

            // swap uiFrom and uiTo //

            size_t uiTmp = xPointers[ uiFrom ];
            xPointers[ uiFrom ] = xPointers[ uiTo ];
            xPointers[ uiTo ] = uiTmp;
        }

#ifndef NDEBUG
        // sanity check: every (relevant) value in xPointers exists only once
        std::set<size_t> xCheck;
        for( size_t uiFrom = 0; uiFrom < uiN; uiFrom++ )
        {
            assert( xPointers.count( uiFrom ) == 1 );
            assert( xCheck.count( xPointers[ uiFrom ] ) == 0 );
            xCheck.insert( xPointers[ uiFrom ] );
        }
#endif

        for( size_t uiI = 0; uiI < uiN; uiI++ )
            fDo( xPointers[ uiI ] );
    }

    static coordinate_t sampleInterval( const points_t& vPoints, pos_t uiFrom, pos_t uiTo,
                                        typename points_t::Entry& xSortedPoints, size_t uiD, bool bStart )
    {
        coordinate_t uiPos = bStart ? uiFrom[ uiD ] : uiTo[ uiD ];
        while( ( bStart && uiPos < uiTo[ uiD ] ) || ( !bStart && uiPos > uiFrom[ uiD ] ) )
        {
            if( !bStart ) // doing this here to avoid unsigned integer underflow after zero
                --uiPos;
            const size_t uiD2 = uiD != 0 ? 0 : 1;
            bool bFound = false;
            vPoints.forEqualRange(
                // fBefore
                [ & ]( const point_t& rP ) {
                    // ignore all points before the picked coordinate
                    if( rP.vPos[ uiD ] < uiPos )
                        return true;
                    if( rP.vPos[ uiD ] > uiPos )
                        return false;

                    // ignore all points before and after the box in the second dimensions
                    if( rP.vPos[ uiD2 ] < uiFrom[ uiD2 ] )
                        return true;
                    if( rP.vPos[ uiD2 ] >= uiTo[ uiD2 ] )
                        return false;

                    // point should be considered
                    return false;
                },
                // fAfter
                [ & ]( const point_t& rP ) {
                    // once we have found a single point we can stop
                    if( bFound )
                        return true;

                    // once we reached a point not on the picked coordinate we can stop (due to sorting order)
                    if( rP.vPos[ uiD ] != uiPos )
                        return true;

                    // once we have reached a point past the box in the second coordinate we can stop (due to
                    // sorting order)
                    if( rP.vPos[ uiD2 ] >= uiTo[ uiD2 ] )
                        return true;

                    // point should be considered
                    return false;
                },
                // fDo
                [ & ]( const point_t& rP ) {
                    // check if considered point is fully inside
                    for( size_t uiD = 0; uiD < D; uiD++ )
                        if( rP.vPos[ uiD ] < uiFrom[ uiD ] || rP.vPos[ uiD ] >= uiTo[ uiD ] )
                            return; // point not inside -> abort
                    // point inside
                    bFound = true;
                }, //
                xSortedPoints );
            if( bFound )
                break; // found the first index
            if( bStart )
                ++uiPos;
        }
        return uiPos;
    }

    static pos_t sampleIntervalSize( const points_t& vPoints, pos_t uiFrom, pos_t uiTo,
                                     std::array<typename points_t::Entry, D> xSortedPoints )
    {
        pos_t uiNumDistinct;
        for( size_t uiD = 0; uiD < D; uiD++ )
        {
            coordinate_t uiEnd = sampleInterval( vPoints, uiFrom, uiTo, xSortedPoints[ uiD ], uiD, false );
            coordinate_t uiStart = sampleInterval( vPoints, uiFrom, uiTo, xSortedPoints[ uiD ], uiD, true );
            if( uiEnd < uiStart )
                uiNumDistinct[ uiD ] = 0;
            else
                uiNumDistinct[ uiD ] = uiEnd - uiStart;
        }

        return uiNumDistinct;
    }

    static coordinate_t sampleIntervalSize( const points_t& vPoints, pos_t uiFrom, pos_t uiTo,
                                            typename points_t::Entry& xSortedPoints, size_t uiD )
    {
        coordinate_t uiEnd = sampleInterval( vPoints, uiFrom, uiTo, xSortedPoints, uiD, false );
        coordinate_t uiStart = sampleInterval( vPoints, uiFrom, uiTo, xSortedPoints, uiD, true );
        if( uiEnd < uiStart )
            return 0;
        else
            return uiEnd - uiStart;
    }

    static coordinate_t sampleNumDistinct( const points_t& vPoints, pos_t uiFrom, pos_t uiTo,
                                           const uint64_t uiNumSamples, typename points_t::Entry& xSortedPoints,
                                           size_t uiD )
    {
        const size_t uiN = std::min( uiNumSamples, uiTo[ uiD ] - uiFrom[ uiD ] );
        if( uiN > 0 )
        {
            size_t uiNumFound = 0;
            pickRandomDistinct( uiN, uiTo[ uiD ] - uiFrom[ uiD ], [ & ]( size_t uiI ) {
                const coordinate_t uiPickedPos = uiI + uiFrom[ uiD ];
                const size_t uiD2 = uiD != 0 ? 0 : 1;
                bool bFound = false;
                // search in points for one point that is on the picked coordinate inside the box
                vPoints.forEqualRange(
                    // fBefore
                    [ & ]( const point_t& rP ) {
                        // ignore all points before the picked coordinate
                        if( rP.vPos[ uiD ] < uiPickedPos )
                            return true;
                        if( rP.vPos[ uiD ] > uiPickedPos )
                            return false;

                        // ignore all points before and after the box in the second dimensions
                        if( rP.vPos[ uiD2 ] < uiFrom[ uiD2 ] )
                            return true;
                        if( rP.vPos[ uiD2 ] >= uiTo[ uiD2 ] )
                            return false;

                        // point should be considered
                        return false;
                    },
                    // fAfter
                    [ & ]( const point_t& rP ) {
                        // once we have found a single point we can stop
                        if( bFound )
                            return true;

                        // once we reached a point not on the picked coordinate we can stop (due to sorting order)
                        if( rP.vPos[ uiD ] != uiPickedPos )
                            return true;

                        // once we have reached a point past the box in the second coordinate we can stop (due to
                        // sorting order)
                        if( rP.vPos[ uiD2 ] >= uiTo[ uiD2 ] )
                            return true;

                        // point should be considered
                        return false;
                    },
                    // fDo
                    [ & ]( const point_t& rP ) {
                        // check if considered point is fully inside
                        for( size_t uiD = 0; uiD < D; uiD++ )
                            if( rP.vPos[ uiD ] < uiFrom[ uiD ] || rP.vPos[ uiD ] >= uiTo[ uiD ] )
                                return; // point not inside -> abort
                        // point inside
                        bFound = true;
                    }, //
                    xSortedPoints );
                // there is at least one point inside
                if( bFound )
                    ++uiNumFound;
            } );
            // scale up to actual number of coordinates within this box
            return ( ( uiTo[ uiD ] - uiFrom[ uiD ] ) * uiNumFound ) / uiN;
        }
        else
            return 0;
    }

    static pos_t sampleNumDistinct( const points_t& vPoints, pos_t uiFrom, pos_t uiTo, const uint64_t uiNumSamples,
                                    std::array<typename points_t::Entry, D> xSortedPoints )
    {
        pos_t uiNumDistinct;
        for( size_t uiD = 0; uiD < D; uiD++ )
            uiNumDistinct[ uiD ] = sampleNumDistinct( vPoints, uiFrom, uiTo, uiNumSamples, xSortedPoints[ uiD ], uiD );

        return uiNumDistinct;
    }


#pragma GCC diagnostic push
// uiCoordinateSizes not used with UNIFORM_OVERLAY_GRID == false
// vvCoords not used with UNIFORM_OVERLAY_GRID == true
// vvDepDim not used with DEPENDANT_DIMENSION = false
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
    static std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>
    estimateOverlay( const points_t& vPoints, std::array<typename points_t::Entry, D> xSortedPoints,
                     std::array<std::vector<coordinate_t>, D> vvCoords, std::vector<std::vector<coordinate_t>> vvDepDim,
                     size_t, pos_t uiCoordinateSizes, pos_t uiNumOverlays, uint64_t, pos_t uiP, pos_t uiMinPos,
                     const uint64_t uiNumPointSamples )
    {
        uint64_t uiNumPSInternal = 1;
        uint64_t uiNumPSOverlay = 0;
        uint64_t uiNumLookUpTablesInternal = 0;
        uint64_t uiNumLookUpTablesOverlay = 0;
        pos_t uiFromGlob, uiToGlob;
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            if constexpr( UNIFORM_OVERLAY_GRID )
            {
                uiFromGlob[ uiI ] =
                    uiP[ uiI ] * ( 1 + ( uiCoordinateSizes[ uiI ] - 1 ) / uiNumOverlays[ uiI ] ) + uiMinPos[ uiI ];
                uiToGlob[ uiI ] = ( uiP[ uiI ] + 1 ) * ( 1 + ( uiCoordinateSizes[ uiI ] - 1 ) / uiNumOverlays[ uiI ] ) +
                                  uiMinPos[ uiI ];
            }
            else
            {
                if( DEPENDANT_DIMENSION && uiI == 1 )
                {
                    if( uiP[ 0 ] >= vvDepDim.size( ) )
                    {
                        uiFromGlob[ 1 ] = 0;
                        uiToGlob[ 1 ] = 0;
                    }
                    else
                    {
                        uiFromGlob[ 1 ] =
                            vvDepDim[ uiP[ 0 ] ][ ( uiP[ 1 ] * vvDepDim[ uiP[ 0 ] ].size( ) ) / uiNumOverlays[ 1 ] ];
                        size_t uiEnd = ( ( 1 + uiP[ 1 ] ) * vvDepDim[ uiP[ 0 ] ].size( ) ) / uiNumOverlays[ 1 ];
                        if( uiEnd < vvDepDim[ uiP[ 0 ] ].size( ) )
                            uiToGlob[ 1 ] = vvDepDim[ uiP[ 0 ] ][ uiEnd ];
                        else
                            uiToGlob[ 1 ] = vvDepDim[ uiP[ 0 ] ].back( ) + 1;
                    }
                }
                else
                {
                    uiFromGlob[ uiI ] =
                        vvCoords[ uiI ][ ( uiP[ uiI ] * vvCoords[ uiI ].size( ) ) / uiNumOverlays[ uiI ] ];
                    size_t uiEnd = ( ( 1 + uiP[ uiI ] ) * vvCoords[ uiI ].size( ) ) / uiNumOverlays[ uiI ];
                    if( uiEnd < vvCoords[ uiI ].size( ) )
                        uiToGlob[ uiI ] = vvCoords[ uiI ][ uiEnd ];
                    else
                        uiToGlob[ uiI ] = vvCoords[ uiI ].back( ) + 1;
                }
            }
        }
        pos_t uiNumDistinctInOverlay =
            sampleNumDistinct( vPoints, uiFromGlob, uiToGlob, uiNumPointSamples, xSortedPoints );
        pos_t uiSampledIntervalSize = sampleIntervalSize( vPoints, uiFromGlob, uiToGlob, xSortedPoints );


        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            if( uiP[ uiI ] > 0 )
            {
                uint64_t uiNumPSOverlayCurr = 1;
                for( size_t uiJ = 0; uiJ < D; uiJ++ )
                    if( uiJ != uiI )
                    {
                        pos_t uiFrom, uiTo;
                        for( size_t uiK = 0; uiK < D; uiK++ )
                        {
                            if( uiK == uiI )
                                uiTo[ uiK ] = uiFromGlob[ uiK ];
                            else
                                uiTo[ uiK ] = uiToGlob[ uiK ];
                            if( uiK != uiJ )
                                uiFrom[ uiK ] = uiMinPos[ uiK ];
                            else
                                uiFrom[ uiK ] = uiFromGlob[ uiK ];
                        }
                        uiNumPSOverlayCurr *=
                            ( uiP[ uiJ ] > 0 ? 1 : 0 ) +
                            sampleNumDistinct( vPoints, uiFrom, uiTo, uiNumPointSamples, xSortedPoints[ uiJ ], uiJ );

                        uiNumLookUpTablesOverlay +=
                            ( uiP[ uiJ ] > 0 ? 1 : 0 ) +
                            sampleIntervalSize( vPoints, uiFrom, uiTo, xSortedPoints[ uiJ ], uiJ );
                    }
                uiNumPSOverlay += uiNumPSOverlayCurr;
            }

            // prefix sums in overlay
            // use the coupons collectors problem to adjust for several points lying in the same row or column of
            // the overlay the size of the prefix sum matrix is determined merely by the number of non-empty columns
            // and rows. so if multiple points lie in the same spot they do not need both be considered. translated
            // into the coupons collectors problem: we have num_sparse_coords / num_overlays_in_axis many different
            // coupons in total and we draw uiNumAvgPointsInOverlay many times. returned we get how many different
            // coupons we have drawn, which correspons to the number of rows that actually have at least one point.
            uiNumPSInternal *= uiNumDistinctInOverlay[ uiI ];

            uiNumLookUpTablesInternal += uiSampledIntervalSize[ uiI ];
        }
        return std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>(
            uiNumPSOverlay, uiNumPSInternal, uiNumLookUpTablesOverlay, uiNumLookUpTablesInternal );
    }
#pragma GCC diagnostic pop

    static std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>
    estimateDataStructureElements( points_t& vPoints, std::array<typename points_t::Entry, D> xSortedPoints,
                                   std::array<std::vector<coordinate_t>, D> vvCoords, pos_t uiNumOverlays,
                                   pos_t uiCoordinateSizes, size_t uiNumPoints, pos_t uiMinPos,
                                   const uint64_t uiNumOverlaySamples, const uint64_t uiNumPointSamples )
    {
        uint64_t uiNumOverlaysTotal = 1;
        for( size_t uiI = 0; uiI < D; uiI++ )
            uiNumOverlaysTotal *= uiNumOverlays[ uiI ];
        std::vector<std::vector<coordinate_t>> vvDepDim;
        if constexpr( !UNIFORM_OVERLAY_GRID && DEPENDANT_DIMENSION )
        {
            typename points_t::Entry xTmp = vPoints.copyEntry( xSortedPoints[ 0 ] );

            size_t uiCurrOverlay = 0;
            vPoints.iterate(
                [ & ]( point_t& rP ) {
                    size_t uiEnd = ( ( 1 + uiCurrOverlay ) * vvCoords[ 0 ].size( ) ) / uiNumOverlays[ 0 ];
                    if( uiEnd < vvCoords[ 0 ].size( ) && vvCoords[ 0 ][ uiEnd ] >= rP.vPos[ 0 ] )
                        ++uiCurrOverlay;
#ifndef NBEBUG
                    uiEnd = ( ( 1 + uiCurrOverlay ) * vvCoords[ 0 ].size( ) ) / uiNumOverlays[ 0 ];
                    assert( uiEnd >= vvCoords[ 0 ].size( ) || rP.vPos[ 0 ] < vvCoords[ 0 ][ uiEnd ] );
#endif
                    rP.vPos[ 0 ] = uiCurrOverlay;
                },
                xTmp );

            vPoints.sortByDim( 0, 1, xTmp );

            coordinate_t uiLastX = std::numeric_limits<coordinate_t>::max( );
            vPoints.iterate(
                [ & ]( const point_t& rP ) {
                    if( vvDepDim.size( ) == 0 || ( vvDepDim.back( ).size( ) > 0 && uiLastX != rP.vPos[ 0 ] ) )
                        vvDepDim.emplace_back( );
                    uiLastX = rP.vPos[ 0 ];
                    if( vvDepDim.back( ).size( ) == 0 || rP.vPos[ 1 ] != vvDepDim.back( ).back( ) )
                        vvDepDim.back( ).push_back( rP.vPos[ 1 ] );
                },
                xTmp );

            vPoints.popEntry( xTmp );
        }

        std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> tTotal{ };
        {
            ThreadPool xPool( vPoints.THREADSAVE ? std::thread::hardware_concurrency( ) : 0 );

            std::mutex xResMutex;

            for( size_t uiI = 0; uiI < uiNumOverlaySamples; uiI++ )
                xPool.enqueue( [ & ]( size_t ) {
                    pos_t uiP{ };
                    for( size_t uiI = 0; uiI < D; uiI++ )
                    {
                        std::uniform_int_distribution<uint64_t> xDis( 0, uiNumOverlays[ uiI ] - 1 );
                        uiP[ uiI ] = xDis( xGen );
                    }
                    auto tCurr =
                        estimateOverlay( vPoints, xSortedPoints, vvCoords, vvDepDim, uiNumPoints, uiCoordinateSizes,
                                         uiNumOverlays, uiNumOverlaysTotal, uiP, uiMinPos, uiNumPointSamples );


                    std::lock_guard<std::mutex> xGuard( xResMutex );
                    std::get<0>( tTotal ) += std::get<0>( tCurr );
                    std::get<1>( tTotal ) += std::get<1>( tCurr );
                    std::get<2>( tTotal ) += std::get<2>( tCurr );
                    std::get<3>( tTotal ) += std::get<3>( tCurr );
                } );
        } // scope for xPool


        uint64_t uiNumBlockLookup = 0;
        if constexpr( !UNIFORM_OVERLAY_GRID )
        {
            for( size_t uiI = 0; uiI < D; uiI++ )
                uiNumBlockLookup += uiCoordinateSizes[ uiI ];
            if constexpr( DEPENDANT_DIMENSION )
                uiNumBlockLookup += ( uiCoordinateSizes[ 1 ] ) * ( uiNumOverlays[ 0 ] - 1 );
        }


        uint64_t uiNumPSTotal =
            ( ( uiNumOverlaysTotal * ( std::get<0>( tTotal ) + std::get<1>( tTotal ) ) ) / uiNumOverlaySamples ) *
            sizeof( sps_t );

        uint64_t uiNumLookupTotal =
            ( ( uiNumOverlaysTotal * ( std::get<2>( tTotal ) + std::get<3>( tTotal ) ) ) / uiNumOverlaySamples +
              uiNumBlockLookup ) *
            sizeof( coordinate_t );

        uint64_t uiSizeOverlaysOverhead = uiNumOverlaysTotal * sizeof( overlay_t );

        // full
        uint64_t uiSizeFull = uiNumPSTotal + uiNumLookupTotal + uiSizeOverlaysOverhead + sizeof( Dataset );

        return std::make_tuple( ( uiNumOverlaysTotal * std::get<0>( tTotal ) ) / uiNumOverlaySamples,
                                ( uiNumOverlaysTotal * std::get<1>( tTotal ) ) / uiNumOverlaySamples,
                                ( uiNumOverlaysTotal * std::get<2>( tTotal ) ) / uiNumOverlaySamples,
                                ( uiNumOverlaysTotal * std::get<3>( tTotal ) ) / uiNumOverlaySamples,
                                uiNumOverlaysTotal,
                                uiNumBlockLookup,
                                uiSizeFull );
    }

    static pos_t toNumbers( std::array<double, D> vRatios, double fFactor )
    {
        pos_t uiRet;
        for( size_t uiI = 0; uiI < D; uiI++ )
            uiRet[ uiI ] = std::max( (coordinate_t)1, (coordinate_t)( vRatios[ uiI ] * fFactor ) );

        return uiRet;
    }

    static std::array<double, D> toNumRatios( pos_t uiCoordinateSizes )
    {
        std::array<double, D> vNumRatios;
        double fTotal = 0;

        for( size_t uiI = 0; uiI < D; uiI++ )
            fTotal += (double)uiCoordinateSizes[ uiI ];

        for( size_t uiI = 0; uiI < D; uiI++ )
            vNumRatios[ uiI ] = (double)uiCoordinateSizes[ uiI ] / fTotal;

        return vNumRatios;
    }

    static double toFactor( pos_t uiCoordinateSizes, uint64_t uiNumOverlays )
    {
        auto vRatios = toNumRatios( uiCoordinateSizes );
        double fStartPos = 0;
        double fEndPos = 1000.0;
        while( fEndPos - fStartPos > 0.01 )
        {
            double fCenter = ( fStartPos + fEndPos ) / 2;
            auto vNums = toNumbers( vRatios, fCenter );

            uint64_t uiNumCurr = 1;
            for( size_t uiI = 0; uiI < vNums.size( ); uiI++ )
                uiNumCurr *= vNums[ uiI ];

            if( uiNumCurr == uiNumOverlays )
                return fCenter;
            else if( uiNumCurr > uiNumOverlays )
                fEndPos = fCenter;
            else
                fStartPos = fCenter;
        }
        return ( fStartPos + fEndPos ) / 2;
    }

    static double toFactor( points_t& vPoints, const typename points_t::Entry xPoints, uint64_t uiNumOverlays )
    {
        pos_t uiCoordinateSizes = generateCoordSizes( vPoints, xPoints )[ 0 ];
        return toFactor( uiCoordinateSizes, uiNumOverlays );
    }

    static std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>>
    estimateDataStructureElements( points_t& vPoints, const typename points_t::Entry xPoints, std::vector<double> vFac,
                                   const uint64_t uiNumOverlaySamples, const uint64_t uiNumPointSamples )
    {
        std::array<typename points_t::Entry, D> xSortedPoints;
        std::array<std::vector<coordinate_t>, D> vvCoords;
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            xSortedPoints[ uiI ] = vPoints.copyEntry( xPoints );

            vPoints.sortByDim( uiI, uiI != 0 ? 0 : 1, xSortedPoints[ uiI ] );

            if constexpr( !UNIFORM_OVERLAY_GRID )
                vPoints.iterate(
                    [ & ]( const point_t& rP ) {
                        if( vvCoords[ uiI ].size( ) == 0 || rP.vPos[ uiI ] != vvCoords[ uiI ].back( ) )
                            vvCoords[ uiI ].push_back( rP.vPos[ uiI ] );
                    },
                    xSortedPoints[ uiI ] );
        }

        auto xA = generateCoordSizes( vPoints, xPoints );
        pos_t uiCoordinateSizes = xA[ 0 ];
        pos_t uiMinPos = xA[ 2 ];
        std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>> vRet;
        for( float fFac : vFac )
            vRet.push_back( estimateDataStructureElements( vPoints, xSortedPoints, vvCoords,
                                                           toNumbers( toNumRatios( uiCoordinateSizes ), fFac ),
                                                           uiCoordinateSizes, xPoints.uiEndIndex - xPoints.uiStartIndex,
                                                           uiMinPos, uiNumOverlaySamples, uiNumPointSamples ) );

        // remove additional entries again
        for( size_t uiI = 0; uiI < D; uiI++ )
            vPoints.popEntry( xSortedPoints[ D - uiI - 1 ] );

        return vRet;
    }

    static coordinate_t toAmount( pos_t vNums )
    {
        coordinate_t uiRet = 1;
        for( size_t uiI = 0; uiI < D; uiI += 1 )
            uiRet *= vNums[ uiI ];

        return uiRet;
    }


    static double pickOverlayFactor( points_t& vPoints, const typename points_t::Entry xPoints, pos_t uiCoordinateSizes,
                                     size_t uiNumPoints, pos_t uiMinPos, const uint64_t uiNumOverlaySamples,
                                     const uint64_t uiNumPointSamples, progress_stream_t xProg )
    {
        std::array<typename points_t::Entry, D> xSortedPoints;
        std::array<std::vector<coordinate_t>, D> vvCoords;
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            xSortedPoints[ uiI ] = vPoints.copyEntry( xPoints );

            vPoints.sortByDim( uiI, uiI != 0 ? 0 : 1, xSortedPoints[ uiI ] );

            if constexpr( !UNIFORM_OVERLAY_GRID )
                vPoints.iterate(
                    [ & ]( const point_t& rP ) {
                        if( vvCoords[ uiI ].size( ) == 0 || rP.vPos[ uiI ] != vvCoords[ uiI ].back( ) )
                            vvCoords[ uiI ].push_back( rP.vPos[ uiI ] );
                    },
                    xSortedPoints[ uiI ] );
        }


        std::array<double, D> vNumRatios = toNumRatios( uiCoordinateSizes );

        double fEnd = 100;
        double fStart = 0;


        // make sure to place end past the minimum
        while( true )
        {
            xProg << Verbosity( 0 ) << "searching factors [" << fStart << ", " << fEnd << ")\n";
            uint64_t uiEnd = std::get<6>( estimateDataStructureElements(
                vPoints, xSortedPoints, vvCoords, toNumbers( vNumRatios, fEnd ), uiCoordinateSizes, uiNumPoints,
                uiMinPos, uiNumOverlaySamples, uiNumPointSamples ) );
            uint64_t uiCenter = std::get<6>( estimateDataStructureElements(
                vPoints, xSortedPoints, vvCoords, toNumbers( vNumRatios, fEnd / 2 ), uiCoordinateSizes, uiNumPoints,
                uiMinPos, uiNumOverlaySamples, uiNumPointSamples ) );
            if( uiEnd > uiCenter + 10000 )
                break;
            fEnd *= 2;
        }

        // fund the minumum between 0 and fEnd
        const double fSampleSteps = 10;
        uint64_t uiMin = std::numeric_limits<uint64_t>::max( );
        while( toAmount( toNumbers( vNumRatios, fStart ) ) != toAmount( toNumbers( vNumRatios, fEnd ) ) )
        {
            xProg << Verbosity( 0 ) << "searching factors [" << fStart << ", " << fEnd << ")\n";
            double fMin = 0;
            uiMin = std::numeric_limits<uint64_t>::max( );
            for( double fPos = fStart; fPos <= fEnd; fPos += ( fEnd - fStart ) / fSampleSteps )
            {
                uint64_t uiCurr = std::get<6>( estimateDataStructureElements(
                    vPoints, xSortedPoints, vvCoords, toNumbers( vNumRatios, fPos ), uiCoordinateSizes, uiNumPoints,
                    uiMinPos, uiNumOverlaySamples, uiNumPointSamples ) );
                if( uiCurr < uiMin )
                {
                    uiMin = uiCurr;
                    fMin = fPos;
                }
            }
            fStart = std::max( 0.0, fMin - ( fEnd - fStart ) / fSampleSteps );
            fEnd = fMin + ( fEnd - fStart ) / fSampleSteps;
        }


        // remove additional entries again
        for( size_t uiI = 0; uiI < D; uiI++ )
            vPoints.popEntry( xSortedPoints[ D - uiI - 1 ] );

        xProg << Verbosity( 0 ) << "picked factor " << fStart << " -> expected index size is "
              << uiMin / (double)std::pow( 10, 9 ) << "Gb\n";

        return fStart;
    }

    static coordinate_t pickNumOverlays( points_t& vPoints, const typename points_t::Entry xPoints,
                                         pos_t uiCoordinateSizes, size_t uiNumPoints, pos_t uiMinPos,
                                         const uint64_t uiNumOverlaySamples, const uint64_t uiNumPointSamples,
                                         progress_stream_t xProg )
    {
        auto vNums = toNumbers( toNumRatios( uiCoordinateSizes ),
                                pickOverlayFactor( vPoints, xPoints, uiCoordinateSizes, uiNumPoints, uiMinPos,
                                                   uiNumOverlaySamples, uiNumPointSamples, xProg ) );

        return toAmount( vNums );
    }

    static coordinate_t pickNumOverlays( points_t& vPoints, const typename points_t::Entry xPoints,
                                         const uint64_t uiNumOverlaySamples, const uint64_t uiNumPointSamples,
                                         progress_stream_t xProg )
    {
        auto xA = generateCoordSizes( vPoints, xPoints );
        pos_t uiCoordinateSizes = xA[ 0 ];
        pos_t uiMinPos = xA[ 2 ];
        return pickNumOverlays( vPoints, xPoints, uiCoordinateSizes, xPoints.uiEndIndex - xPoints.uiStartIndex,
                                uiMinPos, uiNumOverlaySamples, uiNumPointSamples, xProg );
    }

    static pos_t pickOverlayNumbers( points_t& vPoints, const typename points_t::Entry xPoints, pos_t uiCoordinateSizes,
                                     size_t uiNumPoints, pos_t uiMinPos, double fFac,
                                     const uint64_t uiNumOverlaySamples, const uint64_t uiNumPointSamples,
                                     progress_stream_t xProg )
    {
        std::array<double, D> vNumRatios = toNumRatios( uiCoordinateSizes );
        if( fFac >= 0 )
            return toNumbers( vNumRatios, fFac );
        return toNumbers( vNumRatios,
                          pickOverlayFactor( vPoints, xPoints, uiCoordinateSizes, uiNumPoints, uiMinPos,
                                             uiNumOverlaySamples, uiNumPointSamples, xProg ) );
    }

    static std::array<pos_t, 3> generateCoordSizes( points_t& vPoints, const typename points_t::Entry xPoints )
    {
        pos_t uiMaxCoords;
        pos_t uiMinCoords;
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            uiMinCoords[ uiI ] = std::numeric_limits<coordinate_t>::max( );
            uiMaxCoords[ uiI ] = 0;
        }

        vPoints.iterate(
            [ & ]( const point_t& xPoint ) {
                for( size_t uiI = 0; uiI < D; uiI++ )
                {
                    uiMinCoords[ uiI ] = std::min( uiMinCoords[ uiI ], xPoint.vPos[ uiI ] );
                    uiMaxCoords[ uiI ] = std::max( uiMaxCoords[ uiI ], xPoint.vPos[ uiI ] + 1 );
                }
            },
            xPoints );

        pos_t uiCoordinateSizes;
        for( size_t uiI = 0; uiI < D; uiI++ )
            uiCoordinateSizes[ uiI ] = uiMaxCoords[ uiI ] - uiMinCoords[ uiI ];
        return std::array<pos_t, 3>{ uiCoordinateSizes, uiMaxCoords, uiMinCoords };
    }

#pragma GCC diagnostic push
// vPos not used with DEPENDANT_DIMENSION == false
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
    void generateOverlayCoords( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords, points_t& vPoints,
                                typename points_t::Entry xPoints, double fFac, progress_stream_t xProg,
                                const uint64_t uiNumOverlaySamples, const uint64_t uiNumPointSamples )
    {
        auto xA = generateCoordSizes( vPoints, xPoints );
        pos_t uiCoordinateSizes = xA[ 0 ];
        uiMaxCoords = xA[ 1 ];
        uiMinCoords = xA[ 2 ];

        xProg << Verbosity( 0 ) << "picking number of overlay boxes.\n";
        pos_t uiNumOverlaysPerDim =
            pickOverlayNumbers( vPoints, xPoints, uiCoordinateSizes, xPoints.uiEndIndex - xPoints.uiStartIndex,
                                uiMinCoords, fFac, uiNumOverlaySamples, uiNumPointSamples, xProg );

        xProg << Verbosity( 0 ) << "generating overlay grid.\n";
        if constexpr( UNIFORM_OVERLAY_GRID )
        {
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
                // 'round up' size to make sure that no points are past the last overlay
                uiSizeOverlays[ uiI ] = 1 + ( uiCoordinateSizes[ uiI ] - 1 ) / uiNumOverlaysPerDim[ uiI ];
                xProg << "generating " << uiNumOverlaysPerDim[ uiI ] << " overlays in dimension " << uiI << "\n";
            }

            xOverlays = rOverlays.template add<D, true, true>( uiNumOverlaysPerDim );
        }
        else
        {
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

                            auto xCurr =
                                makeSparseCoords( rSparseCoords, vPoints, xCurrPoints, 1, uiNumOverlaysPerDim[ uiI ],
                                                  xProg, uiAbsoluteStart, uiAbsoluteEnd );
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
                vSparseCoords[ uiI ] =
                    makeSparseCoords( rSparseCoords, vPoints, xPoints, uiI, uiNumOverlaysPerDim[ uiI ], xProg );
            }
            xProg << Verbosity( 2 );
            if( xProg.active( ) )
            {
                vSparseCoords[ 0 ].stream( std::cout << "sparse coords 0: ", rSparseCoords ) << std::endl;
                if constexpr( DEPENDANT_DIMENSION )
                    xSparseCoordsDependantDimension.stream( std::cout << "sparse coords 1: ", rSparseCoords )
                        << std::endl;
                else
                    vSparseCoords[ 1 ].stream( std::cout << "sparse coords 1: ", rSparseCoords ) << std::endl;
                for( size_t uiI = 2; uiI < D; uiI++ )
                    vSparseCoords[ uiI ].stream( std::cout << "sparse coords " << uiI << ": ", rSparseCoords )
                        << std::endl;
            }


            // generate overlay grid
            xOverlays = rOverlays.template add<D, true, true>( rSparseCoords.axisSizes( vSparseCoords ) );
        }

        // sort points so that they match the overlay grid order
        xProg << Verbosity( 0 ) << "sorting points into overlays.\n";
        /* Weirdly stxxl::sort seems to be parallel eventhough the access to stxxl::vector is required to be sequential.
         * This does not seem to be a problem for the vector that is currently sorted.
         * However, my PointsBinComperator accesses another stxxl::vector and therefore has to be guarded with a mutex.
         */
        std::mutex xLock;
        sort_points_bin( vPoints.vData.begin( ) + xPoints.uiStartIndex, vPoints.vData.begin( ) + xPoints.uiEndIndex,
                         PointsBinComperator( *this, rOverlays, rSparseCoords, xLock ) );

        // generate all overlays
        coordinate_t uiNumTotal = rOverlays.sizeOf( xOverlays ); // @todo count considering non existant entries

        // @todo vector does not need to be stored here, running variable is enough
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

            rOverlays.vData[ uiI + xOverlays.uiStartIndex ].xPoints = vSplitPoints[ uiI ];
        }
        assert( vSplitPoints[ uiNumTotal - 1 ].uiEndIndex == xPoints.uiEndIndex );
    }
#pragma GCC diagnostic pop

    Dataset( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords, prefix_sum_grid_t& rPrefixSums,
             points_t& vPoints, typename points_t::Entry xPoints, double fFac, progress_stream_t xProg,
             const uint64_t uiNumOverlaySamples, const uint64_t uiNumPointSamples )
        : vSparseCoords( ), xSparseCoordsDependantDimension( )
    {
        if( xPoints.uiEndIndex == xPoints.uiStartIndex )
            return;
        // generate the overall sparse coordinates
        generateOverlayCoords( rOverlays, rSparseCoords, vPoints, xPoints, fFac, xProg, uiNumOverlaySamples,
                               uiNumPointSamples );

        xProg << Verbosity( 0 ) << "generating internal sparse coords.\n";
        coordinate_t uiSizeInternalPrefixSums =
            generateInternalSparseCoords( rOverlays, rSparseCoords, vPoints, xProg );

        xProg << Verbosity( 0 ) << "generating overlay sparse coords.\n";
        coordinate_t uiSizeOverlayPrefixSums = generateOverlaySparseCoords( rOverlays, rSparseCoords, vPoints, xProg );

        xProg << Verbosity( 0 ) << "generating internal prefix sums.\n";
        generateInternalPrefixSums( rOverlays, rSparseCoords, rPrefixSums, vPoints, uiSizeInternalPrefixSums, xProg );

        xProg << Verbosity( 0 ) << "generating overlay prefix sums.\n";
        generateOverlayPrefixSums( rOverlays, rSparseCoords, rPrefixSums, vPoints, uiSizeOverlayPrefixSums, xProg );

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
        if constexpr( UNIFORM_OVERLAY_GRID )
        {
            pos_t vRet;
            for( size_t uiI = 0; uiI < D; uiI++ )
                vRet[ uiI ] = vPos[ uiI ] * uiSizeOverlays[ uiI ] + uiMinCoords[ uiI ];
            return vRet;
        }
        else
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
        if constexpr( UNIFORM_OVERLAY_GRID )
        {
            pos_t vRet;
            for( size_t uiI = 0; uiI < D; uiI++ )
                if( vPos[ uiI ] < uiMinCoords[ uiI ] )
                    vRet[ uiI ] = std::numeric_limits<coordinate_t>::max( );
                else
                    vRet[ uiI ] = std::min( xOverlays.vAxisSizes[ uiI ] - 1,
                                            ( vPos[ uiI ] - uiMinCoords[ uiI ] ) / uiSizeOverlays[ uiI ] );
            return vRet;
        }
        else
        {
            pos_t vRet = rSparseCoords.sparse( vPos, vSparseCoords );

            // fix dependant dimension
            if constexpr( DEPENDANT_DIMENSION )
            {
                if( vRet[ 0 ] != std::numeric_limits<coordinate_t>::max( ) )
                    vRet[ 1 ] = rSparseCoords.replace(
                        vPos[ 1 ], sparse_coord_t::at( xSparseCoordsDependantDimension, vRet[ 0 ] ) );
                else
                    vRet[ 1 ] =
                        rSparseCoords.replace( vPos[ 1 ],
                                               sparse_coord_t::at( xSparseCoordsDependantDimension,
                                                                   rSparseCoords.axisSize( vSparseCoords[ 0 ] ) - 1 ) );
            }

            return vRet;
        }
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

    val_t get( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
               const prefix_sum_grid_t& rPrefixSums, const pos_t& vPos, size_t uiCornerIdx
#if GET_PROG_PRINTS
               ,
               progress_stream_t& xProg
#endif
    ) const
    {
        auto vSparsePos = overlayCoord( rSparseCoords, vPos );
#if GET_PROG_PRINTS
        xProg << Verbosity( 2 );
        if( xProg.active( ) )
            xProg << "\t" << vPos << " -> " << vSparsePos << "; that's overlay "
                  << rOverlays.indexOf( vSparsePos, xOverlays ) << "\n";
#endif
        for( size_t uiI = 0; uiI < D; uiI++ )
            if( vSparsePos[ uiI ] == std::numeric_limits<coordinate_t>::max( ) )
                return val_t{ };
        return rOverlays.get( vSparsePos, xOverlays )
            .get( rSparseCoords, rPrefixSums, vPos, actualFromGridPos( rSparseCoords, vSparsePos ), uiCornerIdx
#if GET_PROG_PRINTS
                  ,
                  xProg
#endif
            );
    }

    sps_t getAll( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
                  const prefix_sum_grid_t& rPrefixSums, const pos_t& vPos, progress_stream_t& xProg ) const
    {
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
            .getAll( rSparseCoords, rPrefixSums, vPos, actualFromGridPos( rSparseCoords, vSparsePos ), xProg );
    }

    coordinate_t getNumInternalPrefixSums( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords ) const
    {
        coordinate_t uiRet = 0;
        for( coordinate_t uiI = 0; uiI < rOverlays.sizeOf( xOverlays ); uiI++ )
            uiRet += rOverlays.vData[ uiI + xOverlays.uiStartIndex ].getNumInternalPrefixSums( rSparseCoords );
        return uiRet;
    }

    coordinate_t getNumOverlayPrefixSums( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords ) const
    {
        coordinate_t uiRet = 0;
        for( coordinate_t uiI = 0; uiI < rOverlays.sizeOf( xOverlays ); uiI++ )
            uiRet += rOverlays.vData[ uiI + xOverlays.uiStartIndex ].getNumOverlayPrefixSums( rSparseCoords );
        return uiRet;
    }

    coordinate_t getNumInternalSparseCoords( const overlay_grid_t& rOverlays ) const
    {
        coordinate_t uiRet = 0;
        for( coordinate_t uiI = 0; uiI < rOverlays.sizeOf( xOverlays ); uiI++ )
            uiRet += rOverlays.vData[ uiI + xOverlays.uiStartIndex ].getNumInternalSparseCoords( );
        return uiRet;
    }

    coordinate_t getNumOverlaySparseCoords( const overlay_grid_t& rOverlays ) const
    {
        coordinate_t uiRet = 0;
        for( coordinate_t uiI = 0; uiI < rOverlays.sizeOf( xOverlays ); uiI++ )
            uiRet += rOverlays.vData[ uiI + xOverlays.uiStartIndex ].getNumOverlaySparseCoords( );
        return uiRet;
    }

    coordinate_t getNumGlobalSparseCoords( ) const
    {
        coordinate_t uiRet = 0;

        if constexpr( !UNIFORM_OVERLAY_GRID )
        {
            for( size_t uiI = 0; uiI < D; uiI++ )
                if( uiI != 1 || !DEPENDANT_DIMENSION )
                    uiRet += 1 + vSparseCoords[ uiI ].uiEndCord - vSparseCoords[ uiI ].uiStartCord;
            if constexpr( DEPENDANT_DIMENSION )
                uiRet += xSparseCoordsDependantDimension.uiNum * ( 1 + xSparseCoordsDependantDimension.uiEndCord -
                                                                   xSparseCoordsDependantDimension.uiStartCord );
        }
        return uiRet;
    }

    coordinate_t getNumOverlays( ) const
    {
        return overlay_grid_t::sizeOf( xOverlays );
    }

    coordinate_t getSize( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords ) const
    {
        uint64_t uiNumPSTotal = ( getNumInternalPrefixSums( rOverlays, rSparseCoords ) +
                                  getNumOverlayPrefixSums( rOverlays, rSparseCoords ) ) *
                                sizeof( sps_t );

        uint64_t uiNumLookupTotal = ( getNumInternalSparseCoords( rOverlays ) + getNumOverlaySparseCoords( rOverlays ) +
                                      getNumGlobalSparseCoords( ) ) *
                                    sizeof( coordinate_t );

        uint64_t uiSizeOverlaysOverhead = getNumOverlays( ) * sizeof( overlay_t );

        // full
        return uiNumPSTotal + uiNumLookupTotal + uiSizeOverlaysOverhead + sizeof( Dataset );
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

template <typename type_defs> std::mt19937 Dataset<type_defs>::xGen = std::mt19937( std::random_device( )( ) );

template <typename type_defs>
typename Dataset<type_defs>::points_sort_func_t
    Dataset<type_defs>::sort_points_distinct = typename Dataset<type_defs>::points_sort_func_t( );

} // namespace sps