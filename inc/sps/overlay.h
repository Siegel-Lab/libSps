#pragma once

#include "sps/nd_grid.h"
#include "sps/points.h"
#include "sps/sparse_coordinate.h"
#include "sps/thread_pool.h"
#include <cassert>
#include <functional>
#include <string>


namespace sps
{

template <typename type_defs> class Dataset;

template <typename type_defs> class Overlay
{
    EXTRACT_TYPE_DEFS; // macro call

    using sparse_coord_t = SparseCoord<type_defs>;
  public:
    using prefix_sum_grid_t = NDGrid<type_defs, sps_t, prefix_sums_tmpl_vec_generator_t>;
    using overlay_grid_t = NDGrid<type_defs, AlignedPower2<Overlay>, overlays_tmpl_vec_generator_t>;
  private:


    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;
    using desc_t = Desc<type_defs>;
    using dataset_t = Dataset<type_defs>;

    using cord_it_t = typename sparse_coord_t::EntryIterator;
    using point_it_t = typename points_t::EntryIterator;

    using red_pos_t = std::array<coordinate_t, D - 1>;

    using entry_arr_t = std::array<typename sparse_coord_t::Entry, D>;
    using red_entry_arr_t = std::array<typename sparse_coord_t::Entry, D - 1>;

    class MergableIterator
    {
      public:
        virtual void operator++( )
        {}
        virtual const coordinate_t operator*( ) const
        {
            return 0;
        }
        virtual bool operator!=( const std::shared_ptr<MergableIterator> /*pOther*/ ) const
        {
            return true;
        }

        friend std::ostream& operator<<( std::ostream& os, const MergableIterator& /*rIt*/ )
        {
            return os;
        }
    };

    class CordIterator : public MergableIterator
    {
        using cord_it_t = typename sparse_coord_t::EntryIterator;

        cord_it_t xIt;

      public:
        CordIterator( cord_it_t xIt ) : xIt( xIt )
        {}

        void operator++( )
        {
            ++xIt;
        }

        const coordinate_t operator*( ) const
        {
            return ( *xIt ).first;
        }

        bool operator!=( const std::shared_ptr<MergableIterator> pOther ) const
        {
            auto pCasted = std::dynamic_pointer_cast<CordIterator>( pOther );
            if( pOther == nullptr )
                return true;
            return xIt != pCasted->xIt;
        }

        bool operator!=( const CordIterator& rOther ) const
        {
            return xIt != rOther.xIt;
        }

        friend std::ostream& operator<<( std::ostream& os, const CordIterator& rIt )
        {
            return os << rIt.xIt;
        }
    };
    class MergeIterator
    {
        std::vector<std::shared_ptr<MergableIterator>> vBegin;
        std::vector<std::shared_ptr<MergableIterator>> vEnd;
        coordinate_t uiEnd;

      public:
        MergeIterator( std::vector<std::shared_ptr<MergableIterator>> vBegin,
                       std::vector<std::shared_ptr<MergableIterator>> vEnd, coordinate_t uiEnd )
            : vBegin( vBegin ), vEnd( vEnd ), uiEnd( uiEnd )
        {}

        size_t getSmallestValid( ) const
        {
            size_t uiSmallestValid = 0;
            while( uiSmallestValid < vBegin.size( ) && !( *vBegin[ uiSmallestValid ] != vEnd[ uiSmallestValid ] ) )
                ++uiSmallestValid;
            assert( uiSmallestValid < vBegin.size( ) );
            for( size_t uiI = uiSmallestValid + 1; uiI < vBegin.size( ); uiI++ )
                if( *vBegin[ uiI ] != vEnd[ uiI ] && **vBegin[ uiI ] < **vBegin[ uiSmallestValid ] )
                    uiSmallestValid = uiI;
            return uiSmallestValid;
        }

        void operator++( )
        {
            size_t uiSmallestValid = getSmallestValid( );

            for( size_t uiI = 0; uiI < vBegin.size( ); uiI++ )
                if( uiI != uiSmallestValid && *vBegin[ uiI ] != vEnd[ uiI ] &&
                    !( **vBegin[ uiSmallestValid ] < **vBegin[ uiI ] ) )
                    ++*vBegin[ uiI ];
            ++*vBegin[ uiSmallestValid ];
        }

        const coordinate_t operator*( ) const
        {
            return **vBegin[ getSmallestValid( ) ];
        }

        bool operator!=( const MergeIterator& rOther ) const
        {
            for( size_t uiI = 0; uiI < vBegin.size( ); uiI++ )
                if( *vBegin[ uiI ] != rOther.vBegin[ uiI ] )
                {
                    if( **this >= uiEnd )
                        return false;
                    return true;
                }
            return false;
        }

        friend std::ostream& operator<<( std::ostream& os, const MergeIterator& rIt )
        {
            return os << rIt.vBegin;
        }
    };


    class PointIterator
    {
        point_it_t xIt;
        size_t uiDim;

      public:
        PointIterator( point_it_t xIt, size_t uiDim ) : xIt( xIt ), uiDim( uiDim )
        {}

        void operator++( )
        {
            coordinate_t uiLast = **this;
            while( !xIt.eof( ) && uiLast == **this )
                ++xIt;
        }

        const coordinate_t operator*( ) const
        {
            return xIt->vPos[ uiDim ];
        }

        bool operator!=( const PointIterator& rOther ) const
        {
            return xIt != rOther.xIt;
        }

        friend std::ostream& operator<<( std::ostream& os, const PointIterator& rIt )
        {
            return os;
        }
    };
    class MergeVecIt : public MergableIterator
    {
        using it_t = typename std::vector<coordinate_t>::iterator;

        it_t xIt;

      public:
        MergeVecIt( it_t xIt ) : xIt( xIt )
        {}

        void operator++( )
        {
            ++xIt;
        }

        const coordinate_t operator*( ) const
        {
            return *xIt;
        }

        bool operator!=( const std::shared_ptr<MergableIterator> pOther ) const
        {
            auto pCasted = std::dynamic_pointer_cast<MergeVecIt>( pOther );
            if( pOther == nullptr )
                return true;
            return xIt != pCasted->xIt;
        }

        friend std::ostream& operator<<( std::ostream& os, const MergeVecIt& rIt )
        {
            return os;
        }
    };


  public:
    std::array<red_entry_arr_t, D> vSparseCoordsOverlay;
    entry_arr_t vSparseCoordsInternal;
    std::array<typename prefix_sum_grid_t::template Entry<D - 1>, D> vOverlayEntries;
    typename prefix_sum_grid_t::template Entry<D> xInternalEntires;
    typename points_t::Entry xPoints;

    Overlay( ) : vSparseCoordsOverlay{ }, vSparseCoordsInternal{ }, vOverlayEntries{ }, xInternalEntires{ }, xPoints{ }
    {}

    void iterate( coordinate_t uiEnd, std::function<void( coordinate_t )> fDo ) const
    {
        for( coordinate_t uiI = 0; uiI < uiEnd; uiI++ )
            fDo( uiI );
    }

    template <size_t I, size_t N>
    inline void iterateHelper( const std::array<coordinate_t, N>& rEnds,
                               std::function<void( const std::array<coordinate_t, N>& )>
                                   fDo,
                               std::array<coordinate_t, N>& rCurr ) const
    {
        if constexpr( I == N )
            fDo( rCurr );
        else
            iterate( rEnds[ I ], [ & ]( coordinate_t uiCurr ) {
                rCurr[ I ] = uiCurr;
                iterateHelper<I + 1, N>( rEnds, fDo, rCurr );
            } );
    }

    template <size_t N>
    void iterate( const std::array<coordinate_t, N>& rEnds,
                  std::function<void( const std::array<coordinate_t, N>& )> fDo ) const
    {
        std::array<coordinate_t, N> rCurr;
        iterateHelper<0, N>( rEnds, fDo, rCurr );
    }

//#ifndef NDEBUG
#pragma GCC diagnostic push
// multiple unused variabels in release mode
#pragma GCC diagnostic ignored "-Wunused-parameter"
    //#endif
    void generate( const overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords, prefix_sum_grid_t& rPrefixSums,
                   points_t& vPoints, typename points_t::Entry xPoints,
                   std::array<std::vector<coordinate_t>, D> vPredecessors, pos_t vMyBottomLeft, pos_t vPosTopRight,
                   dataset_t* pDataset, size_t uiOverlaysNow, size_t uiOverlaysTotal, progress_stream_t& xProg,
                   Profiler& xProfiler, ThreadPool& rPool )
    {
        this->xPoints = xPoints;
#ifndef NDEBUG
        xProfiler.step( "overlay coord construction" );
        // construct sparse coordinates for each dimension
        xProg << Verbosity( 1 ) << "constructing sparse coordinates for overlay bottom left= " << vMyBottomLeft
              << " top right " << vPosTopRight << "\n";
#endif
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
#ifndef NDEBUG
            xProg << Verbosity( 2 ) << "dim " << uiI << "\n";
#endif
            for( size_t uiJ = 0; uiJ < D - 1; uiJ++ )
            {
                size_t uiJAct = uiJ + ( uiJ >= uiI ? 1 : 0 );
#ifndef NDEBUG
                xProg << Verbosity( 2 ) << "sub dim " << uiJAct << " (" << uiJ << ")"
                      << "\n";
#endif
                std::vector<std::shared_ptr<MergableIterator>> vBegin{ };
                std::vector<std::shared_ptr<MergableIterator>> vEnd{ };
                std::vector<coordinate_t> vCollectedCoords;

                // add coordinates from previous overlays to the overlay entries
                for( size_t uiI2 = 0; uiI2 < D; uiI2++ )
                    if( uiJAct != uiI2 )
                        for( coordinate_t uiPred : vPredecessors[ uiI2 ] )
                        {
                            const Overlay* pPred = &rOverlays.vData[ uiPred ];
                            vBegin.push_back( std::make_shared<CordIterator>(
                                rSparseCoords.cbegin( pPred->vSparseCoordsOverlay[ uiI ][ uiJ ] ) ) );
                            vEnd.push_back( std::make_shared<CordIterator>(
                                rSparseCoords.cend( pPred->vSparseCoordsOverlay[ uiI ][ uiJ ] ) ) );
#ifndef NDEBUG
                            if( xProg.active( ) )
                                pPred->vSparseCoordsOverlay[ uiI ][ uiJ ].stream(
                                    std::cout << "from predecessor in dim " << uiI2 << " overlay: ", rSparseCoords )
                                    << std::endl;
#endif
                            if constexpr( DEPENDANT_DIMENSION )
                                if( uiI == 1 && uiI2 != 1 /* <- cause this is done below anyways */ )
                                // predecessors in dim 1 could reach below the start of this overlay ->
                                // in that case their overlay coords are not sufficient & we also have to
                                // use their internal coords (taking only the relevant coords from the points)
                                {
                                    vPoints.iterate(
                                        [ & ]( const point_t& rP ) {
                                            if( rP.vPos[ uiJAct ] >= vMyBottomLeft[ uiJAct ] &&
                                                rP.vPos[ uiJAct ] < vPosTopRight[ uiJAct ] )
                                                vCollectedCoords.push_back( rP.vPos[ uiJAct ] );
                                        },
                                        pPred->xPoints );
                                }
                        }

                if constexpr( DEPENDANT_DIMENSION )
                {
                    std::sort( vCollectedCoords.begin( ), vCollectedCoords.end( ) );
                    vBegin.push_back( std::make_shared<MergeVecIt>( vCollectedCoords.begin( ) ) );
                    vEnd.push_back( std::make_shared<MergeVecIt>( vCollectedCoords.end( ) ) );
                }

                // add coordinates from the points of the previous overlay to the overlay entries
                for( coordinate_t uiPred : vPredecessors[ uiI ] )
                {
                    const Overlay* pPred = &rOverlays.vData[ uiPred ];
                    vBegin.push_back( std::make_shared<CordIterator>(
                        rSparseCoords.cbegin( pPred->vSparseCoordsInternal[ uiJAct ] ) ) );
                    vEnd.push_back( std::make_shared<CordIterator>(
                        rSparseCoords.cend( pPred->vSparseCoordsInternal[ uiJAct ] ) ) );
#ifndef NDEBUG
                    if( xProg.active( ) )
                        pPred->vSparseCoordsInternal[ uiJAct ].stream( std::cout << "from internal: ", rSparseCoords )
                            << std::endl;
#endif
                }

#ifndef NDEBUG
                xProg << "from bottom left: { " << vMyBottomLeft[ uiJAct ] - 1 << " }\n";
#endif

                MergeIterator xBegin( vBegin, vEnd, vPosTopRight[ uiJAct ] );
                MergeIterator xEnd( vEnd, vEnd, vPosTopRight[ uiJAct ] );

                {
                    auto xPartialLock = rSparseCoords.xLockable.partialLock( );
                    // skip over positions that are before and after this overlay
                    while( xBegin != xEnd && *xBegin < vMyBottomLeft[ uiJAct ] )
                        ++xBegin;
                } // end of scope xPartialLock

                {
                    auto xFullLock = rSparseCoords.xLockable.fullLock( );

                    if( vPredecessors[ uiJAct ].size( ) > 0 )
                        vSparseCoordsOverlay[ uiI ][ uiJ ] =
                            rSparseCoords.addStart( xBegin, xEnd, vMyBottomLeft[ uiJAct ] - 1 );
                    else
                        vSparseCoordsOverlay[ uiI ][ uiJ ] = rSparseCoords.add( xBegin, xEnd );

#ifndef NDEBUG
                    xProg << Verbosity( 2 );
                    if( xProg.active( ) )
                        vSparseCoordsOverlay[ uiI ][ uiJ ].stream( std::cout << "result: ", rSparseCoords )
                            << std::endl;
#endif
                } // end of scope xFullLock
            }
        }

        if( xPoints.size( ) > 0 )
        {
#ifndef NDEBUG
            xProfiler.step( "internal coord construction" );
            xProg << Verbosity( 1 ) << "constructing sparse coordinates for points\n";
#endif
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
#ifndef NDEBUG
                xProg << Verbosity( 2 ) << "dim " << uiI << "\n";
#endif
                vPoints.sortByDim( uiI, xPoints );

#ifndef NDEBUG
                if( xProg.active( ) )
                    xPoints.stream( std::cout << "from points: ", vPoints ) << std::endl;
#endif

                auto xFullLock = rSparseCoords.xLockable.fullLock( );
                vSparseCoordsInternal[ uiI ] = rSparseCoords.add( PointIterator( vPoints.cbegin( xPoints ), uiI ),
                                                                  PointIterator( vPoints.cend( xPoints ), uiI ) );
#ifndef NDEBUG
                if( xProg.active( ) )
                    vSparseCoordsInternal[ uiI ].stream( std::cout << "result: ", rSparseCoords ) << std::endl;
#endif
                // end of scope xFullLock
            }

            // construct internal grid
#ifndef NDEBUG
            xProg << Verbosity( 1 ) << "constructing internal grid\n";
#endif
            pos_t vInternalAxisSizes;
            {
                auto xPartialLock = rSparseCoords.xLockable.partialLock( );
                vInternalAxisSizes = rSparseCoords.axisSizes( vSparseCoordsInternal );
            }
#ifndef NDEBUG
            xProg << Verbosity( 2 ) << "axis sizes: " << vInternalAxisSizes << "\n";
#endif
            {
                auto xFullLock = rPrefixSums.xLockable.fullLock( );
                xInternalEntires = rPrefixSums.add( vInternalAxisSizes );
            }

            coordinate_t uiNumTotal = prefix_sum_grid_t::sizeOf( xInternalEntires );
            coordinate_t uiNumDone = 0;
            std::vector<sps_t> vTmp( rPrefixSums.THREADSAVE ? 0 : uiNumTotal, sps_t{ } );

            {
                auto xPartialLock1 = rSparseCoords.xLockable.partialLock( );
                auto xPartialLock2 = rPrefixSums.xLockable.partialLock( );
                vPoints.iterate(
                    [ & ]( const point_t& xPoint ) {
                        size_t uiIdx = prefix_sum_grid_t::indexOf(
                            rSparseCoords.sparse( xPoint.vPos, vSparseCoordsInternal ), xInternalEntires );
                        if constexpr( rPrefixSums.THREADSAVE )
                            xPoint.addTo( rPrefixSums.vData[ uiIdx ] );
                        else
                        {
                            uiIdx -= xInternalEntires.uiStartIndex;
                            assert( uiIdx < vTmp.size( ) );
                            xPoint.addTo( vTmp[ uiIdx ] );
                        }
                    },
                    xPoints );
            }

#ifndef NDEBUG
            xProg << "vSparseCoordsOverlay " << vSparseCoordsOverlay << "\n";
            xProg << "vSparseCoordsInternal " << vSparseCoordsInternal << "\n";
            xProg << "rSparseCoords " << rSparseCoords << "\n";


            xProfiler.step( "prefix sum computation" );
#endif
            // compute internal prefix sum
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
#ifndef NDEBUG
                xProg << Verbosity( 3 ) << "computing prefix sums over dimension " << uiI << "\n";
#endif
                red_pos_t vRelevantInternalAxisSizes = relevant( vInternalAxisSizes, uiI );
                std::mutex xLock;
                std::condition_variable xCv;
                size_t uiEnqueuedTasks = 0;

                auto xPartialLock2 = rPrefixSums.xLockable.partialLock( );

                std::vector<red_pos_t> vvTos;
                std::vector<coordinate_t> vuiTos;
                iterate<D - 1>( vRelevantInternalAxisSizes, [ & ]( const red_pos_t& vTo ) {
                    {
                        std::unique_lock xGuard( xLock );
                        ++uiEnqueuedTasks;
                    }
                    rPool.enqueue(
                        [ & ]( size_t uiTid, red_pos_t vTo, size_t uiI ) {
                            pos_t vFullTo = expand( vTo, uiI );
                            sps_t uiPrefixSum = { };

#ifndef NDEBUG
                            if( uiTid == 0 )
                                xProg << Verbosity( 3 ) << "starting...: " << vFullTo << ": " << uiPrefixSum << "\n";
#endif

                            coordinate_t uiNumDoneLocal = 0;
                            iterate( vInternalAxisSizes[ uiI ], [ & ]( coordinate_t uiTo ) {
                                vFullTo[ uiI ] = uiTo;
                                size_t uiIdx = prefix_sum_grid_t::indexOf( vFullTo, xInternalEntires );
                                if constexpr( !rPrefixSums.THREADSAVE )
                                {
                                    uiIdx -= xInternalEntires.uiStartIndex;
                                    assert( uiIdx < vTmp.size( ) );
                                }

                                if constexpr( rPrefixSums.THREADSAVE )
                                    uiPrefixSum += rPrefixSums.vData[ uiIdx ];
                                else
                                    uiPrefixSum += vTmp[ uiIdx ];

#ifndef NDEBUG
                                if( uiTid == 0 )
                                    xProg << Verbosity( 3 ) << vFullTo << ": " << uiPrefixSum << "\n";
#endif

                                if constexpr( rPrefixSums.THREADSAVE )
                                    rPrefixSums.vData[ uiIdx ] = uiPrefixSum;
                                else
                                    vTmp[ uiIdx ] = uiPrefixSum;

                                ++uiNumDoneLocal;
                            } );

                            std::unique_lock xGuard( xLock );
                            uiNumDone += uiNumDoneLocal;
                            --uiEnqueuedTasks;
                            if( uiEnqueuedTasks == 0 )
                                xCv.notify_one( );
#ifndef NDEBUG
                            if( uiTid == 0 && xProg.printAgain( ) )
                                xProg << Verbosity( 0 ) << "computed " << uiOverlaysNow << " out of " << uiOverlaysTotal
                                      << " overlays, thats "
                                      << 100.0 * ( (double)uiOverlaysNow / (double)uiOverlaysTotal ) << "%. "
                                      << uiNumDone << " out of " << uiNumTotal * D << " prefix sums, thats "
                                      << 100.0 * ( (double)uiNumDone / (double)( uiNumTotal * D ) ) << "%.\n";
#endif
                        },
                        vTo, uiI );
                    xCv.notify_one( );
                } );
                {
                    std::unique_lock xGuard( xLock );
                    if( uiEnqueuedTasks > 0 )
                        xCv.wait( xGuard );
                    assert( uiEnqueuedTasks == 0 );
                }
            }
            if constexpr( !rPrefixSums.THREADSAVE )
            {
#ifndef NDEBUG
                xProfiler.step( "copying prefix sums" );
#endif
                auto xPartialLock = rPrefixSums.xLockable.partialLock( );
                for( size_t uiIdx = 0; uiIdx < uiNumTotal; uiIdx++ )
                    rPrefixSums.vData[ uiIdx + xInternalEntires.uiStartIndex ] = vTmp[ uiIdx ];
            }
        }

        // construct overlay sum grid
#ifndef NDEBUG
        xProg << Verbosity( 1 ) << "constructing overlay sum grid\n";
        xProfiler.step( "filling overlay" );
#endif
        for( size_t uiI = 0; uiI < D; uiI++ )
            if( vPredecessors[ uiI ].size( ) > 0 )
            {
#ifndef NDEBUG
                xProg << Verbosity( 2 ) << "dim " << uiI << "\n";
#endif

                red_pos_t vAxisSizes;
                {
                    auto xPartialLock1 = rSparseCoords.xLockable.partialLock( );
                    vAxisSizes = rSparseCoords.axisSizes( vSparseCoordsOverlay[ uiI ] );
                }
                {
                    auto xFullLock = rPrefixSums.xLockable.fullLock( );
                    vOverlayEntries[ uiI ] = rPrefixSums.add( vAxisSizes );
                }

                assert( vMyBottomLeft[ uiI ] > 0 );
                if( prefix_sum_grid_t::sizeOf( vOverlayEntries[ uiI ] ) > 0 )
                {
                    auto xPartialLock1 = rSparseCoords.xLockable.partialLock( );
                    auto xPartialLock2 = rPrefixSums.xLockable.partialLock( );

                    rSparseCoords.template iterate<D - 1>(
                        [ & ]( const red_pos_t& vFrom, const red_pos_t& vTo ) {
                            pos_t vFullFrom = expand( vFrom, uiI );
                            vFullFrom[ uiI ] = vMyBottomLeft[ uiI ] - 1;

#ifndef NDEBUG
                            xProg << Verbosity( 3 ) << "query " << vFullFrom << "\n";
#endif
                            auto uiRet = pDataset->get( rOverlays, rSparseCoords, rPrefixSums, vFullFrom, xProg );
#ifndef NDEBUG
                            xProg << Verbosity( 3 ) << "query " << vFullFrom << ": " << uiRet << "\n";
#endif

                            rPrefixSums.get( vTo, vOverlayEntries[ uiI ] ) = uiRet;
                        },
                        vSparseCoordsOverlay[ uiI ] );
                }
            }
#ifndef NDEBUG
        xProg << Verbosity( 1 ) << "done\n";
#endif
    }
//#ifndef NDEBUG
#pragma GCC diagnostic pop
    //#endif


    sps_t get( const sparse_coord_t& rSparseCoords, const prefix_sum_grid_t& rPrefixSums, pos_t vCoords,
               pos_t vMyBottomLeft, progress_stream_t& xProg ) const
    {
        sps_t uiRet{ };

        xProg << Verbosity( 2 ) << "\tquerying overlay for " << vCoords << "...\n";

        for( size_t uiI = 0; uiI < D; uiI++ )
            --vMyBottomLeft[ uiI ]; // will turn zero values into max values
        xProg << Verbosity( 3 ) << "\t\tvMyBottomLeft " << vMyBottomLeft << "\n";

        forAllCombinations<pos_t>(
            [ & ]( size_t, pos_t vPos, size_t uiDistToTo ) {
                if( uiDistToTo == 0 )
                    return;
                size_t uiI = 0;
                while( uiI < D &&
                       ( vPos[ uiI ] != vMyBottomLeft[ uiI ] ||
                         vSparseCoordsOverlay[ uiI ][ 0 ].uiStartIndex == std::numeric_limits<coordinate_t>::max( ) ) )
                    ++uiI;

                if( uiI == D )
                    return;

                xProg << Verbosity( 3 ) << "\t\tquery: " << vPos << " in overlay " << uiI << "\n";

                red_pos_t vRelevant = relevant( vPos, uiI );
                red_pos_t vSparse = rSparseCoords.sparse( vRelevant, vSparseCoordsOverlay[ uiI ] );

                xProg << "\t\trelevant: " << vRelevant << " sparse: " << vSparse << "\n";

                sps_t uiCurr = rPrefixSums.get( vSparse, vOverlayEntries[ uiI ] );

                xProg << "\t\tis " << ( uiDistToTo % 2 == 0 ? "-" : "+" ) << uiCurr << "\n";
                uiRet += uiCurr * ( uiDistToTo % 2 == 0 ? -1 : 1 );
            },
            vMyBottomLeft,
            vCoords,
            []( coordinate_t uiPos ) { return uiPos != std::numeric_limits<coordinate_t>::max( ); } );

        if( rPrefixSums.sizeOf( xInternalEntires ) > 0 )
        {
            auto vSparseCoords = rSparseCoords.sparse( vCoords, vSparseCoordsInternal );
            auto uiCurr = rPrefixSums.get( vSparseCoords, xInternalEntires );
            xProg << Verbosity( 2 ) << "\tquerying internal " << vCoords << " -> " << vSparseCoords << ": +" << uiCurr
                  << "\n";
            uiRet += uiCurr;
        }
        else
            xProg << Verbosity( 2 ) << "\tno internal entries\n";


        return uiRet;
    }

    template <typename T> std::array<T, D - 1> relevant( const std::array<T, D>& vAllEntries, size_t uiI ) const
    {
        std::array<T, D - 1> vRelevantEntries;
        for( size_t uiJ = 0; uiJ < uiI; uiJ++ )
            vRelevantEntries[ uiJ ] = vAllEntries[ uiJ ];
        for( size_t uiJ = uiI + 1; uiJ < D; uiJ++ )
            vRelevantEntries[ uiJ - 1 ] = vAllEntries[ uiJ ];
        return vRelevantEntries;
    }

    template <typename T> std::array<T, D> expand( const std::array<T, D - 1>& vCompressed, size_t uiI ) const
    {
        std::array<T, D> vAllEntries;
        for( size_t uiJ = 0; uiJ < uiI; uiJ++ )
            vAllEntries[ uiJ ] = vCompressed[ uiJ ];
        for( size_t uiJ = uiI + 1; uiJ < D; uiJ++ )
            vAllEntries[ uiJ ] = vCompressed[ uiJ - 1 ];
        return vAllEntries;
    }

    friend std::ostream& operator<<( std::ostream& os, const Overlay& rOverlay )
    {
        os << "<" << std::endl;
        os << "\tvSparseCoordsOverlay: ";
        os << rOverlay.vSparseCoordsOverlay << std::endl;

        os << "\tvSparseCoordsInternal: ";
        os << rOverlay.vSparseCoordsInternal << std::endl;

        os << "\tvOverlayEntries: ";
        os << rOverlay.vOverlayEntries << std::endl;

        os << "\txInternalEntires: ";
        os << rOverlay.xInternalEntires << std::endl;

        os << "\txPoints: ";
        os << rOverlay.xPoints << std::endl;

        os << ">";

        return os;
    }

    std::ostream& stream( std::ostream& os, pos_t vGridPos, const sparse_coord_t& rSparseCoords,
                          const prefix_sum_grid_t& rPrefixSums, const dataset_t& rData ) const
    {
        os << "<";
        if( rData.exists( rSparseCoords, vGridPos ) )
        {
            os << std::endl;
            auto vPosActual = rData.actualFromGridPos( rSparseCoords, vGridPos );
            os << "\tbottom left: " << vPosActual << "\n";
            auto vPosTR = rData.actualTopRightFromGridPos( rSparseCoords, vGridPos );
            os << "\ttop right: " << vPosTR << "\n";
            os << "\tvSparseCoordsOverlay: ";
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
                os << uiI << ":<";
                for( size_t uiJ = 0; uiJ < D - 1; uiJ++ )
                {
                    os << " " << ( uiJ + ( uiJ >= uiI ? 1 : 0 ) ) << ":";
                    vSparseCoordsOverlay[ uiI ][ uiJ ].stream( os, rSparseCoords ) << " ";
                }
                os << "> ";
            }
            os << std::endl;

            os << "\tvSparseCoordsInternal: ";
            for( size_t uiI = 0; uiI < D; uiI++ )
                vSparseCoordsInternal[ uiI ].stream( os, rSparseCoords ) << " ";
            os << std::endl;

            os << "\tvOverlayEntries: ";
            for( size_t uiI = 0; uiI < D; uiI++ )
                vOverlayEntries[ uiI ].streamOp( os, rPrefixSums ) << " ";
            os << std::endl;

            os << "\txInternalEntires: ";
            xInternalEntires.streamOp( os, rPrefixSums ) << std::endl;
        }
        else
            os << "n/a";

        return os;
    }
    std::ostream& stream( std::ostream& os, pos_t vGridPos, const sparse_coord_t& rSparseCoords,
                          const prefix_sum_grid_t& rPrefixSums, const dataset_t& rData, const points_t& vPoints ) const
    {
        this->stream( os, vGridPos, rSparseCoords, rPrefixSums, rData );

        if( rData.exists( rSparseCoords, vGridPos ) )
        {
            os << "\txPoints: ";
            xPoints.stream( os, vPoints ) << std::endl;
        }

        os << ">";

        return os;
    }

    std::ostream& stream( std::ostream& os, pos_t vGridPos, const sparse_coord_t& rSparseCoords,
                          const prefix_sum_grid_t& rPrefixSums, const dataset_t& rData, const points_t& vPoints,
                          const desc_t& vDesc ) const
    {
        this->stream( os, vGridPos, rSparseCoords, rPrefixSums, rData );

        if( rData.exists( rSparseCoords, vGridPos ) )
        {
            os << "\txPoints: ";
            xPoints.stream( os, vPoints, vDesc ) << std::endl;
        }

        os << ">";

        return os;
    }
};


} // namespace sps
