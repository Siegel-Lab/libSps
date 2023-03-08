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

#if WITH_PYTHON
#include <pybind11/pybind11.h>
#endif

namespace sps
{


template <typename type_defs> class Dataset
{
    EXTRACT_TYPE_DEFS; // macro call

    using sparse_coord_t = SparseCoord<type_defs>;
    using overlay_t = AlignedPower2<Overlay<type_defs>>;
    using overlay_grid_t = typename Overlay<type_defs>::overlay_grid_t;
    using prefix_sum_grid_t = typename Overlay<type_defs>::prefix_sum_grid_t;
    using corner_t = AlignedPower2<Corner<type_defs>>;
    using corners_t = Corners<type_defs>;
    using desc_t = Desc<type_defs>;

    typename overlay_grid_t::template Entry<D> xOverlays;

    pos_t uiSizeOverlays;
    pos_t uiMaxCoords;
    pos_t uiMinCoords;

    static std::mt19937 xGen;

    class DimIterator
    {
        typename corners_t::EntryIterator xIt, xItEnd;
        size_t uiDimension;
        size_t uiBlockSize;

      public:
        DimIterator( typename corners_t::EntryIterator xIt, typename corners_t::EntryIterator xItEnd,
                     size_t uiDimension, size_t uiBlockSize )
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

    struct PointsBinComperator
    {
        const Dataset& rDataset;
        overlay_grid_t& rOverlays;
        std::mutex& xLock;

        PointsBinComperator( const Dataset& rDataset, overlay_grid_t& rOverlays, std::mutex& xLock )
            : rDataset( rDataset ), rOverlays( rOverlays ), xLock( xLock )
        {}

        bool comp( const corner_t& a, const corner_t& b ) const
        {
            coordinate_t uiA = rDataset.overlayIndex( rOverlays, a.vPos );
            coordinate_t uiB = rDataset.overlayIndex( rOverlays, b.vPos );
            return uiA < uiB;
        }

        bool operator( )( const corner_t& a, const corner_t& b ) const
        {
            if constexpr( sparse_coord_t::THREADSAVE )
                return comp( a, b );
            else
            {
                std::lock_guard<std::mutex> xGuard( xLock );
                return comp( a, b );
            }
        }
    };

    struct PointsComperator
    {
        const size_t uiDim;

        PointsComperator( size_t uiDim ) : uiDim( uiDim )
        {}

        bool operator( )( const corner_t& a, const corner_t& b ) const
        {
            return a.vPos[ uiDim ] < b.vPos[ uiDim ];
        }
    };

    struct PointsComperatorDistinct
    {
        bool operator( )( const corner_t& a, const corner_t& b ) const
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

        corner_t min_value( ) const
        {
            corner_t xRet{ };
            return xRet;
        };

        corner_t max_value( ) const
        {
            corner_t xRet{ };
            for( size_t uiI = 0; uiI < D; uiI++ )
                xRet.vPos[ uiI ] = std::numeric_limits<coordinate_t>::max( );
            return xRet;
        };
    };


  public:
    Dataset( )
    {}

    typename sparse_coord_t::Entry
    makeSparseCoords( sparse_coord_t& rSparseCoords, corners_t& vCorners, typename corners_t::Entry xCorners,
                      size_t uiDim, size_t uiNumBlocks, progress_stream_t xProg,
                      coordinate_t uiFixedStart = std::numeric_limits<coordinate_t>::max( ),
                      coordinate_t uiFixedEnd = std::numeric_limits<coordinate_t>::max( ) ) const
    {
        size_t uiNumCoords = 0;
        coordinate_t uiLast = std::numeric_limits<coordinate_t>::max( );
        coordinate_t uiFirst = std::numeric_limits<coordinate_t>::max( );
        std::sort<typename std::vector<corner_t>::iterator, PointsComperator>(
            vCorners.vData.begin( ) + xCorners.uiStartIndex,
            vCorners.vData.begin( ) + xCorners.uiEndIndex,
            PointsComperator( uiDim ) );
        vCorners.iterate(
            [ & ]( const corner_t& xCorner ) {
                coordinate_t uiCurr = xCorner.vPos[ uiDim ];
                if( uiFirst == std::numeric_limits<coordinate_t>::max( ) )
                    uiFirst = uiCurr;
                if( uiCurr != uiLast )
                {
                    uiLast = uiCurr;
                    ++uiNumCoords;
                }
            },
            xCorners );
        size_t uiBlockSize = 1 + ( uiNumCoords - 1 ) / uiNumBlocks;
        xProg << "generating " << uiNumBlocks << " overlays in dimension " << uiDim << " for " << uiNumCoords
              << " different coordinates.\n";
        auto xStart = DimIterator( vCorners.cbegin( xCorners ), vCorners.cend( xCorners ), uiDim, uiBlockSize );
        auto xEnd = DimIterator( vCorners.cend( xCorners ), vCorners.cend( xCorners ), uiDim, uiBlockSize );
        if( uiFixedStart == std::numeric_limits<coordinate_t>::max( ) )
        {
            uiFixedEnd = uiLast;
            uiFixedStart = uiFirst;
        }
        return rSparseCoords.template addStartEnd<true>( xStart, xEnd, uiFixedStart, uiFixedEnd );
    }

    std::array<std::vector<coordinate_t>, D> getPredecessor( const overlay_grid_t& rOverlays,
                                                             const sparse_coord_t& rSparseCoords, pos_t vGridPos,
                                                             pos_t vPosBottomLeftActual, progress_stream_t xProg ) const
    {

        auto vbHasPredecessor = hasPredecessor( rSparseCoords, vGridPos );
        std::array<std::vector<coordinate_t>, D> vPredecessors{ };
        for( size_t uiD = 0; uiD < D; uiD++ )
        {
            if( vbHasPredecessor[ uiD ] )
            {
                --vPosBottomLeftActual[ uiD ];

                vPredecessors[ uiD ].push_back( overlayIndex( rOverlays, vPosBottomLeftActual ) );

                xProg << Verbosity( 3 ) << "predecessor dim " << uiD << " pos " << vPosBottomLeftActual << " index "
                      << vPredecessors[ uiD ].back( ) << "\n";

                ++vPosBottomLeftActual[ uiD ];
            }
            else
                xProg << Verbosity( 2 ) << "predecessor dim " << uiD << " nonexistant\n";
        }
        return vPredecessors;
    }


    coordinate_t generateSparseCoords( overlay_grid_t& rOverlays,
                                       sparse_coord_t& rSparseCoords,
                                       corners_t& vCorners,
                                       progress_stream_t xProg,
                                       std::function<coordinate_t( size_t )>
                                           fDo )
    {
        rSparseCoords.reserve( rSparseCoords.vData.size( ) + ( 2 << 8 ) );

        std::vector<coordinate_t> vPrefixSumSize;
        size_t uiFrom = xOverlays.uiStartIndex;
        size_t uiTo = uiFrom + overlay_grid_t::sizeOf( xOverlays );
        {
            ThreadPool xPool( rOverlays.THREADSAVE && rSparseCoords.THREADSAVE && vCorners.THREADSAVE
                                  ? std::min( (size_t)overlay_grid_t::sizeOf( xOverlays ),
                                              (size_t)std::thread::hardware_concurrency( ) )
                                  : 0 );
            vPrefixSumSize.resize( std::max( xPool.numThreads( ), (size_t)1 ) );

            std::mutex xPrintMutex;


            for( size_t uiOverlayId = uiFrom; uiOverlayId < uiTo; uiOverlayId++ )
                xPool.enqueue(
                    [ & ]( size_t uiTid, size_t uiOverlayId ) {
#if WITH_PYTHON
                        if( PyErr_CheckSignals( ) != 0 ) // allow Ctrl-C cancel
                            throw pybind11::error_already_set( );
#endif
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
                                               corners_t& vCorners,
                                               std::vector<typename corners_t::Entry>& vSplitPoints,
                                               progress_stream_t xProg )
    {
        return generateSparseCoords( rOverlays, rSparseCoords, vCorners, xProg, [ & ]( size_t uiOverlayId ) {
            return rOverlays.vData[ uiOverlayId ].generateInternalSparseCoords(
                rSparseCoords, vCorners, vSplitPoints[ uiOverlayId - xOverlays.uiStartIndex ], xProg );
        } );
    }

    static coordinate_t maxSpaceForOverlayCoords( pos_t uiNumCoordsPerDim, pos_t uiNumOverlaysPerDim )
    {
        return maxSpaceForInternalCoords( uiNumCoordsPerDim, uiNumOverlaysPerDim );
        return ( 1 + maxSpaceForInternalCoords( uiNumCoordsPerDim, uiNumOverlaysPerDim ) ) * ( D - 1 );
    }


    coordinate_t generateOverlaySparseCoords( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords,
                                              corners_t& vCorners, progress_stream_t xProgIn )
    {
        rSparseCoords.reserve( rSparseCoords.vData.size( ) + ( 2 << 8 ) );

        auto xIterator = rOverlays.template genIterator<D>( xOverlays, successors );


        size_t uiNumThreads =
            rSparseCoords.THREADSAVE && rOverlays.THREADSAVE && vCorners.THREADSAVE
                ? std::min( (size_t)overlay_grid_t::sizeOf( xOverlays ), (size_t)std::thread::hardware_concurrency( ) )
                : 0;
        std::vector<coordinate_t> vPrefixSumSize( std::max( (size_t)1, uiNumThreads ) );
        // actually process the overlays
        xIterator.process( uiNumThreads, xProgIn, [ & ]( size_t uiTid, size_t uiOverlayId ) {
#if WITH_PYTHON
            if( PyErr_CheckSignals( ) != 0 ) // allow Ctrl-C cancel
                throw pybind11::error_already_set( );
#endif
            auto xGuard = rSparseCoords.getCapacityGuard( );
            progress_stream_t xProg = xProgIn;
            // get bottom left position (compressed)
            pos_t vGridPos = rOverlays.posOf( uiOverlayId, xOverlays );

            pos_t vPosBottomLeftActual = actualFromGridPos( vGridPos );
            pos_t vPosTopRightActual = actualTopRightFromGridPos( vGridPos );

#ifndef NDEBUG
            xProg << Verbosity( 2 ) << "overlay anchor is " << vGridPos << " actual pos is " << vPosBottomLeftActual
                  << "\n";
#endif

            // collect direct predecessor overlays for each dimension
            std::array<std::vector<coordinate_t>, D> vPredecessors =
                getPredecessor( rOverlays, rSparseCoords, vGridPos, vPosBottomLeftActual, xProg );


            vPrefixSumSize[ uiTid ] += rOverlays.vData[ uiOverlayId ].generateOverlaySparseCoords(
                rOverlays, rSparseCoords, vPredecessors, vPosBottomLeftActual, vPosTopRightActual, xProg );
        } );
        rSparseCoords.shrink_to_fit( );

        coordinate_t uiRet = 0;
        for( coordinate_t uiX : vPrefixSumSize )
            uiRet += uiX;
        return uiRet;
    }

    void generateInternalPrefixSums( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords,
                                     prefix_sum_grid_t& rPrefixSums, corners_t& vCorners,
                                     std::vector<typename corners_t::Entry>& vSplitPoints,
                                     coordinate_t uiSizeInternalPrefixSums, progress_stream_t xProg )
    {
        rPrefixSums.reserve( rPrefixSums.vData.size( ) + uiSizeInternalPrefixSums );

        size_t uiFrom = xOverlays.uiStartIndex;
        size_t uiTo = uiFrom + overlay_grid_t::sizeOf( xOverlays );
        {
            ThreadPool xPool( rSparseCoords.THREADSAVE && rOverlays.THREADSAVE && vCorners.THREADSAVE
                                  ? std::min( (size_t)overlay_grid_t::sizeOf( xOverlays ),
                                              (size_t)std::thread::hardware_concurrency( ) )
                                  : 0 );

            for( size_t uiOverlayId = uiFrom; uiOverlayId < uiTo; uiOverlayId++ )
                xPool.enqueue(
                    [ & ]( size_t uiTid, size_t uiOverlayId, progress_stream_t xProgCpy ) {
#if WITH_PYTHON
                        if( PyErr_CheckSignals( ) != 0 ) // allow Ctrl-C cancel
                            throw pybind11::error_already_set( );
#endif
                        rOverlays.vData[ uiOverlayId ].generateInternalPrefixSums(
                            rOverlays, rSparseCoords, rPrefixSums, vCorners, vSplitPoints[ uiOverlayId - uiFrom ],
                            xProgCpy );

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
                                           const typename overlay_grid_t::template Entry<D>& xEntry )
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

    void generateOverlayPrefixSums( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords,
                                    prefix_sum_grid_t& rPrefixSums, corners_t& vCorners,
                                    coordinate_t uiSizeOverlayPrefixSums, progress_stream_t xProgIn )
    {
        rPrefixSums.reserve( rPrefixSums.vData.size( ) + uiSizeOverlayPrefixSums );

        auto xIterator = rOverlays.template genIterator<D>( xOverlays, successors );
        // actually process the overlays
        xIterator.process(
            rSparseCoords.THREADSAVE && rOverlays.THREADSAVE && rPrefixSums.THREADSAVE
                ? std::min( (size_t)overlay_grid_t::sizeOf( xOverlays ), (size_t)std::thread::hardware_concurrency( ) )
                : 0,
            xProgIn,
            [ & ]( size_t /* uiTid */, size_t uiOverlayId ) {
#if WITH_PYTHON
                if( PyErr_CheckSignals( ) != 0 ) // allow Ctrl-C cancel
                    throw pybind11::error_already_set( );
#endif
                progress_stream_t xProg = xProgIn;
                pos_t vGridPos = rOverlays.posOf( uiOverlayId, xOverlays );

                pos_t vPosBottomLeftActual = actualFromGridPos( vGridPos );

                // collect direct predecessor overlays for each dimension
                std::array<std::vector<coordinate_t>, D> vPredecessors =
                    getPredecessor( rOverlays, rSparseCoords, vGridPos, vPosBottomLeftActual, xProg );

                rOverlays.vData[ uiOverlayId ].generateOverlayPrefixSums(
                    rOverlays, rSparseCoords, rPrefixSums, vCorners, vPredecessors, vPosBottomLeftActual, this, xProg );
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

    static coordinate_t sampleInterval( const corners_t& vCorners, pos_t uiFrom, pos_t uiTo,
                                        typename corners_t::Entry& xSortedPoints, size_t uiD, bool bStart )
    {
        coordinate_t uiPos = bStart ? uiFrom[ uiD ] : uiTo[ uiD ];
        while( ( bStart && uiPos < uiTo[ uiD ] ) || ( !bStart && uiPos > uiFrom[ uiD ] ) )
        {
            if( !bStart ) // doing this here to avoid unsigned integer underflow after zero
                --uiPos;
            const size_t uiD2 = uiD != 0 ? 0 : 1;
            bool bFound = false;
            vCorners.forEqualRange(
                // fBefore
                [ & ]( const corner_t& rP ) {
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
                [ & ]( const corner_t& rP ) {
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
                [ & ]( const corner_t& rP ) {
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

    static pos_t sampleIntervalSize( const corners_t& vCorners, pos_t uiFrom, pos_t uiTo,
                                     std::array<typename corners_t::Entry, D> xSortedPoints )
    {
        pos_t uiNumDistinct;
        for( size_t uiD = 0; uiD < D; uiD++ )
        {
            coordinate_t uiEnd = sampleInterval( vCorners, uiFrom, uiTo, xSortedPoints[ uiD ], uiD, false );
            coordinate_t uiStart = sampleInterval( vCorners, uiFrom, uiTo, xSortedPoints[ uiD ], uiD, true );
            if( uiEnd < uiStart )
                uiNumDistinct[ uiD ] = 0;
            else
                uiNumDistinct[ uiD ] = uiEnd - uiStart;
        }

        return uiNumDistinct;
    }

    static coordinate_t sampleIntervalSize( const corners_t& vCorners, pos_t uiFrom, pos_t uiTo,
                                            typename corners_t::Entry& xSortedPoints, size_t uiD )
    {
        coordinate_t uiEnd = sampleInterval( vCorners, uiFrom, uiTo, xSortedPoints, uiD, false );
        coordinate_t uiStart = sampleInterval( vCorners, uiFrom, uiTo, xSortedPoints, uiD, true );
        if( uiEnd < uiStart )
            return 0;
        else
            return uiEnd - uiStart;
    }

    static coordinate_t sampleNumDistinct( const corners_t& vCorners, pos_t uiFrom, pos_t uiTo,
                                           const uint64_t uiNumSamples, typename corners_t::Entry& xSortedPoints,
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
                vCorners.forEqualRange(
                    // fBefore
                    [ & ]( const corner_t& rP ) {
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
                    [ & ]( const corner_t& rP ) {
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
                    [ & ]( const corner_t& rP ) {
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

    static pos_t sampleNumDistinct( const corners_t& vCorners, pos_t uiFrom, pos_t uiTo, const uint64_t uiNumSamples,
                                    std::array<typename corners_t::Entry, D> xSortedPoints )
    {
        pos_t uiNumDistinct;
        for( size_t uiD = 0; uiD < D; uiD++ )
            uiNumDistinct[ uiD ] = sampleNumDistinct( vCorners, uiFrom, uiTo, uiNumSamples, xSortedPoints[ uiD ], uiD );

        return uiNumDistinct;
    }


    static std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>
    estimateOverlay( const corners_t& vCorners, std::array<typename corners_t::Entry, D> xSortedPoints,
                     pos_t uiCoordinateSizes, pos_t uiNumOverlays, uint64_t, pos_t uiP, pos_t uiMinPos,
                     const uint64_t uiNumPointSamples )
    {
        uint64_t uiNumPSInternal = 1;
        uint64_t uiNumPSOverlay = 0;
        uint64_t uiNumLookUpTablesInternal = 0;
        uint64_t uiNumLookUpTablesOverlay = 0;
        pos_t uiFromGlob, uiToGlob;
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            uiFromGlob[ uiI ] =
                uiP[ uiI ] * ( 1 + ( uiCoordinateSizes[ uiI ] - 1 ) / uiNumOverlays[ uiI ] ) + uiMinPos[ uiI ];
            uiToGlob[ uiI ] =
                ( uiP[ uiI ] + 1 ) * ( 1 + ( uiCoordinateSizes[ uiI ] - 1 ) / uiNumOverlays[ uiI ] ) + uiMinPos[ uiI ];
        }
        pos_t uiNumDistinctInOverlay =
            sampleNumDistinct( vCorners, uiFromGlob, uiToGlob, uiNumPointSamples, xSortedPoints );
        pos_t uiSampledIntervalSize;
        if constexpr( BINARY_SEARCH_BASED_SPARSE )
            uiSampledIntervalSize =
                sampleNumDistinct( vCorners, uiFromGlob, uiToGlob, uiNumPointSamples, xSortedPoints );
        else
            uiSampledIntervalSize = sampleIntervalSize( vCorners, uiFromGlob, uiToGlob, xSortedPoints );


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
                            sampleNumDistinct( vCorners, uiFrom, uiTo, uiNumPointSamples, xSortedPoints[ uiJ ], uiJ );

                        if constexpr( BINARY_SEARCH_BASED_SPARSE )
                            uiNumLookUpTablesOverlay += ( uiP[ uiJ ] > 0 ? 1 : 0 ) +
                                                        sampleNumDistinct( vCorners, uiFrom, uiTo, uiNumPointSamples,
                                                                           xSortedPoints[ uiJ ], uiJ );
                        else
                            uiNumLookUpTablesOverlay +=
                                ( uiP[ uiJ ] > 0 ? 1 : 0 ) +
                                sampleIntervalSize( vCorners, uiFrom, uiTo, xSortedPoints[ uiJ ], uiJ );
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

    static std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>
    estimateDataStructureElements( corners_t& vCorners, std::array<typename corners_t::Entry, D> xSortedPoints,
                                   pos_t uiNumOverlays, pos_t uiCoordinateSizes, pos_t uiMinPos )
    {
        uint64_t uiNumOverlaysTotal = 1;
        for( size_t uiI = 0; uiI < D; uiI++ )
            uiNumOverlaysTotal *= uiNumOverlays[ uiI ];

        uint64_t uiNumOverlaySamples = std::max( 1ul, (uint64_t)std::log2( uiNumOverlaysTotal ) );
        uint64_t uiNumPointSamples =
            std::max( 1ul, 5 * (uint64_t)std::log2( xSortedPoints[ 0 ].uiEndIndex - xSortedPoints[ 0 ].uiStartIndex ) );

        std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> tTotal{ };
        {
            ThreadPool xPool( vCorners.THREADSAVE ? std::min( (size_t)std::thread::hardware_concurrency( ),
                                                              (size_t)uiNumOverlaySamples )
                                                  : 0 );

            std::mutex xResMutex;

            for( size_t uiI = 0; uiI < uiNumOverlaySamples; uiI++ )
                xPool.enqueue( [ & ]( size_t ) {
                    pos_t uiP{ };
                    for( size_t uiI = 0; uiI < D; uiI++ )
                    {
                        std::uniform_int_distribution<uint64_t> xDis( 0, uiNumOverlays[ uiI ] - 1 );
                        uiP[ uiI ] = xDis( xGen );
                    }
                    auto tCurr = estimateOverlay( vCorners, xSortedPoints, uiCoordinateSizes, uiNumOverlays,
                                                  uiNumOverlaysTotal, uiP, uiMinPos, uiNumPointSamples );


                    std::lock_guard<std::mutex> xGuard( xResMutex );
                    std::get<0>( tTotal ) += std::get<0>( tCurr );
                    std::get<1>( tTotal ) += std::get<1>( tCurr );
                    std::get<2>( tTotal ) += std::get<2>( tCurr );
                    std::get<3>( tTotal ) += std::get<3>( tCurr );
                } );
        } // scope for xPool


        uint64_t uiNumPSTotal =
            ( ( uiNumOverlaysTotal * ( std::get<0>( tTotal ) + std::get<1>( tTotal ) ) ) / uiNumOverlaySamples ) *
            sizeof( sps_t );

        uint64_t uiNumLookupTotal =
            ( ( uiNumOverlaysTotal * ( std::get<2>( tTotal ) + std::get<3>( tTotal ) ) ) / uiNumOverlaySamples ) *
            sizeof( coordinate_t );

        uint64_t uiSizeOverlaysOverhead = uiNumOverlaysTotal * sizeof( overlay_t );

        // full
        uint64_t uiSizeFull = uiNumPSTotal + uiNumLookupTotal + uiSizeOverlaysOverhead + sizeof( Dataset );

        return std::make_tuple( ( uiNumOverlaysTotal * std::get<0>( tTotal ) ) / uiNumOverlaySamples,
                                ( uiNumOverlaysTotal * std::get<1>( tTotal ) ) / uiNumOverlaySamples,
                                ( uiNumOverlaysTotal * std::get<2>( tTotal ) ) / uiNumOverlaySamples,
                                ( uiNumOverlaysTotal * std::get<3>( tTotal ) ) / uiNumOverlaySamples,
                                uiNumOverlaysTotal,
                                0, // @todo remove from tuple
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

    static double toFactor( corners_t& vCorners, const typename corners_t::Entry xCorners, uint64_t uiNumOverlays )
    {
        pos_t uiCoordinateSizes = generateCoordSizes( vCorners, xCorners )[ 0 ];
        return toFactor( uiCoordinateSizes, uiNumOverlays );
    }

    static std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>>
    estimateDataStructureElements( corners_t& vCorners, const typename corners_t::Entry xCorners,
                                   std::vector<double> vFac )
    {
        std::array<typename corners_t::Entry, D> xSortedPoints;
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            xSortedPoints[ uiI ] = vCorners.copyEntry( xCorners );

            vCorners.sortByDim( uiI, uiI != 0 ? 0 : 1, xSortedPoints[ uiI ] );
        }

        auto xA = generateCoordSizes( vCorners, xCorners );
        pos_t uiCoordinateSizes = xA[ 0 ];
        pos_t uiMinPos = xA[ 2 ];
        std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>> vRet;
        for( float fFac : vFac )
            vRet.push_back( estimateDataStructureElements( vCorners, xSortedPoints,
                                                           toNumbers( toNumRatios( uiCoordinateSizes ), fFac ),
                                                           uiCoordinateSizes, uiMinPos ) );

        // remove additional entries again
        for( size_t uiI = 0; uiI < D; uiI++ )
            vCorners.popEntry( xSortedPoints[ D - uiI - 1 ] );

        return vRet;
    }

    static coordinate_t toAmount( pos_t vNums )
    {
        coordinate_t uiRet = 1;
        for( size_t uiI = 0; uiI < D; uiI += 1 )
            uiRet *= vNums[ uiI ];

        return uiRet;
    }


    static double pickOverlayFactor( corners_t& vCorners, const typename corners_t::Entry xCorners,
                                     pos_t uiCoordinateSizes, pos_t uiMinPos, progress_stream_t xProg )
    {
        coordinate_t uiArea = 1;

        for( size_t uiI = 0; uiI < D; uiI++ )
            uiArea *= uiCoordinateSizes[ uiI ];

        if( uiArea < 1000 )
        {
            xProg << Verbosity( 0 ) << "data-area smaller than 1,000. Picking factor 0.\n";
            return 0;
        }

        std::array<typename corners_t::Entry, D> xSortedPoints;
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            xSortedPoints[ uiI ] = vCorners.copyEntry( xCorners );

            vCorners.sortByDim( uiI, uiI != 0 ? 0 : 1, xSortedPoints[ uiI ] );
        }


        std::array<double, D> vNumRatios = toNumRatios( uiCoordinateSizes );

        double fEnd = 100;
        double fStart = 0;


        // make sure to place end past the minimum
        while( true )
        {
            xProg << Verbosity( 0 ) << "searching factors [" << fStart << ", " << fEnd << ")\n";
            uint64_t uiEnd = std::get<6>( estimateDataStructureElements(
                vCorners, xSortedPoints, toNumbers( vNumRatios, fEnd ), uiCoordinateSizes, uiMinPos ) );
            uint64_t uiCenter = std::get<6>( estimateDataStructureElements(
                vCorners, xSortedPoints, toNumbers( vNumRatios, fEnd / 2 ), uiCoordinateSizes, uiMinPos ) );
            if( uiEnd > uiCenter + 10000 )
                break;
            fEnd *= 2;
        }

        // fund the minumum between 0 and fEnd
        const double fSampleSteps = 10;
        uint64_t uiMin = std::numeric_limits<uint64_t>::max( );
        for( size_t uiI = 0;
             uiI < 10 && toAmount( toNumbers( vNumRatios, fStart ) ) != toAmount( toNumbers( vNumRatios, fEnd ) );
             uiI++ )
        {
            xProg << Verbosity( 0 ) << "searching factors [" << fStart << ", " << fEnd << ")\n";
            double fMin = 0;
            uiMin = std::numeric_limits<uint64_t>::max( );
            for( double fPos = fStart; fPos <= fEnd; fPos += ( fEnd - fStart ) / fSampleSteps )
            {
                uint64_t uiCurr = std::get<6>( estimateDataStructureElements(
                    vCorners, xSortedPoints, toNumbers( vNumRatios, fPos ), uiCoordinateSizes, uiMinPos ) );
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
            vCorners.popEntry( xSortedPoints[ D - uiI - 1 ] );

        xProg << Verbosity( 0 ) << "picked factor " << fStart << " -> expected index size is "
              << uiMin / (double)std::pow( 10, 9 ) << "Gb\n";

        return ( fStart + fEnd ) / 2;
    }

    static coordinate_t pickNumOverlays( corners_t& vCorners, const typename corners_t::Entry xCorners,
                                         pos_t uiCoordinateSizes, pos_t uiMinPos, progress_stream_t xProg )
    {
        auto vNums = toNumbers( toNumRatios( uiCoordinateSizes ),
                                pickOverlayFactor( vCorners, xCorners, uiCoordinateSizes, uiMinPos, xProg ) );

        return toAmount( vNums );
    }

    static coordinate_t pickNumOverlays( corners_t& vCorners, const typename corners_t::Entry xCorners,
                                         progress_stream_t xProg )
    {
        auto xA = generateCoordSizes( vCorners, xCorners );
        pos_t uiCoordinateSizes = xA[ 0 ];
        pos_t uiMinPos = xA[ 2 ];
        return pickNumOverlays( vCorners, xCorners, uiCoordinateSizes, uiMinPos, xProg );
    }

    static uint64_t getShekelyanEtAlNumBoxes( corners_t& vCorners, const typename corners_t::Entry xCorners )
    {
        pos_t vSizes = generateCoordSizes( vCorners, xCorners )[ 0 ];
        uint64_t uiRet = 1;
        for( coordinate_t uiX : vSizes )
            uiRet *= uiX;
        return uiRet;
    }

    static pos_t pickOverlayNumbers( corners_t& vCorners, const typename corners_t::Entry xCorners,
                                     pos_t uiCoordinateSizes, pos_t uiMinPos, double fFac, progress_stream_t xProg )
    {
        std::array<double, D> vNumRatios = toNumRatios( uiCoordinateSizes );
        if( fFac >= 0 )
        {
            xProg << Verbosity( 0 ) << "Fixed overlay factor.\n";
            return toNumbers( vNumRatios, fFac );
        }
        else if( fFac == -1 )
        {
            xProg << Verbosity( 0 ) << "Trying to predict optimal overlay factor.\n";
            return toNumbers( vNumRatios, pickOverlayFactor( vCorners, xCorners, uiCoordinateSizes, uiMinPos, xProg ) );
        }
        else if( fFac == -2 )
        {
            xProg << Verbosity( 0 ) << "Picking overlay factor based on Shekelyan et. al.'s formula.\n";
            return toNumbers( vNumRatios,
                              toFactor( vCorners, xCorners, getShekelyanEtAlNumBoxes( vCorners, xCorners ) ) );
        }
        else
            throw std::runtime_error(
                "overlay number factor can only be a positive value (incl. 0), -1 (predict optimum using our "
                "approach), or -2 (calculate optimum using Shekelyan et. al.)" );
    }

    static std::array<pos_t, 3> generateCoordSizes( corners_t& vCorners, const typename corners_t::Entry xCorners )
    {
        pos_t uiMaxCoords;
        pos_t uiMinCoords;
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            uiMinCoords[ uiI ] = std::numeric_limits<coordinate_t>::max( );
            uiMaxCoords[ uiI ] = 0;
        }

        vCorners.iterate(
            [ & ]( const corner_t& xCorner ) {
                for( size_t uiI = 0; uiI < D; uiI++ )
                {
                    uiMinCoords[ uiI ] = std::min( uiMinCoords[ uiI ], xCorner.vPos[ uiI ] );
                    uiMaxCoords[ uiI ] = std::max( uiMaxCoords[ uiI ], xCorner.vPos[ uiI ] + 1 );
                }
            },
            xCorners );

        pos_t uiCoordinateSizes;
        for( size_t uiI = 0; uiI < D; uiI++ )
            uiCoordinateSizes[ uiI ] = uiMaxCoords[ uiI ] - uiMinCoords[ uiI ];
        return std::array<pos_t, 3>{ uiCoordinateSizes, uiMaxCoords, uiMinCoords };
    }

    void generateOverlayCoords( overlay_grid_t& rOverlays, corners_t& vCorners, typename corners_t::Entry xCorners,
                                std::vector<typename corners_t::Entry>& vSplitPoints, double fFac,
                                progress_stream_t xProg )
    {
        auto xA = generateCoordSizes( vCorners, xCorners );
        pos_t uiCoordinateSizes = xA[ 0 ];
        uiMaxCoords = xA[ 1 ];
        uiMinCoords = xA[ 2 ];

        xProg << Verbosity( 0 ) << "picking number of overlay boxes.\n";

        pos_t uiNumOverlaysPerDim =
            pickOverlayNumbers( vCorners, xCorners, uiCoordinateSizes, uiMinCoords, fFac, xProg );

        xProg << Verbosity( 0 ) << "generating overlay grid.\n";
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            // 'round up' size to make sure that no points are past the last overlay
            // also make overlays at least size 2
            uiSizeOverlays[ uiI ] = std::max(
                (coordinate_t)2, (coordinate_t)( 1 + ( uiCoordinateSizes[ uiI ] - 1 ) / uiNumOverlaysPerDim[ uiI ] ) );
            xProg << "generating " << uiNumOverlaysPerDim[ uiI ] << " overlays in dimension " << uiI << "\n";
        }

        xOverlays = rOverlays.template add<D, true, true>( uiNumOverlaysPerDim );

        // sort points so that they match the overlay grid order
        xProg << Verbosity( 0 ) << "sorting points into overlays.\n";
        /* Weirdly stxxl::sort seems to be parallel eventhough the access to stxxl::vector is required to be sequential.
         * This does not seem to be a problem for the vector that is currently sorted.
         * However, my PointsBinComperator accesses another stxxl::vector and therefore has to be guarded with a mutex.
         */
        std::mutex xLock;

        std::sort<typename std::vector<corner_t>::iterator, PointsBinComperator>(
            vCorners.vData.begin( ) + xCorners.uiStartIndex, vCorners.vData.begin( ) + xCorners.uiEndIndex,
            PointsBinComperator( *this, rOverlays, xLock ) );

        // generate all overlays
        coordinate_t uiNumTotal = rOverlays.sizeOf( xOverlays );
        vSplitPoints.resize( uiNumTotal );
        for( coordinate_t uiI = 0; uiI < uiNumTotal; uiI++ )
        {
            if( PyErr_CheckSignals( ) != 0 ) // allow Ctrl-C canc
                throw pybind11::error_already_set( );
            vSplitPoints[ uiI ].uiStartIndex = uiI > 0 ? vSplitPoints[ uiI - 1 ].uiEndIndex : xCorners.uiStartIndex;
            vSplitPoints[ uiI ].uiEndIndex = vSplitPoints[ uiI ].uiStartIndex;
            // collect points for overlay uiI
            while( vSplitPoints[ uiI ].uiEndIndex < xCorners.uiEndIndex &&
                   overlayIndex( rOverlays, vCorners.vData[ vSplitPoints[ uiI ].uiEndIndex ].vPos ) ==
                       uiI + xOverlays.uiStartIndex )
                ++vSplitPoints[ uiI ].uiEndIndex;
        }
        assert( vSplitPoints[ uiNumTotal - 1 ].uiEndIndex == xCorners.uiEndIndex );
    }

    Dataset( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords, prefix_sum_grid_t& rPrefixSums,
             corners_t& vCorners, typename corners_t::Entry xCorners, double fFac, progress_stream_t xProg )
    {
        if( xCorners.uiEndIndex == xCorners.uiStartIndex )
        {
            xProg << Verbosity( 0 ) << "done generating empty index.\n";
            return;
        }

        std::vector<typename corners_t::Entry> vSplitPoints;

        // generate the overall sparse coordinates
        generateOverlayCoords( rOverlays, vCorners, xCorners, vSplitPoints, fFac, xProg );

        xProg << Verbosity( 0 ) << "generating internal sparse coords.\n";
        coordinate_t uiSizeInternalPrefixSums =
            generateInternalSparseCoords( rOverlays, rSparseCoords, vCorners, vSplitPoints, xProg );

        xProg << Verbosity( 0 ) << "generating overlay sparse coords.\n";
        coordinate_t uiSizeOverlayPrefixSums = generateOverlaySparseCoords( rOverlays, rSparseCoords, vCorners, xProg );

        xProg << Verbosity( 0 ) << "generating internal prefix sums.\n";
        generateInternalPrefixSums( rOverlays, rSparseCoords, rPrefixSums, vCorners, vSplitPoints,
                                    uiSizeInternalPrefixSums, xProg );

        xProg << Verbosity( 0 ) << "generating overlay prefix sums.\n";
        generateOverlayPrefixSums( rOverlays, rSparseCoords, rPrefixSums, vCorners, uiSizeOverlayPrefixSums, xProg );
        xProg << Verbosity( 0 ) << "done.\n";
    }

    std::array<bool, D> hasPredecessor( const sparse_coord_t& /*rSparseCoords*/, pos_t vPos ) const
    {
        std::array<bool, D> vRet;

        for( size_t uiI = 0; uiI < D; uiI++ )
            vRet[ uiI ] = vPos[ uiI ] > 0;

        return vRet;
    }


    inline coordinate_t actualFromGridPos( coordinate_t uiPos, size_t uiD ) const
    {
        return uiPos * uiSizeOverlays[ uiD ] + uiMinCoords[ uiD ];
    }

    inline pos_t actualFromGridPos( pos_t vPos ) const
    {
        pos_t vRet;
        for( size_t uiI = 0; uiI < D; uiI++ )
            vRet[ uiI ] = actualFromGridPos( vPos[ uiI ], uiI );
        return vRet;
    }

    inline pos_t actualTopRightFromGridPos( pos_t vPos ) const
    {
        for( size_t uiI = 0; uiI < D; uiI++ )
            ++vPos[ uiI ];

        pos_t vRet = actualFromGridPos( vPos );

        return vRet;
    }

    coordinate_t overlayCoord( const coordinate_t uiDataCoord, size_t uiD ) const
    {
        assert( uiSizeOverlays[ uiD ] != 0 );
        if( uiDataCoord < uiMinCoords[ uiD ] )
            return std::numeric_limits<coordinate_t>::max( );
        else
            return std::min( xOverlays.vAxisSizes[ uiD ] - 1,
                             ( uiDataCoord - uiMinCoords[ uiD ] ) / uiSizeOverlays[ uiD ] );
    }

    pos_t overlayCoord( const pos_t& vPos ) const
    {
        pos_t vRet;
        for( size_t uiI = 0; uiI < D; uiI++ )
            vRet[ uiI ] = overlayCoord( vPos[ uiI ], uiI );
        return vRet;
    }

    pos_t overlaySize( const pos_t& vOverlayPos, const pos_t& vSize ) const
    {
        pos_t vRet;
        for( size_t uiI = 0; uiI < D; uiI++ )
            vRet[ uiI ] =
                std::min( xOverlays.vAxisSizes[ uiI ] - 1, vOverlayPos[ uiI ] + vSize[ uiI ] / uiSizeOverlays[ uiI ] );
        return vRet;
    }

    struct OverlayInfo
    {
        pos_t vBottomLeft, vTopRight;
        pos_t vGridPos;
        coordinate_t uiIdx;
        std::array<std::vector<coordinate_t>, D> vPredIds;
        std::vector<pos_t> vvCorners{ };
    };

    std::vector<OverlayInfo> getOverlayInfo( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
                                             const corners_t& /*vCorners*/ ) const
    {
        std::vector<OverlayInfo> vRet;

        for( coordinate_t uiI = 0; uiI < rOverlays.sizeOf( xOverlays ); uiI++ )
        {
            vRet.emplace_back( );
            vRet.back( ).uiIdx = uiI;
            vRet.back( ).vGridPos = rOverlays.posOf( uiI, xOverlays );
            vRet.back( ).vBottomLeft = actualFromGridPos( vRet.back( ).vGridPos );
            vRet.back( ).vTopRight = actualTopRightFromGridPos( vRet.back( ).vGridPos );
            vRet.back( ).vPredIds =
                getPredecessor( rOverlays, rSparseCoords, vRet.back( ).vGridPos, vRet.back( ).vBottomLeft,
                                typename type_defs::progress_stream_t( 0 ) );
            // vCorners.iterate( [ & ]( const corner_t& xP ) { vRet.back( ).vvCorners.push_back( xP.vPos ); },
            //                  rOverlays.vData[ uiI ].xCorners  );
        }

        return vRet;
    }

    std::vector<std::array<pos_t, 3>> getOverlayGrid( const overlay_grid_t& rOverlays ) const
    {
        std::vector<std::array<pos_t, 3>> vRet;

        for( coordinate_t uiI = 0; uiI < rOverlays.sizeOf( xOverlays ); uiI++ )
        {
            vRet.emplace_back( );
            pos_t vGridPos = rOverlays.posOf( uiI, xOverlays );
            vRet.back( )[ 0 ] = vGridPos;
            vRet.back( )[ 1 ] = actualFromGridPos( vGridPos );
            vRet.back( )[ 2 ] = actualTopRightFromGridPos( vGridPos );
        }
        return vRet;
    }

    coordinate_t overlayIndex( const overlay_grid_t& rOverlays, const pos_t& vPos ) const
    {
        pos_t vPosOverlay = overlayCoord( vPos );
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
        auto vSparsePos = overlayCoord( vPos );
#if GET_PROG_PRINTS
        xProg << Verbosity( 2 );
        if( xProg.active( ) )
            xProg << "\t" << vPos << " -> " << vSparsePos << "; that's overlay "
                  << rOverlays.indexOf( vSparsePos, xOverlays ) << "\n";
#endif
        for( size_t uiI = 0; uiI < D; uiI++ )
            if( vSparsePos[ uiI ] == std::numeric_limits<coordinate_t>::max( ) )
                return val_t{ };
        return rOverlays.template get<D, false, false>( vSparsePos, xOverlays )
            .get( rSparseCoords, rPrefixSums, vPos, actualFromGridPos( vSparsePos ), uiCornerIdx
#if GET_PROG_PRINTS
                  ,
                  xProg
#endif
            );
    }

    struct OverlayBounds
    {
        size_t uiGridFrom;
        size_t uiGridTo;
        size_t uiOverlayIdx;
        coordinate_t uiBottomLeft;
        coordinate_t uiTopRight;

        friend std::ostream& operator<<( std::ostream& os, const OverlayBounds& rBounds )
        {
            os << "< grid_from: " << rBounds.uiGridFrom << " grid_to: " << rBounds.uiGridTo
               << " overlay_idx: " << rBounds.uiOverlayIdx << " bottom_left: " << rBounds.uiBottomLeft
               << " top_right: " << rBounds.uiTopRight << " >";

            return os;
        }
    };

    using OverlayBoundsGrid = std::array<std::vector<OverlayBounds>, D>;

    using grid_ret_t = NDGrid<type_defs, val_t, RamVecGenerator>;
    // @todo this needs not to consider orthotope dimensions
    using grid_ret_entry_t = typename grid_ret_t::template Entry<D>;


    class CollectGridValues
    {
      public:
        typename Overlay<type_defs>::grid_ret_t xRet;
        std::array<std::vector<coordinate_t>, D> vvSparsePoss;
        pos_t vNum;
        const typename Overlay<type_defs>::grid_ret_entry_t xRetEntryInternal;
        const OverlayBoundsGrid rOverlayBounds;
        const std::array<typename Overlay<type_defs>::grid_ret_entry_t, D> vRetEntriesOverlay;

      private:
        const overlay_grid_t& rOverlays;
        const typename overlay_grid_t::template Entry<D>& xOverlays;
        const sparse_coord_t& rSparseCoords;
        const prefix_sum_grid_t& rPrefixSums;
        const std::array<std::vector<coordinate_t>, D>& vGrid;
#if GET_PROG_PRINTS
        progress_stream_t& xProg;
#endif

        // no need to initialize
        pos_t vGridFrom;
        pos_t vGridTo;
        pos_t vGridFromOverlay;
        pos_t vGridToOverlay;
        pos_t vOverlayIdx;
        pos_t vOverlayIdxInGrid;
        pos_t vOverlayBottomLeft;
        pos_t vOverlayTopRight;

        template <size_t N> inline void callGridOnOverlays( )
        {
            if constexpr( N < D )
            {
                for( size_t uiI = 0; uiI < rOverlayBounds[ N ].size( ); uiI++ )
                {
                    const OverlayBounds& rBounds = rOverlayBounds[ N ][ uiI ];
                    vGridFrom[ N ] = rBounds.uiGridFrom;
                    vGridTo[ N ] = rBounds.uiGridTo;
                    vGridFromOverlay[ N ] = rBounds.uiGridFrom + uiI;
                    vGridToOverlay[ N ] = rBounds.uiGridTo + uiI + 1;
                    vOverlayIdx[ N ] = rBounds.uiOverlayIdx;
                    vOverlayIdxInGrid[ N ] = uiI;
                    vOverlayBottomLeft[ N ] = rBounds.uiBottomLeft;
                    vOverlayTopRight[ N ] = rBounds.uiTopRight;
                    callGridOnOverlays<N + 1>( );
                }
            }
            else
            {
#if GET_PROG_PRINTS
                xProg << Verbosity( 2 ) << "\tQuerying " //
                      << "overlay " << vOverlayIdx //
                      << "\n\t\tvOverlayIdxInGrid " << vOverlayIdxInGrid //
                      << "\n\t\tvOverlayBottomLeft " << vOverlayBottomLeft //
                      << "\n\t\tvOverlayTopRight " << vOverlayTopRight //
                      << "\n";
#endif

                rOverlays.template get<D, false, false>( vOverlayIdx, xOverlays )
                    .grid( xRet, vvSparsePoss, xRetEntryInternal, vRetEntriesOverlay, rSparseCoords, rPrefixSums,
                           vOverlayBottomLeft, vOverlayTopRight, vGridFrom, vGridTo, vGridFromOverlay, vGridToOverlay,
                           vGrid, vOverlayIdxInGrid
#if GET_PROG_PRINTS
                           ,
                           xProg
#endif
                    );
            }
        }

        static pos_t initVNum( const std::array<std::vector<coordinate_t>, D>& vGrid )
        {
            pos_t vNum;
            for( size_t uiI = 0; uiI < D; uiI++ )
                vNum[ uiI ] = vGrid[ uiI ].size( );
            return vNum;
        }

        static std::array<typename Overlay<type_defs>::grid_ret_entry_t, D>
        initVRetEntriesOverlay( typename Overlay<type_defs>::grid_ret_t& xRet, pos_t vNum,
                                const OverlayBoundsGrid& rOverlayBounds
#if GET_PROG_PRINTS
                                ,
                                [[maybe_unused]] progress_stream_t& xProg
#endif
        )
        {
            std::array<typename Overlay<type_defs>::grid_ret_entry_t, D> vRetEntriesOverlay;
            for( size_t uiI = 0; uiI < D; uiI++ )
                vNum[ uiI ] += rOverlayBounds[ uiI ].size( );
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
                // the overlay values are flat in one dimension
                // instead of removing this flat dimensions, we keep it and have it at a size according to the number of
                // boxes in this dimension. Hence each box can store its overlay values.
                // here, "number of boxes" refers to the number of overlay bounds,
                // i.e. if boxes are skipped they do not get an entry.
                // Also boxes that do not have a query-cell overlapping their bounds do not get an entry.
                size_t uiPrev = vNum[ uiI ];
                vNum[ uiI ] = rOverlayBounds[ uiI ].size( );

                vRetEntriesOverlay[ uiI ] = xRet.template add<D, true, true>( vNum );

#if GET_PROG_PRINTS
                xProg << Verbosity( 4 ) << "\tinitVRetEntriesOverlay size " << vNum << " dim " << uiI << " entry "
                      << vRetEntriesOverlay[ uiI ] << "\n";
#endif
                vNum[ uiI ] = uiPrev;
            }
            return vRetEntriesOverlay;
        }

        static OverlayBoundsGrid initOverlayBounds( const Dataset& rDataset,
                                                    const std::array<std::vector<coordinate_t>, D>& vGrid )
        {
            OverlayBoundsGrid vRet;
            OverlayBounds xCurr;
            for( size_t uiD = 0; uiD < D; uiD++ )
            {
                vRet[ uiD ].reserve( vGrid[ uiD ].size( ) );
                xCurr.uiGridFrom = 0;

                // skip entries before start of index (i.e. before the first point)
                while( xCurr.uiGridFrom < vGrid[ uiD ].size( ) &&
                       ( vGrid[ uiD ][ xCurr.uiGridFrom ] == std::numeric_limits<coordinate_t>::max( ) ||
                         rDataset.overlayCoord( vGrid[ uiD ][ xCurr.uiGridFrom ], uiD ) ==
                             std::numeric_limits<coordinate_t>::max( ) ) )
                    ++xCurr.uiGridFrom;

                while( xCurr.uiGridFrom < vGrid[ uiD ].size( ) )
                {
                    xCurr.uiOverlayIdx = rDataset.overlayCoord( vGrid[ uiD ][ xCurr.uiGridFrom ], uiD );
                    xCurr.uiBottomLeft = rDataset.actualFromGridPos( xCurr.uiOverlayIdx, uiD );
                    xCurr.uiTopRight = rDataset.actualFromGridPos( xCurr.uiOverlayIdx + 1, uiD );
                    xCurr.uiGridTo = xCurr.uiGridFrom + 1;
                    while( xCurr.uiGridTo < vGrid[ uiD ].size( ) &&
                           xCurr.uiOverlayIdx == rDataset.overlayCoord( vGrid[ uiD ][ xCurr.uiGridTo ], uiD ) )
                        ++xCurr.uiGridTo;
                    vRet[ uiD ].push_back( xCurr );
                    xCurr.uiGridFrom = xCurr.uiGridTo;
                }
            }
            return vRet;
        }

      public:
        CollectGridValues( const Dataset& rDataset,
                           const overlay_grid_t& rOverlays,
                           const typename overlay_grid_t::template Entry<D>& xOverlays,
                           const sparse_coord_t& rSparseCoords,
                           const prefix_sum_grid_t& rPrefixSums,
                           const std::array<std::vector<coordinate_t>, D>& vGrid
#if GET_PROG_PRINTS
                           ,
                           progress_stream_t& xProg
#endif
                           )
            : xRet( ".tmp", true ),
              vvSparsePoss{ },
              vNum( initVNum( vGrid ) ),
              xRetEntryInternal( xRet.template add<D, true, true>( vNum ) ),
              rOverlayBounds( initOverlayBounds( rDataset, vGrid ) ),
              vRetEntriesOverlay( initVRetEntriesOverlay( xRet, vNum, rOverlayBounds
#if GET_PROG_PRINTS
                                                          ,
                                                          xProg
#endif
                                                          ) ),
              rOverlays( rOverlays ),
              xOverlays( xOverlays ),
              rSparseCoords( rSparseCoords ),
              rPrefixSums( rPrefixSums ),
              vGrid( vGrid )
#if GET_PROG_PRINTS
              ,
              xProg( xProg )
#endif
        {
#if GET_PROG_PRINTS
            xProg << Verbosity( 3 ) << "collecting grid values in grid of size " << vNum << "\n";
            xProg << "with overlay borders " << rOverlayBounds << "\n";
            xProg << "box size " << rDataset.uiSizeOverlays << "\n";
            xProg << "min coords " << rDataset.uiMinCoords << "\n";
#endif
            for( size_t uiI = 0; uiI < D; uiI++ )
                vvSparsePoss[ uiI ].reserve( vGrid[ uiI ].size( ) );

            callGridOnOverlays<0>( );

#if GET_PROG_PRINTS
            xProg << Verbosity( 3 ) << "grid internal values: " << xRetEntryInternal << "; ";
            xRetEntryInternal.streamOp( xProg, xRet );
            xProg << "\n";
            xProg << "grid overlay values:\n";
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
                xProg << "\tdim " << uiI << ": " << vRetEntriesOverlay[ uiI ] << "; ";
                vRetEntriesOverlay[ uiI ].streamOp( xProg, xRet );
                xProg << "\n";
            }
#endif
        }
    };

    class CollectCellValues
    {
      public:
        grid_ret_t xRet;

      private:
        const std::array<std::vector<coordinate_t>, D>& vGrid;
        const pos_t vNumMin1;
        const grid_ret_entry_t rRetEntry;

        const OverlayBoundsGrid& vOverlayBounds;
        struct ExtractedOverlayBounds
        {
            size_t uiBottomLeftIndex;
            size_t uiIndex;
        };
        using ExtractedOverlayBoundsGrid = std::array<std::vector<ExtractedOverlayBounds>, D>;

        const ExtractedOverlayBoundsGrid vExtractedOverlayBounds;
        const typename Overlay<type_defs>::grid_ret_t& xCollectedVals;
        const typename Overlay<type_defs>::grid_ret_entry_t& xRetEntryInternal;
        const isect_arr_t& vInterTypes;
        const size_t uiFac;

        const std::array<typename Overlay<type_defs>::grid_ret_entry_t, D>& vRetEntriesOverlay;

        size_t uiNumOverlayBorders;

#if GET_PROG_PRINTS
        progress_stream_t& xProg;
#endif

        // no need to initialize
        const typename Overlay<type_defs>::grid_ret_entry_t* pCurrEntry;
        pos_t xGridIdx;
        pos_t xGridIdxBottomLeft;
        // running variable for addGridValuesToCurrCellItr
        pos_t xCollectedIdx;
        val_t* puiCurrCellValue;
        size_t uiCurrOverlay;


        template <size_t uiD, size_t uiN, size_t uiDistToTo> inline void addGridValuesToCurrCellItr( )
        {
            if constexpr( uiD < D )
            {
                xCollectedIdx[ uiD ] = xGridIdx[ uiD ];
                if( xCollectedIdx[ uiD ] != std::numeric_limits<size_t>::max( ) )
                    addGridValuesToCurrCellItr<uiD + 1, uiN, uiDistToTo + 1>( );
                xCollectedIdx[ uiD ] = xGridIdx[ uiD ] + 1;
                addGridValuesToCurrCellItr<uiD + 1, uiN + ( 1 << ( D - ( uiD + 1 ) ) ), uiDistToTo>( );
            }
            else
            {
#if GET_PROG_PRINTS
                xProg << Verbosity( 3 ) << "\t\t\taddGridValuesToCurrCellItr grid_pos: " << xCollectedIdx
                      << " uiN: " << uiN << " uiDistToTo: " << uiDistToTo;
#endif
                const sps_t& uiCurrArr = xCollectedVals.template get<D, false>( xCollectedIdx, *pCurrEntry );
                val_t uiCurr;
                if constexpr( IS_ORTHOTOPE )
                    uiCurr =
                        uiCurrArr[ Overlay<type_defs>::template intersectionTypeToCornerIndex<uiN>( vInterTypes ) ];
                else
                    uiCurr = uiCurrArr;

#if GET_PROG_PRINTS
                xProg << " uiCurr: " << uiCurr << " uiFac: " << ( uiFac == 1 ? "1" : "-1" )
                      << " fac2: " << ( uiDistToTo % 2 == 0 ? "1" : "-1" ) << "\n";
#endif

                *puiCurrCellValue += ( uiDistToTo % 2 == 0 ? 1 : -1 ) * uiFac * uiCurr;
            }
        }

        void addGridValuesToCurrCell( )
        {
            addGridValuesToCurrCellItr<0, 0, 0>( );
        }

        /*
         * @todo instead of nested loops have one loop that sets the coordinates while iterating...
         * for each dimension there can be a divisor and modulator. so we iterate over the raw 1d array.
         * then for each cell compute the coordinates by division & modulo.
         * This would be an almost if-less implementation (there is one if left in the loop termination...)
         *
         * To figure out the box-borders, there should then be a bool vector for each dimension.
         */
        template <size_t uiD> inline void execOnAllCells( )
        {
            if constexpr( uiD < D )
                for( size_t uiX = 0; uiX < vNumMin1[ uiD ]; uiX++ )
                {
                    xGridIdx[ uiD ] = uiX;
                    execOnAllCells<uiD + 1>( );
                }
            else
            {
#if GET_PROG_PRINTS
                xProg << Verbosity( 3 ) << "\t\texecOnAllCells " << xGridIdx << "\n";
#endif
                puiCurrCellValue = &xRet.template get<D, false, false>( xGridIdx, rRetEntry );
                pCurrEntry = &xRetEntryInternal;
                addGridValuesToCurrCell( );
            }
        }

        template <size_t uiO, size_t uiD> inline void execOnAllBorderOverlappingCells( )
        {
            if constexpr( uiD == uiO )
            {
#if GET_PROG_PRINTS
                xProg << Verbosity( 4 ) << "\t\texecOnAllOverlayBorders_4 N: " << uiO << " uiD: " << uiD
                      << " size: " << vOverlayBounds[ uiD ].size( ) << " cell_pos: " << xGridIdx << "\n";
#endif
                if( vOverlayBounds[ uiD ].size( ) > 0 )
                {
                    // first overlay bounds might be the transition from
                    // before the first overlay-box into the first overlay-box.
                    // in that case we need the overlay entries for these grid-cells
                    // why are they not all zero?:
                    // Answer: e.g. the box at x: 0 (<- first box) y: 3 (<- not first box) will have overlay entries,
                    // below it (in y-direction) that are non empty eventhough all the entries to the left (x-direction)
                    // are empty.
                    if( vOverlayBounds[ uiD ][ 0 ].uiGridFrom > 0 )
                    {
                        uiCurrOverlay = std::numeric_limits<size_t>::max( );
                        // -1 since we need the grid-cell that overlaps the two overlay-boxes
                        xGridIdx[ uiD ] = vOverlayBounds[ uiD ][ 0 ].uiGridFrom - 1;
#if GET_PROG_PRINTS
                        xProg << Verbosity( 4 ) << "\t\texecOnAllOverlayBorders_3 N: " << uiO << " uiD: " << uiD
                              << " cell_pos: " << xGridIdx << " uiCurrOverlay: " << uiCurrOverlay << "\n";
#endif
                        execOnAllBorderOverlappingCells<uiO, uiD + 1>( );
                    }
                    for( size_t uiI = 1; uiI < vOverlayBounds[ uiD ].size( ); uiI++ )
                    {
                        uiCurrOverlay = vOverlayBounds[ uiD ][ uiI - 1 ].uiOverlayIdx + uiI;
                        xGridIdx[ uiD ] = vOverlayBounds[ uiD ][ uiI ].uiGridFrom - 1; // see above
#if GET_PROG_PRINTS
                        xProg << Verbosity( 4 ) << "\t\texecOnAllOverlayBorders_2 N: " << uiO << " uiD: " << uiD
                              << " cell_pos: " << xGridIdx << " uiCurrOverlay: " << uiCurrOverlay << "\n";
#endif
                        execOnAllBorderOverlappingCells<uiO, uiD + 1>( );
                    }
                }
            }
            else if constexpr( uiD < D )
                for( size_t uiX = 0; uiX < vNumMin1[ uiD ]; uiX++ )
                {
                    xGridIdx[ uiD ] = uiX;
#if GET_PROG_PRINTS
                    xProg << Verbosity( 4 ) << "\t\texecOnAllOverlayBorders_1 N: " << uiO << " uiD: " << uiD
                          << " cell_pos: " << xGridIdx << "\n";
#endif
                    execOnAllBorderOverlappingCells<uiO, uiD + 1>( );
                }
            else
            {
#if GET_PROG_PRINTS
                xProg << Verbosity( 3 ) << "\t\texecOnAllBorderOverlappingCells cell_pos: " << xGridIdx
                      << " uiO: " << uiO << " uiCurrOverlay: " << uiCurrOverlay;
#endif
                puiCurrCellValue = &xRet.template get<D, false, false>( xGridIdx, rRetEntry );
                coordinate_t uiPrev = xGridIdx[ uiO ];
                xGridIdx[ uiO ] = uiCurrOverlay;
#if GET_PROG_PRINTS
                xProg << Verbosity( 3 ) << " xGridIdx: " << xGridIdx << "\n";
#endif
                pCurrEntry = &vRetEntriesOverlay[ uiO ];
                addGridValuesToCurrCell( );
                xGridIdx[ uiO ] = uiPrev;
            }
        }

        template <size_t N> inline void execOnAllOverlayBorders( )
        {
            if constexpr( N < D )
            {
#if GET_PROG_PRINTS
                xProg << Verbosity( 3 ) << "\texecOnAllOverlayBorders N: " << N << "\n";
#endif
                execOnAllBorderOverlappingCells<N, 0>( );
                execOnAllOverlayBorders<N + 1>( );
            }
        }

        static pos_t initVNumMin1( const std::array<std::vector<coordinate_t>, D>& vGrid )
        {
            pos_t vNumMin1;
            for( size_t uiI = 0; uiI < D; uiI++ )
                vNumMin1[ uiI ] = vGrid[ uiI ].size( ) - 1;
            return vNumMin1;
        }

        static ExtractedOverlayBoundsGrid
        initExtractedOverlayBoundsGrid( [[maybe_unused]] const OverlayBoundsGrid& vOverlayBounds )
        {
            ExtractedOverlayBoundsGrid vRet;
#if 0
            for( size_t uiD = 0; uiD < D; uiI++ )
                if( vOverlayBounds[ uiD ].size( ) > 0 )
                {
                    vRet[uiD].reserve(vOverlayBounds[uiD].back().uiGridTo);
                    for( size_t uiX = 0; uiX < vOverlayBounds[ uiD ][ 0 ].uiGridFrom; uiX++ )
                        vRet[uiD].push_back( ExtractedOverlayBounds{
                            .uiBottomLeftIndex = std::numeric_limits<size_t>::max( ),
                            .uiIndex = std::numeric_limits<size_t>::max( ),
                        } );
                    for( size_t uiI = 0; uiI < vOverlayBounds[ uiD ].size( ); uiI++ )
                    {
                        vRet[uiD].push_back( ExtractedOverlayBounds{
                            .uiBottomLeftIndex = vOverlayBounds[ uiD ][ uiI ].uiGridFrom - 1,
                            .uiIndex = vOverlayBounds[ uiD ][ uiI ].uiGridFrom + uiI, // - 1 + 1
                        } );
                        for(size_t uiX = rBounds.uiGridFrom; uiX < rBounds.uiGridTo; uiX++)
                            vRet[uiD].push_back( ExtractedOverlayBounds{
                                .uiBottomLeftIndex = std::numeric_limits<size_t>::max( ),
                                .uiIndex = uiI + uiX, // - 1 + 1
                            } );
                    }
                }
#endif
            return vRet;
        }

      public:
        CollectCellValues( const std::array<std::vector<coordinate_t>, D>& vGrid,
                           const OverlayBoundsGrid& vOverlayBounds,
                           const typename Overlay<type_defs>::grid_ret_t& xCollectedVals,
                           const typename Overlay<type_defs>::grid_ret_entry_t& xRetEntryInternal,
                           const std::array<typename Overlay<type_defs>::grid_ret_entry_t, D>& vRetEntriesOverlay,
                           const isect_arr_t& vInterTypes
#if GET_PROG_PRINTS
                           ,
                           progress_stream_t& xProg
#endif
                           )
            : xRet( ".tmp", true ),
              vGrid( vGrid ),
              vNumMin1( initVNumMin1( vGrid ) ),
              rRetEntry( xRet.template add<D, false, true>( vNumMin1 ) ),
              vOverlayBounds( vOverlayBounds ),
              vExtractedOverlayBounds( initExtractedOverlayBoundsGrid( vOverlayBounds ) ),
              xCollectedVals( xCollectedVals ),
              xRetEntryInternal( xRetEntryInternal ),
              vInterTypes( vInterTypes ),
              uiFac( Overlay<type_defs>::intersectionTypeToFactor( vInterTypes ) ),
              vRetEntriesOverlay( vRetEntriesOverlay ),
              uiNumOverlayBorders( 0 )
#if GET_PROG_PRINTS
              ,
              xProg( xProg )
#endif
        {
#if GET_PROG_PRINTS
            xProg << Verbosity( 3 ) << "collecting cell values in grid of size " << vNumMin1 << "\n";
#endif
            execOnAllCells<0>( );
            execOnAllOverlayBorders<0>( );

#if GET_PROG_PRINTS
            xProg << Verbosity( 3 ) << "cell values: " << rRetEntry << "; ";
            rRetEntry.streamOp( xProg, xRet );
            xProg << "\n";
#endif
        }
    };


    std::vector<val_t> grid( const overlay_grid_t& rOverlays,
                             const sparse_coord_t& rSparseCoords,
                             const prefix_sum_grid_t& rPrefixSums,
                             const std::array<std::vector<coordinate_t>, D>& vGrid,
                             const isect_arr_t& vInterTypes
#if GET_PROG_PRINTS
                             ,
                             progress_stream_t& xProg
#endif
    ) const
    {

#if GET_PROG_PRINTS
        xProg << Verbosity( 2 ) << "\t"
              << " vGrid " << vGrid //
              << "\n";
#endif

        // work is done in constructor
        CollectGridValues xGridValues( *this, rOverlays, xOverlays, rSparseCoords, rPrefixSums, vGrid
#if GET_PROG_PRINTS
                                       ,
                                       xProg
#endif
        );
        // work is done in constructor
        CollectCellValues xCellValues( vGrid, xGridValues.rOverlayBounds, xGridValues.xRet,
                                       xGridValues.xRetEntryInternal, xGridValues.vRetEntriesOverlay, vInterTypes
#if GET_PROG_PRINTS
                                       ,
                                       xProg
#endif
        );


        return xCellValues.xRet.vData;
    }

    sps_t getAll( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
                  const prefix_sum_grid_t& rPrefixSums, const pos_t& vPos, progress_stream_t& xProg ) const
    {
        for( size_t uiI = 0; uiI < D; uiI++ )
            if( vPos[ uiI ] == std::numeric_limits<coordinate_t>::max( ) )
                return sps_t{ };

        auto vSparsePos = overlayCoord( vPos );
#if GET_PROG_PRINTS
        xProg << Verbosity( 2 );
        if( xProg.active( ) )
            xProg << "\t" << vPos << " -> " << vSparsePos << "; that's overlay "
                  << rOverlays.indexOf( vSparsePos, xOverlays ) << "\n";
#endif
        for( size_t uiI = 0; uiI < D; uiI++ )
            if( vSparsePos[ uiI ] == std::numeric_limits<coordinate_t>::max( ) )
                return sps_t{ };
        return rOverlays.template get<D, false, false>( vSparsePos, xOverlays )
            .getAll( rSparseCoords, rPrefixSums, vPos, actualFromGridPos( vSparsePos ), xProg );
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

    coordinate_t getNumOverlays( ) const
    {
        return overlay_grid_t::sizeOf( xOverlays );
    }

    coordinate_t getSize( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords ) const
    {
        uint64_t uiNumPSTotal = ( getNumInternalPrefixSums( rOverlays, rSparseCoords ) +
                                  getNumOverlayPrefixSums( rOverlays, rSparseCoords ) ) *
                                sizeof( sps_t );

        uint64_t uiNumLookupTotal =
            ( getNumInternalSparseCoords( rOverlays ) + getNumOverlaySparseCoords( rOverlays ) ) *
            sizeof( coordinate_t );

        uint64_t uiSizeOverlaysOverhead = getNumOverlays( ) * sizeof( overlay_t );

        // full
        return uiNumPSTotal + uiNumLookupTotal + uiSizeOverlaysOverhead + sizeof( Dataset );
    }

    friend std::ostream& operator<<( std::ostream& os, const Dataset& rDataset )
    {
        os << "<" << std::endl;

        os << "\txOverlayGrid: ";
        os << rDataset.xOverlays << std::endl;
        os << ">";

        return os;
    }

    std::ostream& stream( std::ostream& os, const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
                          const prefix_sum_grid_t& rPrefixSums, const corners_t& vCorners ) const
    {
        os << "<" << std::endl;

        os << "\txOverlayGridShape: " << xOverlays.vAxisSizes << "\n";
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            os << "\txOverlayGridCoords: dim: " << uiI << " coords: {";
            for( size_t uiX = 0; uiX < xOverlays.vAxisSizes[ uiI ]; uiX++ )
            {
                if( uiX > 0 )
                    os << ", ";
                os << actualFromGridPos( uiX, uiI );
            }
            os << "}\n";
        }
        os << "\txOverlayGrid: ";
        xOverlays.stream( os, rOverlays, rSparseCoords, rPrefixSums, *this, vCorners ) << std::endl;

        os << ">";

        return os;
    }
};

template <typename type_defs> std::mt19937 Dataset<type_defs>::xGen = std::mt19937( std::random_device( )( ) );


} // namespace sps