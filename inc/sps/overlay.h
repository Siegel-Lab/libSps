#pragma once

#include "sps/corners.h"
#include "sps/nd_grid.h"
#include "sps/sparse_coordinate.h"
#include "sps/thread_pool.h"
#include <bitset>
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
    using corner_t = Corner<type_defs>;
    using corners_t = Corners<type_defs>;
    using desc_t = Desc<type_defs>;
    using dataset_t = Dataset<type_defs>;

    using cord_it_t = typename sparse_coord_t::EntryIterator;
    using point_it_t = typename corners_t::EntryIterator;

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
        virtual const coordinate_t subtract( const std::shared_ptr<MergableIterator> /*pOther*/ ) const
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
        virtual std::shared_ptr<MergableIterator> copy( ) const
        {
            return nullptr;
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

        const coordinate_t subtract( const std::shared_ptr<MergableIterator> pOther ) const
        {
            auto pCasted = std::dynamic_pointer_cast<CordIterator>( pOther );
            if( pOther == nullptr )
                return 0;
            return xIt.subtract( pCasted->xIt );
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

        std::shared_ptr<MergableIterator> copy( ) const
        {
            return std::make_shared<CordIterator>( cord_it_t( xIt ) );
        }
    };

    class MergeIterator
    {
        std::vector<std::shared_ptr<MergableIterator>> vBegin;
        std::vector<std::shared_ptr<MergableIterator>> vEnd;
        coordinate_t uiEnd;
        size_t uiSmallestValid = std::numeric_limits<size_t>::max( );

      public:
        MergeIterator( std::vector<std::shared_ptr<MergableIterator>> vBegin,
                       std::vector<std::shared_ptr<MergableIterator>> vEnd, coordinate_t uiEnd )
            : vBegin( vBegin ), vEnd( vEnd ), uiEnd( uiEnd )
        {
            setSmallestValid( );
        }

        void setSmallestValid( )
        {
            if( vBegin.size( ) > 0 )
            {
                uiSmallestValid = 0;
                while( uiSmallestValid < vBegin.size( ) && !( *vBegin[ uiSmallestValid ] != vEnd[ uiSmallestValid ] ) )
                    ++uiSmallestValid;
                for( size_t uiI = uiSmallestValid + 1; uiI < vBegin.size( ); uiI++ )
                    if( *vBegin[ uiI ] != vEnd[ uiI ] && **vBegin[ uiI ] < **vBegin[ uiSmallestValid ] )
                        uiSmallestValid = uiI;
            }
        }

        void operator++( )
        {
            if( uiSmallestValid < vBegin.size( ) )
            {
                for( size_t uiI = 0; uiI < vBegin.size( ); uiI++ )
                    if( uiI != uiSmallestValid && *vBegin[ uiI ] != vEnd[ uiI ] &&
                        !( **vBegin[ uiSmallestValid ] < **vBegin[ uiI ] ) )
                        ++*vBegin[ uiI ];
                ++*vBegin[ uiSmallestValid ];
            }
            setSmallestValid( );
        }

        const coordinate_t operator*( ) const
        {
            return **vBegin[ uiSmallestValid ];
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
    sps_t uiBottomLeftPrefixSum;

    Overlay( ) : vSparseCoordsOverlay{ }, vSparseCoordsInternal{ }, vOverlayEntries{ }, xInternalEntires{ }
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

    template <size_t I, size_t N>
    inline void iterateDiagHelper( const size_t uiDiag, [[maybe_unused]] const size_t uiTotal,
                                   const std::array<coordinate_t, N>& rEnds,
                                   std::function<void( size_t, const std::array<coordinate_t, N>& )> fDo, size_t uiCurr,
                                   std::array<coordinate_t, N>& rCurr ) const
    {
        uiCurr *= rEnds[ I ];
        if constexpr( I >= N )
            throw std::runtime_error( "unreachable statement!" );
        if constexpr( I + 1 == N )
        {
            assert( uiTotal == rEnds[ I ] );
            assert( uiDiag < uiTotal );
            rCurr[ I ] = uiDiag;
            fDo( uiCurr + uiDiag, rCurr );
        }
        else
            for( coordinate_t uiI = uiDiag > uiTotal - rEnds[ I ] ? uiDiag - ( uiTotal - rEnds[ I ] ) : 0;
                 uiI < std::min( uiDiag + 1, rEnds[ I ] );
                 uiI++ )
            {
                rCurr[ I ] = uiI;
                iterateDiagHelper<I + 1, N>( uiDiag - uiI, uiTotal - ( rEnds[ I ] - 1 ), rEnds, fDo, uiCurr + uiI,
                                             rCurr );
            }
    }

    template <size_t N>
    void iterateDiag( const size_t uiDiag, const std::array<coordinate_t, N>& rEnds,
                      std::function<void( size_t, const std::array<coordinate_t, N>& )> fDo ) const
    {
        std::array<coordinate_t, N> rCurr;
        size_t uiCurr = 0;
        coordinate_t uiTotal = 1;
        for( coordinate_t uiX : rEnds )
            uiTotal += uiX - 1;
        iterateDiagHelper<0, N>( uiDiag, uiTotal, rEnds, fDo, uiCurr, rCurr );
    }

    coordinate_t getNumInternalSparseCoords( ) const
    {
        coordinate_t uiNumInternalSparseCoords = 0;
        for( size_t uiI = 0; uiI < D; uiI++ )
            if( vSparseCoordsInternal[ uiI ].uiStartIndex != std::numeric_limits<coordinate_t>::max( ) )
                uiNumInternalSparseCoords += sparse_coord_t::size( vSparseCoordsInternal[ uiI ] );

        return uiNumInternalSparseCoords;
    }

    coordinate_t getNumInternalPrefixSums( const sparse_coord_t& rSparseCoords ) const
    {
        pos_t vInternalAxisSizes = rSparseCoords.axisSizes( vSparseCoordsInternal );
        coordinate_t uiNumInternalPrefixSums = 1;
        for( size_t uiI = 0; uiI < D; uiI++ )
            uiNumInternalPrefixSums *= vInternalAxisSizes[ uiI ];

        return uiNumInternalPrefixSums;
    }

    coordinate_t generateInternalSparseCoords( sparse_coord_t& rSparseCoords, corners_t& vCorners,
                                               typename corners_t::Entry& xCorners,
                                               progress_stream_t&
#ifndef NDEBUG
                                                   xProg
#endif
    )
    {
        if( xCorners.size( ) > 0 )
        {
#ifndef NDEBUG
            xProg << Verbosity( 1 ) << "constructing sparse coordinates for points\n";
#endif
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
#ifndef NDEBUG
                xProg << Verbosity( 2 ) << "dim " << uiI << "\n";
#endif
                vCorners.sortByDim( uiI, xCorners );

#ifndef NDEBUG
                if( xProg.active( ) )
                    xCorners.stream( std::cout << "from points: ", vCorners ) << std::endl;
#endif

                vSparseCoordsInternal[ uiI ] = rSparseCoords.add( PointIterator( vCorners.cbegin( xCorners ), uiI ),
                                                                  PointIterator( vCorners.cend( xCorners ), uiI ) );
#ifndef NDEBUG
                if( xProg.active( ) )
                    vSparseCoordsInternal[ uiI ].stream( std::cout << "result: ", rSparseCoords ) << std::endl;
#endif
                // end of scope xFullLock
            }


            return getNumInternalPrefixSums( rSparseCoords );
        }
        return 0;
    }

    coordinate_t getNumOverlaySparseCoords( ) const
    {
        coordinate_t uiNumOverlaySparseCoords = 0;
        for( size_t uiI = 0; uiI < D; uiI++ )
            for( size_t uiJ = 0; uiJ < D - 1; uiJ++ )
                if( vSparseCoordsOverlay[ uiI ][ uiJ ].uiStartIndex != std::numeric_limits<coordinate_t>::max( ) )
                    uiNumOverlaySparseCoords += sparse_coord_t::size( vSparseCoordsOverlay[ uiI ][ uiJ ] );

        return uiNumOverlaySparseCoords;
    }

    coordinate_t getNumOverlayPrefixSums( const sparse_coord_t& rSparseCoords ) const
    {
        coordinate_t uiTotalOverlayPrefixSumSize = 0;
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            red_pos_t vAxisSizes = rSparseCoords.axisSizes( vSparseCoordsOverlay[ uiI ] );
            coordinate_t uiCurr = 1;
            for( size_t uiI = 0; uiI < D - 1; uiI++ )
                uiCurr *= vAxisSizes[ uiI ];
            uiTotalOverlayPrefixSumSize += uiCurr;
        }

        return uiTotalOverlayPrefixSumSize;
    }

    size_t allSubsetOfOne( std::vector<std::shared_ptr<MergableIterator>> vBegin,
                           std::vector<std::shared_ptr<MergableIterator>>
                               vEnd,
                           coordinate_t uiStart,
                           coordinate_t /*uiEnd*/ )
    {
        assert( vBegin.size( ) == vEnd.size( ) );

        std::vector<std::shared_ptr<MergableIterator>> vCurr;
        for( size_t uiI = 0; uiI < vBegin.size( ); uiI++ )
        {
            vCurr.push_back( vBegin[ uiI ]->copy( ) );
            while( *vCurr[ uiI ] != vEnd[ uiI ] && **vCurr[ uiI ] < uiStart )
                ++( *vCurr[ uiI ] );

            assert( !( *vCurr[ uiI ] != vEnd[ uiI ] ) || **vCurr[ uiI ] >= uiStart );
        }


        size_t uiLargest = 0;
        size_t uiLargestIdx = 0;
        for( size_t uiI = 0; uiI < vBegin.size( ); uiI++ )
        {
            size_t uiL = vEnd[ uiI ]->subtract( vCurr[ uiI ] );
            if( uiL > uiLargest )
            {
                uiLargest = uiL;
                uiLargestIdx = uiI;
            }
        }

        if( uiLargest == 0 )
            return vBegin.size( );


        while( *vCurr[ uiLargestIdx ] != vEnd[ uiLargestIdx ] )
        {
            coordinate_t uiCurr = **vCurr[ uiLargestIdx ];
            // assert( uiCurr < uiEnd );
            //  no iterator can be smaller than uiCurr
            //  otherwise we have found an element that is not in uiLargestIdx
            for( size_t uiI = 0; uiI < vBegin.size( ); uiI++ )
                if( *vCurr[ uiI ] != vEnd[ uiI ] && **vCurr[ uiI ] < uiCurr )
                    return vBegin.size( );

            // inc all iterators equal to uiCurr
            for( size_t uiI = 0; uiI < vBegin.size( ); uiI++ )
                while( *vCurr[ uiI ] != vEnd[ uiI ] && **vCurr[ uiI ] == uiCurr )
                    ++( *vCurr[ uiI ] );
        }

        // all iterators must have reached the end
        // otherwise we have found an element that is not in uiLargestIdx
        for( size_t uiI = 0; uiI < vBegin.size( ); uiI++ )
            if( *vCurr[ uiI ] != vEnd[ uiI ] )
                return vBegin.size( );

        return uiLargestIdx;
    }

    coordinate_t generateOverlaySparseCoords( const overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords,
                                              std::array<std::vector<coordinate_t>, D> vPredecessors,
                                              pos_t vMyBottomLeft, pos_t vPosTopRight,
                                              progress_stream_t&
#ifndef NDEBUG
                                                  xProg
#endif
    )
    {
#ifndef NDEBUG
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
                std::vector<typename sparse_coord_t::Entry> vPrev{ };
                std::vector<coordinate_t> vCollectedCoords;

                // add coordinates from previous overlays to the overlay entries
                for( size_t uiI2 = 0; uiI2 < D; uiI2++ )
                    if( uiJAct != uiI2 )
                        for( coordinate_t uiPred : vPredecessors[ uiI2 ] )
                        {
                            const Overlay* pPred = &rOverlays.vData[ uiPred ];
                            vPrev.push_back( pPred->vSparseCoordsOverlay[ uiI ][ uiJ ] );
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
                        }

                // add coordinates from the points of the previous overlay to the overlay entries
                for( coordinate_t uiPred : vPredecessors[ uiI ] ) // with uniform overlays this is only ever one!
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

                size_t uiSubsetOf = allSubsetOfOne( vBegin, vEnd, vMyBottomLeft[ uiJAct ], vPosTopRight[ uiJAct ] );
                if( uiSubsetOf < vPrev.size( ) )
                {
#ifndef NDEBUG
                    xProg << Verbosity( 2 ) << "using previous sparse coord lookup table " << uiSubsetOf << "\n";
#endif
                    // this exact set of sparse coords does exist
                    vSparseCoordsOverlay[ uiI ][ uiJ ] = vPrev[ uiSubsetOf ];
                }
                else
                { // this exact set of sparse coords does not yet exist
                    MergeIterator xBegin( vBegin, vEnd, vPosTopRight[ uiJAct ] );
                    MergeIterator xEnd( vEnd, vEnd, vPosTopRight[ uiJAct ] );

                    // skip over positions that are before and after this overlay
                    while( xBegin != xEnd && *xBegin < vMyBottomLeft[ uiJAct ] )
                        ++xBegin;

                    if( vPredecessors[ uiJAct ].size( ) > 0 )
                        vSparseCoordsOverlay[ uiI ][ uiJ ] =
                            rSparseCoords.addStartEnd( xBegin, xEnd, vMyBottomLeft[ uiJAct ] - 1 );
                    else
                        vSparseCoordsOverlay[ uiI ][ uiJ ] = rSparseCoords.add( xBegin, xEnd );
                }

#ifndef NDEBUG
                xProg << Verbosity( 2 );
                if( xProg.active( ) )
                    vSparseCoordsOverlay[ uiI ][ uiJ ].stream( std::cout << "result: ", rSparseCoords ) << std::endl;
#endif
            }
        }
        return getNumOverlayPrefixSums( rSparseCoords );
    }

    void generateInternalPrefixSums( const overlay_grid_t&, sparse_coord_t& rSparseCoords,
                                     prefix_sum_grid_t& rPrefixSums, corners_t& vCorners,
                                     typename corners_t::Entry& xCorners,
                                     progress_stream_t&
#ifndef NDEBUG
                                         xProg
#endif
    )
    {
        pos_t vInternalAxisSizes = rSparseCoords.axisSizes( vSparseCoordsInternal );

        coordinate_t uiNumTotal = 1;
        for( coordinate_t uiC : vInternalAxisSizes )
            uiNumTotal *= uiC;
        if( uiNumTotal > 0 )
        {
            std::vector<sps_t> vTmp( uiNumTotal, sps_t{ } );

#ifndef NDEBUG
            xProg << Verbosity( 1 ) << "constructing internal prefix sums: adding corner counts\n";
#endif
            vCorners.iterate(
                [ & ]( const corner_t& xCorner ) {
                    // compute index of prefix sum entry
                    pos_t uiSparsePos = rSparseCoords.sparse( xCorner.vPos, vSparseCoordsInternal );
                    coordinate_t uiIdx = 0;
                    for( size_t uiJ = 0; uiJ < D; uiJ++ )
                    {
                        assert( uiSparsePos[ uiJ ] != std::numeric_limits<coordinate_t>::max( ) );
                        assert( uiSparsePos[ uiJ ] < vInternalAxisSizes[ uiJ ] );
                        uiIdx = uiIdx * vInternalAxisSizes[ uiJ ] + uiSparsePos[ uiJ ];
                    }
                    assert( uiIdx < vTmp.size( ) );

                    // add point to entry
                    xCorner.addTo( vTmp[ uiIdx ] );
                },
                xCorners );

            pos_t vSizes;
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
                size_t uiJ = D - uiI - 1;
                vSizes[ uiJ ] = ( uiJ + 1 < D ? vSizes[ uiJ + 1 ] * vInternalAxisSizes[ uiJ + 1 ] : 1 );
            }
#ifndef NDEBUG
            xProg << Verbosity( 3 ) << "vSizes: " << vSizes << " vInternalAxisSizes: " << vInternalAxisSizes << "\n";
#endif

#if 1 // 1 == old implementation //@todo-low-prio new implementation is buggy!!! - keeping it off for now
      // compute internal prefix sum
            coordinate_t uiNumDone = 0;
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
#ifndef NDEBUG
                xProg << Verbosity( 3 ) << "computing prefix sums over dimension " << uiI << "\n";
#endif
                red_pos_t vRelevantInternalAxisSizes = relevant( vInternalAxisSizes, uiI );

                std::vector<red_pos_t> vvTos;
                std::vector<coordinate_t> vuiTos;
                iterate<D - 1>( vRelevantInternalAxisSizes, [ & ]( const red_pos_t& vTo ) {
                    pos_t vFullTo = expand( vTo, uiI );
                    sps_t uiPrefixSum = sps_t{ };

#ifndef NDEBUG
                    xProg << Verbosity( 3 ) << "starting...: " << vFullTo << ": " << uiPrefixSum << "\n";
#endif

                    iterate( vInternalAxisSizes[ uiI ], [ & ]( coordinate_t uiTo ) {
                        vFullTo[ uiI ] = uiTo;

                        // compute index of prefix sum entry
                        coordinate_t uiIdx = 0;
                        for( size_t uiJ = 0; uiJ < D; uiJ++ )
                        {
                            assert( vFullTo[ uiJ ] != std::numeric_limits<coordinate_t>::max( ) );
                            assert( vFullTo[ uiJ ] < vInternalAxisSizes[ uiJ ] );
                            uiIdx = uiIdx * vInternalAxisSizes[ uiJ ] + vFullTo[ uiJ ];
                        }
                        assert( uiIdx < vTmp.size( ) );

                        uiPrefixSum += vTmp[ uiIdx ];

#ifndef NDEBUG
                        xProg << Verbosity( 3 ) << vFullTo << ": " << uiPrefixSum << "\n";
#endif

                        vTmp[ uiIdx ] = uiPrefixSum;

                        ++uiNumDone;
                    } );
                } );
            }

#else
            coordinate_t uiMaxDiag = 0;
            for( coordinate_t uiC : vInternalAxisSizes )
                uiMaxDiag += uiC - 1;

            // compute internal prefix sum
            for( size_t uiI = 1; uiI <= uiMaxDiag; uiI++ )
            {
#ifndef NDEBUG
                xProg << Verbosity( 3 ) << "constructing internal prefix sums diagonal " << uiI << " of " << uiMaxDiag
                      << "\n";
#endif
                iterateDiag<D>( uiI, vInternalAxisSizes, [ & ]( const size_t uiIdx, const pos_t& vPos ) {
                    sps_t uiPrefixSum = sps_t{ };
                    for( size_t uiI = 0; uiI < D; uiI++ )
                        if( vPos[ uiI ] > 0 )
                        {
                            size_t uiJdx = uiIdx - vSizes[ uiI ];
#ifndef NDEBUG
                            xProg << Verbosity( 3 ) << "idx " << uiIdx << " + idx " << uiJdx << " uiI: " << uiI << "\n";
#endif
                            uiPrefixSum += vTmp[ uiJdx ];
                            for( size_t uiJ = 0; uiJ < uiI; uiJ++ )
                                if( vPos[ uiJ ] > 0 )
                                {
                                    uiPrefixSum -= vTmp[ uiJdx - vSizes[ uiJ ] ];
#ifndef NDEBUG
                                    xProg << Verbosity( 3 ) << "idx " << uiIdx << " - idx " << uiJdx - vSizes[ uiJ ]
                                          << " uiI: " << uiI << "\n";
#endif
                                }
                        }
                    vTmp[ uiIdx ] += uiPrefixSum;
                } );
            }
#endif
#ifndef NDEBUG
            xProg << Verbosity( 1 ) << "constructing internal prefix sums: transferring counts\n";
#endif
            // synchronization hidden within the function
            xInternalEntires = rPrefixSums.template add<D>( vInternalAxisSizes, vTmp );
        }
    }

    void generateOverlayPrefixSums( const overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords,
                                    prefix_sum_grid_t& rPrefixSums, corners_t&,
                                    std::array<std::vector<coordinate_t>, D> vPredecessors, pos_t vMyBottomLeft,
                                    dataset_t* pDataset, progress_stream_t& xProg )
    {
        // construct overlay sum grid
#ifndef NDEBUG
        xProg << Verbosity( 1 ) << "constructing overlay sum grid for overlay with vMyBottomLeft" << vMyBottomLeft
              << "\n";
#endif

        pos_t uiQPos;
        bool uiLargerZero = true;
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            uiQPos[ uiI ] = vMyBottomLeft[ uiI ] - 1;
            uiLargerZero = uiLargerZero && vMyBottomLeft[ uiI ] > 0;
        }
        if( uiLargerZero )
            uiBottomLeftPrefixSum = pDataset->getAll( rOverlays, rSparseCoords, rPrefixSums, uiQPos, xProg );
        else
            uiBottomLeftPrefixSum = sps_t{ };

        for( size_t uiI = 0; uiI < D; uiI++ )
            if( vPredecessors[ uiI ].size( ) > 0 )
            {
#ifndef NDEBUG
                xProg << Verbosity( 2 ) << "dim " << uiI << "\n";
#endif

                red_pos_t vAxisSizes = rSparseCoords.axisSizes( vSparseCoordsOverlay[ uiI ] );

                coordinate_t uiArea = 1;
                for( size_t uiJ = 0; uiJ < D - 1; uiJ++ )
                    uiArea *= vAxisSizes[ uiJ ];

                if( uiArea > 0 )
                {
#ifndef NDEBUG
                    for( size_t uiJ = 0; uiJ < D - 1; uiJ++ )
                        assert( vSparseCoordsOverlay[ uiI ][ uiJ ].uiStartIndex !=
                                std::numeric_limits<coordinate_t>::max( ) );
#endif

                    // no zero initialization needed -> every value will be overridden
                    vOverlayEntries[ uiI ] = rPrefixSums.template add<D - 1, false>( vAxisSizes );

                    assert( vMyBottomLeft[ uiI ] > 0 );
                    // if( prefix_sum_grid_t::sizeOf( vOverlayEntries[ uiI ] ) > 0 )
                    //{

                    rSparseCoords.template iterate<D - 1>(
                        [ & ]( const red_pos_t& vFrom, const red_pos_t& vTo ) {
                            pos_t vFullFrom = expand( vFrom, uiI );
                            vFullFrom[ uiI ] = vMyBottomLeft[ uiI ] - 1;

#ifndef NDEBUG
                            xProg << Verbosity( 3 ) << "query " << vFullFrom << "\n";
#endif
                            sps_t uiRet = pDataset->getAll( rOverlays, rSparseCoords, rPrefixSums, vFullFrom, xProg );
#ifndef NDEBUG
                            xProg << Verbosity( 3 ) << "query " << vFullFrom << ": " << uiRet << "\n";
#endif

                            rPrefixSums.get( vTo, vOverlayEntries[ uiI ] ) = uiRet;
                        },
                        vSparseCoordsOverlay[ uiI ] );
                    //}
                }
            }
    }

#ifdef UNROLL_FOR_ALL_COMBINATIONS
    template <size_t, size_t uiDistToTo, size_t uiFirstZero>
#endif
    struct CombinationsInvariant
    {
        static inline void count(

#ifndef UNROLL_FOR_ALL_COMBINATIONS
            size_t, size_t uiDistToTo, size_t uiFirstZero,
#endif
            pos_t vPos, val_t& uiRet, pos_t& /*vMyBottomLeft*/,
            const std::array<red_entry_arr_t, D>& vSparseCoordsOverlay, const sparse_coord_t& rSparseCoords,
            const prefix_sum_grid_t& rPrefixSums,
            const std::array<typename prefix_sum_grid_t::template Entry<D - 1>, D>& vOverlayEntries,
            [[maybe_unused]] size_t uiCornerIdx, sps_t uiBottomLeftPrefixSum
#if GET_PROG_PRINTS
            ,
            progress_stream_t& xProg
#endif
        )
        {
#if GET_PROG_PRINTS
            xProg << Verbosity( 3 ) << "\t\tONE: query: " << vPos << " in overlay " << uiFirstZero << "\n";
#endif

#ifdef UNROLL_FOR_ALL_COMBINATIONS
            if constexpr( uiDistToTo == 0 || uiFirstZero == D )
#else
            if( uiDistToTo == 0 || uiFirstZero == D )
#endif
                return;
            if( vOverlayEntries[ uiFirstZero ].uiStartIndex == std::numeric_limits<coordinate_t>::max( ) )
                return;

#if GET_PROG_PRINTS
            xProg << Verbosity( 3 ) << "\t\tpast return\n";
#endif
            red_pos_t vRelevant = relevant( vPos, uiFirstZero );
            red_pos_t vSparse =
                rSparseCoords.template sparse<D - 1, false>( vRelevant, vSparseCoordsOverlay[ uiFirstZero ] );

            bool bValid = true;
            for( size_t uiI = 0; uiI < D - 1; uiI++ )
                bValid = bValid && vSparse[ uiI ] != std::numeric_limits<coordinate_t>::max( );

#if GET_PROG_PRINTS
            xProg << "\t\trelevant: " << vRelevant << " sparse: " << vSparse
                  << " valid: " << ( bValid ? "true" : "false" ) << "\n";
#endif
            // in release mode query with sanity=false to avoid sanity checks
            sps_t uiCurrArr = rPrefixSums.template get<D - 1, false>( vSparse, vOverlayEntries[ uiFirstZero ] );

            val_t uiCurr;
            if( bValid )
            {
                if constexpr( IS_ORTHOTOPE )
                    uiCurr = uiCurrArr[ uiCornerIdx ];
                else
                    uiCurr = uiCurrArr;
            }
            else
            {
                if constexpr( IS_ORTHOTOPE )
                    uiCurr = uiBottomLeftPrefixSum[ uiCornerIdx ];
                else
                    uiCurr = uiBottomLeftPrefixSum;
            }

#if GET_PROG_PRINTS
            xProg << "\t\tis " << ( uiDistToTo % 2 == 0 ? "-" : "+" ) << uiCurr << " (" << uiCurrArr << ")\n";
#endif

#ifdef UNROLL_FOR_ALL_COMBINATIONS
            if constexpr( uiDistToTo % 2 == 0 )
#else
            if( uiDistToTo % 2 == 0 )
#endif
                uiRet -= uiCurr;
            else
                uiRet += uiCurr;
        }
    };

#ifdef UNROLL_FOR_ALL_COMBINATIONS
    template <size_t /*uiD*/, size_t uiDistToTo, size_t uiFirstZero>
#endif
    struct CombinationsInvariantAll
    {
        static inline void count(
#ifndef UNROLL_FOR_ALL_COMBINATIONS
            size_t, size_t uiDistToTo, size_t uiFirstZero,
#endif
            pos_t vPos, sps_t& uiRet, pos_t& /*vMyBottomLeft*/,
            const std::array<red_entry_arr_t, D>& vSparseCoordsOverlay, const sparse_coord_t& rSparseCoords,
            const prefix_sum_grid_t& rPrefixSums,
            const std::array<typename prefix_sum_grid_t::template Entry<D - 1>, D>& vOverlayEntries,
            sps_t uiBottomLeftPrefixSum,
            progress_stream_t&
#if GET_PROG_PRINTS
                xProg
#endif
        )
        {
#ifdef UNROLL_FOR_ALL_COMBINATIONS
            if constexpr( uiDistToTo == 0 || uiFirstZero == D )
#else
            if( uiDistToTo == 0 || uiFirstZero == D )
#endif
                return;
            if( vOverlayEntries[ uiFirstZero ].uiStartIndex == std::numeric_limits<coordinate_t>::max( ) )
                return;


#if GET_PROG_PRINTS
            xProg << Verbosity( 3 ) << "\t\tALL: query: " << vPos << " in overlay " << uiFirstZero << "\n";
#endif
            red_pos_t vRelevant = relevant( vPos, uiFirstZero );
            red_pos_t vSparse =
                rSparseCoords.template sparse<D - 1, false>( vRelevant, vSparseCoordsOverlay[ uiFirstZero ] );

            bool bValid = true;
            for( size_t uiI = 0; uiI < D - 1; uiI++ )
                bValid = bValid && vSparse[ uiI ] != std::numeric_limits<coordinate_t>::max( );

#if GET_PROG_PRINTS
            xProg << "\t\trelevant: " << vRelevant << " sparse: " << vSparse
                  << " valid: " << ( bValid ? "true" : "false" ) << "\n";
#endif
            sps_t uiCurr;
            if( bValid )
                uiCurr = rPrefixSums.template get<D - 1, false>( vSparse, vOverlayEntries[ uiFirstZero ] );
            else
                uiCurr = uiBottomLeftPrefixSum;

#if GET_PROG_PRINTS
            xProg << "\t\tis " << ( uiDistToTo % 2 == 0 ? "-" : "+" ) << uiCurr << "\n";
#endif

#ifndef NDEBUG
#if DU_UNREALISTIC_VALUE_CHECK
            if constexpr( IS_ORTHOTOPE )
            {
                for( size_t uiX = 0; uiX < 1 << ORTHOTOPE_DIMS; uiX++ )
                    if( uiCurr[ uiX ] >= std::numeric_limits<val_t>::max( ) / 2 )
                        throw std::runtime_error( "unrealistic value for uiCurr" );
            }
            else
            {
                if( uiCurr >= std::numeric_limits<val_t>::max( ) / 2 )
                    throw std::runtime_error( "unrealistic value for uiCurr" );
            }
#endif
#endif

#ifdef UNROLL_FOR_ALL_COMBINATIONS
            if constexpr( uiDistToTo % 2 == 0 )
#else
            if( uiDistToTo % 2 == 0 )
#endif
                uiRet -= uiCurr;
            else
                uiRet += uiCurr;
        }
    };

    static inline bool getCombinationsCond( coordinate_t uiPos, size_t /*uiD*/, bool /*bIsFrom*/ )
    {
        return uiPos != std::numeric_limits<coordinate_t>::max( );
    }


    val_t get( const sparse_coord_t& rSparseCoords, const prefix_sum_grid_t& rPrefixSums, pos_t vCoords,
               pos_t vMyBottomLeft, size_t uiCornerIdx
#if GET_PROG_PRINTS
               ,
               progress_stream_t& xProg
#endif
    ) const
    {
        val_t uiRet{ };

#if GET_PROG_PRINTS
        xProg << Verbosity( 2 ) << "\tquerying overlay for " << vCoords << " corner idx= " << uiCornerIdx
              << " bottom left prefix sum= " << uiBottomLeftPrefixSum << "...\n";
#endif

        for( size_t uiI = 0; uiI < D; uiI++ )
            --vMyBottomLeft[ uiI ]; // will turn zero values into max values
#if GET_PROG_PRINTS
        xProg << Verbosity( 3 ) << "\t\tvMyBottomLeft " << vMyBottomLeft << "\n";
#endif

        forAllCombinationsTmpl<pos_t, CombinationsInvariant>(
            vMyBottomLeft, vCoords, getCombinationsCond, uiRet, vMyBottomLeft, vSparseCoordsOverlay, rSparseCoords,
            rPrefixSums, vOverlayEntries, uiCornerIdx, uiBottomLeftPrefixSum
#if GET_PROG_PRINTS
            ,
            xProg
#endif

        );

        if( rPrefixSums.sizeOf( xInternalEntires ) > 0 )
        {
            auto vSparseCoords = rSparseCoords.template sparse<D, false>( vCoords, vSparseCoordsInternal );
            auto uiCurrArr = rPrefixSums.template get<D, false>( vSparseCoords, xInternalEntires );
            val_t uiCurr;
            if constexpr( IS_ORTHOTOPE )
                uiCurr = uiCurrArr[ uiCornerIdx ];
            else
                uiCurr = uiCurrArr;
#if GET_PROG_PRINTS
            xProg << Verbosity( 2 ) << "\tquerying internal " << vCoords << " -> " << vSparseCoords << ": +" << uiCurr
                  << " (" << uiCurrArr << ")\n";
#endif
#ifndef NDEBUG
#if DU_UNREALISTIC_VALUE_CHECK
            if( uiCurr >= std::numeric_limits<val_t>::max( ) / 2 )
                throw std::runtime_error( "unrealistic value for uiCurr" );
#endif
#endif

            uiRet += uiCurr;
        }
#if GET_PROG_PRINTS
        else
            xProg << Verbosity( 2 ) << "\tno internal entries\n";
#endif

#ifndef NDEBUG
#if DU_UNREALISTIC_VALUE_CHECK
        if( uiRet >= std::numeric_limits<val_t>::max( ) / 2 )
            throw std::runtime_error( "unrealistic value for uiRet" );
#endif
#endif

        return uiRet;
    }

#define CHECK_BIT( var, pos ) !!( ( var ) & ( 1 << ( pos ) ) )

    // @todo-low-prio this can be computed compiletime and then jsut put into a lookup array
    template <size_t uiD> static size_t intersectionTypeToCornerIndex( [[maybe_unused]] const isect_arr_t& vInterTypes )
    {
        if constexpr( IS_ORTHOTOPE )
        {
            size_t uiRet = 0;
            for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
            {
                size_t uiIsTop;

                switch( vInterTypes[ uiI ] )
                {
                    case IntersectionType::points_only:
                    case IntersectionType::first:
                        uiIsTop = 0;
                        break;
                    case IntersectionType::last:
                        uiIsTop = 1;
                        break;
                    default:
                    case IntersectionType::enclosed:
                    case IntersectionType::encloses:
                        uiIsTop = CHECK_BIT( uiD, D - uiI - 1 );
                        break;
                    case IntersectionType::overlaps:
                        uiIsTop = !CHECK_BIT( uiD, D - uiI - 1 );
                        break;
                }

                uiRet <<= 1;
                uiRet += uiIsTop;
            }
            return uiRet;
        }
        else
            return 0;
    }

    // @todo-low-prio code duplication
    static size_t intersectionTypeToCornerIndex( size_t uiD, [[maybe_unused]] const isect_arr_t& vInterTypes )
    {
        if constexpr( IS_ORTHOTOPE )
        {
            size_t uiRet = 0;
            for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
            {
                size_t uiIsTop;

                switch( vInterTypes[ uiI ] )
                {
                    case IntersectionType::points_only:
                    case IntersectionType::first:
                        uiIsTop = 0;
                        break;
                    case IntersectionType::last:
                        uiIsTop = 1;
                        break;
                    default:
                    case IntersectionType::enclosed:
                    case IntersectionType::encloses:
                        uiIsTop = CHECK_BIT( uiD, D - uiI - 1 );
                        break;
                    case IntersectionType::overlaps:
                        uiIsTop = !CHECK_BIT( uiD, D - uiI - 1 );
                        break;
                }

                uiRet <<= 1;
                uiRet += uiIsTop;
            }
            return uiRet;
        }
        else
            return 0;
    }

    static size_t intersectionTypeToFactor( [[maybe_unused]] const isect_arr_t& vInterTypes )
    {
        size_t uiFac = 1;
        if constexpr( IS_ORTHOTOPE )
            for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
                if( vInterTypes[ uiI ] == IntersectionType::encloses )
                    uiFac *= -1;
        return uiFac;
    }

    using grid_ret_t = NDGrid<type_defs, sps_t, RamVecGenerator>;
    using grid_ret_entry_t = typename grid_ret_t::template Entry<D>;

    template <size_t uiD>
    void addInternalValuesToGrid( grid_ret_t& rRet, pos_t& rSparsePos, pos_t& rGridIdx,
                                  const grid_ret_entry_t& xRetEntry, const pos_t& vGridFrom, const pos_t& vGridTo,
                                  const std::array<std::vector<coordinate_t>, D>& vvSparsePoss,
                                  const prefix_sum_grid_t& rPrefixSums
#if GET_PROG_PRINTS
                                  ,
                                  progress_stream_t& xProg
#endif
    ) const
    {
        if constexpr( uiD < D )
        {
            // iterate over grid-line intersections on the face of dimension uiO
            for( size_t uiX = 0; uiX < vGridTo[ uiD ] - vGridFrom[ uiD ]; uiX++ )
            {
                rGridIdx[ uiD ] = uiX + vGridFrom[ uiD ];
                rSparsePos[ uiD ] = vvSparsePoss[ uiD ][ uiX ];
                addInternalValuesToGrid<uiD + 1>( rRet, rSparsePos, rGridIdx, xRetEntry, vGridFrom, vGridTo,
                                                  vvSparsePoss, rPrefixSums
#if GET_PROG_PRINTS
                                                  ,
                                                  xProg
#endif
                );
            }
        }
        else
        {
#if GET_PROG_PRINTS
            xProg << Verbosity( 2 ) << "\t\tQuerying" //
                  << " sparse position " << rSparsePos //
                  << " grid pos " << rGridIdx //
                  << "\n";
#endif
            rRet.template get<D, false>( rGridIdx, xRetEntry ) =
                rPrefixSums.template get<D, false>( rSparsePos, xInternalEntires );
        }
    }

    class AddOverlayValuesToGrid
    {
        grid_ret_t& rRet;
        red_pos_t rSparsePos;
        pos_t rGridIdx;
        const std::array<grid_ret_entry_t, D>& vRetEntries;
        const sparse_coord_t& rSparseCoords;
        const prefix_sum_grid_t& rPrefixSums;
        const pos_t& vOverlayBottomLeft;
        const pos_t& vOverlayTopRight;
        const pos_t& vGridFrom;
        const pos_t& vGridTo;
        const std::array<std::vector<coordinate_t>, D>& vGrid;
        const pos_t& vOverlayIdxInGrid;
        const std::array<typename prefix_sum_grid_t::template Entry<D - 1>, D>& vOverlayEntries;
        const std::array<red_entry_arr_t, D>& vSparseCoordsOverlay;
        const sps_t& uiBottomLeftPrefixSum;
        std::array<std::vector<coordinate_t>, D>& vvSparsePoss;
#if GET_PROG_PRINTS
        progress_stream_t& xProg;
#endif

        template <size_t uiO, size_t uiD, bool bValid> inline void addValuesInDim( )
        {
            if constexpr( uiD == uiO )
            {
                rGridIdx[ uiD ] = vOverlayIdxInGrid[ uiD ];
                // for dimension uiO we will do the iteration once we have queried the overlay value
                addValuesInDim<uiO, uiD + 1, bValid>( );
            }
            else if constexpr( uiD < D )
                // iterate over grid-line intersections on the face of dimension uiO
                for( size_t uiX = 0; uiX < vGridTo[ uiD ] - vGridFrom[ uiD ]; uiX++ )
                {
                    rGridIdx[ uiD ] = uiX + vGridFrom[ uiD ];
                    if constexpr( bValid )
                    {
                        constexpr size_t uiDmO = uiD - ( uiD >= uiO ? 1 : 0 );
                        rSparsePos[ uiDmO ] = vvSparsePoss[ uiD ][ uiX ];
                        if( rSparsePos[ uiDmO ] != std::numeric_limits<coordinate_t>::max( ) )
                            addValuesInDim<uiO, uiD + 1, true>( );
                        else
                            addValuesInDim<uiO, uiD + 1, false>( );
                    }
                    else
                        addValuesInDim<uiO, uiD + 1, false>( );
                }
            else
            {
                sps_t uiValue;
                if constexpr( bValid )
                    uiValue = rPrefixSums.template get<D - 1, false>( rSparsePos, vOverlayEntries[ uiO ] );
                else
                    uiValue = uiBottomLeftPrefixSum;
#if GET_PROG_PRINTS
                xProg << Verbosity( 4 ) << "\t\t\tgrid_pos: " << rGridIdx << " sparse_pos: " << rSparsePos
                      << " value: " << uiValue << "\n";
#endif
                rRet.template get<D, false>( rGridIdx, vRetEntries[ uiO ] ) = uiValue;
            }
        }

        template <size_t uiD> inline void addValues( )
        {
            if constexpr( uiD < D )
            {
                if( rPrefixSums.sizeOf( vOverlayEntries[ uiD ] ) > 0 )
                {
                    // compute sparse coords for this overlay
                    for( size_t uiI = 0; uiI < D; uiI++ )
                    {
                        vvSparsePoss[ uiI ].clear( );
                        if( uiI != uiD )
                        {
                            // start of overlay needs to be considered
                            vvSparsePoss[ uiI ].push_back( rSparseCoords.template sparse<D - 1, false>(
                                vOverlayBottomLeft[ uiI - ( uiI > uiD ? 1 : 0 ) ], vSparseCoordsOverlay[ uiD ],
                                uiI - ( uiI > uiD ? 1 : 0 ) ) );
                            for( size_t uiX = vGridFrom[ uiI ] + 1; uiX < vGridTo[ uiI ]; uiX++ )
                                vvSparsePoss[ uiI ].push_back( rSparseCoords.template sparse<D - 1, false>(
                                    vGrid[ uiI ][ uiX + 1 ], vSparseCoordsOverlay[ uiD ],
                                    uiI - ( uiI > uiD ? 1 : 0 ) ) );
                        }
                    }

#if GET_PROG_PRINTS
                    xProg << Verbosity( 4 ) << "\t\taddValuesInDim " << uiD << " ret_grid_size: " << vRetEntries[ uiD ]
                          << " sparse_pos: " << vvSparsePoss << "\n";
#endif
                    addValuesInDim<uiD, 0, true>( );
                }
                else
                {
#if GET_PROG_PRINTS
                    xProg << Verbosity( 4 ) << "\t\tno overlay values in dim " << uiD << "\n";
#endif
                    addValuesInDim<uiD, 0, false>( );
                }
                addValues<uiD + 1>( );
            }
        }

      public:
        AddOverlayValuesToGrid( grid_ret_t& rRet, const std::array<grid_ret_entry_t, D>& vRetEntries,
                                const sparse_coord_t& rSparseCoords, const prefix_sum_grid_t& rPrefixSums,
                                const pos_t& vOverlayBottomLeft, const pos_t& vOverlayTopRight, const pos_t& vGridFrom,
                                const pos_t& vGridTo, const std::array<std::vector<coordinate_t>, D>& vGrid,
                                const pos_t& vOverlayIdxInGrid,
                                const std::array<typename prefix_sum_grid_t::template Entry<D - 1>, D>& vOverlayEntries,
                                const std::array<red_entry_arr_t, D>& vSparseCoordsOverlay,
                                const sps_t& uiBottomLeftPrefixSum,
                                std::array<std::vector<coordinate_t>, D>& vvSparsePoss
#if GET_PROG_PRINTS
                                ,
                                progress_stream_t& xProg
#endif
                                )
            : rRet( rRet ),
              vRetEntries( vRetEntries ),
              rSparseCoords( rSparseCoords ),
              rPrefixSums( rPrefixSums ),
              vOverlayBottomLeft( vOverlayBottomLeft ),
              vOverlayTopRight( vOverlayTopRight ),
              vGridFrom( vGridFrom ),
              vGridTo( vGridTo ),
              vGrid( vGrid ),
              vOverlayIdxInGrid( vOverlayIdxInGrid ),
              vOverlayEntries( vOverlayEntries ),
              vSparseCoordsOverlay( vSparseCoordsOverlay ),
              uiBottomLeftPrefixSum( uiBottomLeftPrefixSum ),
              vvSparsePoss( vvSparsePoss )
#if GET_PROG_PRINTS
              ,
              xProg( xProg )
#endif
        {
#if GET_PROG_PRINTS
            xProg << Verbosity( 4 ) << "\tAddOverlayValuesToGrid in overlay " << vOverlayIdxInGrid << " vGridFrom "
                  << vGridFrom << " vGridTo " << vGridTo << "\n";
#endif
            addValues<0>( );
        }
    };


    void grid( grid_ret_t& rRet, std::array<std::vector<coordinate_t>, D>& vvSparsePoss,
               const grid_ret_entry_t& xRetEntryInternal, const std::array<grid_ret_entry_t, D>& vRetEntriesOverlay,
               const sparse_coord_t& rSparseCoords, const prefix_sum_grid_t& rPrefixSums,
               const pos_t& vOverlayBottomLeft, const pos_t& vOverlayTopRight, const pos_t& vGridFrom,
               const pos_t& vGridTo, const pos_t& vGridFromOverlay, const pos_t& vGridToOverlay,
               const std::array<std::vector<coordinate_t>, D>& vGrid, const pos_t& vOverlayIdx
#if GET_PROG_PRINTS
               ,
               progress_stream_t& xProg
#endif
    ) const
    {
        // query overlay values and add them to the grid if necessary
        AddOverlayValuesToGrid( rRet, vRetEntriesOverlay, rSparseCoords, rPrefixSums, vOverlayBottomLeft,
                                vOverlayTopRight, vGridFromOverlay, vGridToOverlay, vGrid, vOverlayIdx, vOverlayEntries,
                                vSparseCoordsOverlay, uiBottomLeftPrefixSum, vvSparsePoss
#if GET_PROG_PRINTS
                                ,
                                xProg
#endif
        );

        const bool bHasInternal = rPrefixSums.sizeOf( xInternalEntires ) > 0;
        if( bHasInternal )
        {
            // compute the sparse coords of internal values
            for( size_t uiD = 0; uiD < D; uiD++ )
            {
                vvSparsePoss[ uiD ].clear( );
                for( size_t uiX = vGridFrom[ uiD ]; uiX < vGridTo[ uiD ]; uiX++ )
                    vvSparsePoss[ uiD ].push_back(
                        rSparseCoords.template sparse<D, false>( vGrid[ uiD ][ uiX ], vSparseCoordsInternal, uiD ) );
            }
#if GET_PROG_PRINTS
            xProg << Verbosity( 3 ) << "\t\tvvSparsePoss: " << vvSparsePoss << "\n";
#endif

            // LOOK UP VALUES FOR EACH CELL OF THE GRID (FOR INTERNAL VALUES)
            pos_t rSparsePos;
            pos_t rGridIdx;
            addInternalValuesToGrid<0>( rRet, rSparsePos, rGridIdx, xRetEntryInternal, vGridFrom, vGridTo, vvSparsePoss,
                                        rPrefixSums
#if GET_PROG_PRINTS
                                        ,
                                        xProg
#endif
            );
        }
#if GET_PROG_PRINTS
        else
            xProg << Verbosity( 3 ) << "\t\tNo internal values in this overlay.\n";
#endif
    }

    sps_t getAll( const sparse_coord_t& rSparseCoords, const prefix_sum_grid_t& rPrefixSums, pos_t vCoords,
                  pos_t vMyBottomLeft, progress_stream_t& xProg ) const
    {
        sps_t uiRet{ };

#if GET_PROG_PRINTS
        xProg << Verbosity( 2 ) << "\tquerying overlay for " << vCoords << "...\n";
#endif

        for( size_t uiI = 0; uiI < D; uiI++ )
            --vMyBottomLeft[ uiI ]; // will turn zero values into max values

#if GET_PROG_PRINTS
        xProg << Verbosity( 3 ) << "\t\tvMyBottomLeft " << vMyBottomLeft << "\n";
#endif
        forAllCombinationsTmpl<pos_t, CombinationsInvariantAll>(
            vMyBottomLeft, vCoords, getCombinationsCond, uiRet, vMyBottomLeft, vSparseCoordsOverlay, rSparseCoords,
            rPrefixSums, vOverlayEntries, uiBottomLeftPrefixSum, xProg

        );

        if( rPrefixSums.sizeOf( xInternalEntires ) > 0 )
        {
            auto vSparseCoords = rSparseCoords.template sparse<D, false>( vCoords, vSparseCoordsInternal );
            // in release mode query with sanity=false to avoid sanity checks
            sps_t uiCurr = rPrefixSums.template get<D, false>( vSparseCoords, xInternalEntires );
#if GET_PROG_PRINTS
            xProg << Verbosity( 2 ) << "\tquerying internal " << vCoords << " -> " << vSparseCoords << ": +" << uiCurr
                  << "\n";
#endif
#ifndef NDEBUG
#if DU_UNREALISTIC_VALUE_CHECK
            if constexpr( IS_ORTHOTOPE )
            {
                for( size_t uiX = 0; uiX < 1 << ORTHOTOPE_DIMS; uiX++ )
                    if( uiCurr[ uiX ] >= std::numeric_limits<val_t>::max( ) / 2 )
                        throw std::runtime_error( "unrealistic value for uiCurr" );
            }
            else
            {
                if( uiCurr >= std::numeric_limits<val_t>::max( ) / 2 )
                    throw std::runtime_error( "unrealistic value for uiCurr" );
            }
#endif
#endif

            uiRet += uiCurr;
        }
#if GET_PROG_PRINTS
        else
            xProg << Verbosity( 2 ) << "\tno internal entries\n";
#endif

#ifndef NDEBUG
#if DU_UNREALISTIC_VALUE_CHECK
        if constexpr( IS_ORTHOTOPE )
        {
            for( size_t uiX = 0; uiX < 1 << ORTHOTOPE_DIMS; uiX++ )
                if( uiRet[ uiX ] >= std::numeric_limits<val_t>::max( ) / 2 )
                    throw std::runtime_error( "unrealistic value for uiRet" );
        }
        else
        {
            if( uiRet >= std::numeric_limits<val_t>::max( ) / 2 )
                throw std::runtime_error( "unrealistic value for uiRet" );
        }
#endif
#endif


        return uiRet;
    }


    template <typename T> static inline std::array<T, D - 1> relevant( const std::array<T, D>& vAllEntries, size_t uiI )
    {
        std::array<T, D - 1> vRelevantEntries;
        for( size_t uiJ = 0; uiJ < D - 1; uiJ++ )
            vRelevantEntries[ uiJ ] = vAllEntries[ uiJ + ( uiJ >= uiI ? 1 : 0 ) ];
        return vRelevantEntries;
    }


    template <typename T> inline std::array<T, D> expand( const std::array<T, D - 1>& vCompressed, size_t uiI ) const
    {
        std::array<T, D> vAllEntries;
        for( size_t uiJ = 0; uiJ < uiI; uiJ++ )
            vAllEntries[ uiJ ] = vCompressed[ uiJ ];
        for( size_t uiJ = uiI + 1; uiJ < D; uiJ++ )
            vAllEntries[ uiJ ] = vCompressed[ uiJ - 1 ];
        vAllEntries[ uiI ] = 0; // init remaining value
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

        os << ">";

        return os;
    }

    std::ostream& stream( std::ostream& os, pos_t vGridPos, const sparse_coord_t& rSparseCoords,
                          const prefix_sum_grid_t& rPrefixSums, const dataset_t& rData ) const
    {
        os << "<";
        os << std::endl;
        auto vPosActual = rData.actualFromGridPos( vGridPos );
        os << "\tbottom left: " << vPosActual << "\n";
        auto vPosTR = rData.actualTopRightFromGridPos( vGridPos );
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

        os << "\tbottom left prefix sum: ";
        os << uiBottomLeftPrefixSum << std::endl;

        return os;
    }
    std::ostream& stream( std::ostream& os, pos_t vGridPos, const sparse_coord_t& rSparseCoords,
                          const prefix_sum_grid_t& rPrefixSums, const dataset_t& rData,
                          const corners_t& /*vCorners*/ ) const
    {
        this->stream( os, vGridPos, rSparseCoords, rPrefixSums, rData );

        os << ">";

        return os;
    }
};


} // namespace sps
