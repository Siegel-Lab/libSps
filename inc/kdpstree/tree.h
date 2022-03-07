#pragma once

#include "kdpstree/overlay_kd_tree.h"
#include "kdpstree/overlay_meta.h"
#include "kdpstree/point.h"
#include "kdpstree/points.h"
#include "kdpstree/type_defs.h"
#include <iostream>
#include <limits>


#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif

namespace kdpstree
{

template <typename type_defs> class Tree
{
    EXTRACT_TYPE_DEFS; // macro call

    using overlay_meta_t = OverlayMeta<type_defs>;
    using overlay_kd_tree_t = OverlayKdTree<type_defs>;
    using overlay_entries_t = OverlayEntries<type_defs>;
    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;
    using desc_t = Desc<type_defs>;

    vec_generator_t<std::pair<coordinate_t, val_t>> axis_vec_generator = vec_generator_t<std::pair<coordinate_t, val_t>>( );
    using axis_vec_t = typeof( axis_vec_generator( ) );

    vec_generator_t<coordinate_t> cord_vec_generator = vec_generator_t<coordinate_t>( );
    using cord_vec_t = typeof( cord_vec_generator( ) );

    vec_generator_t<val_t> prefix_sum_vec_generator = vec_generator_t<val_t>( );
    using prefix_sum_vec_t = typeof( prefix_sum_vec_generator( "" ) );


    overlay_kd_tree_t vTree;
    overlay_entries_t vEntries;
    points_t vPoints;
    desc_t vDesc;

    /**
     * @brief  see https://algorithmica.org/en/eytzinger
     *
     * @tparam axis_vec_t
     * @param vAxisVec
     * @param uiI
     * @param uiK
     * @param uiN
     * @return size_t
     */
    template <typename axis_vec_t>
    size_t constructEytzinger( axis_vec_t& vAxisVecs, size_t uiJ, size_t uiI, size_t uiK, size_t uiN )
    {
        if( uiK <= uiN )
        {
            uiI = constructEytzinger( vAxisVecs, uiJ, uiI, uiK * 2, uiN );
            vEntries.vData[ uiK - 1 + uiJ ] = vAxisVecs[ uiI++ ];
            uiI = constructEytzinger( vAxisVecs, uiJ, uiI, uiK * 2 + 1, uiN );
        }
        return uiI;
    }

    val_t countToZeroFromAllPrecursers(
                                    class_key_t xDatasetId, 
                                    std::array<coordinate_t, d> vOverlayBottomLeft,
                                    std::array<coordinate_t, d> vOverlayTopRight,
                                    std::array<coordinate_t, d>& vCurrStart,
                                    std::array<coordinate_t, d>& vCurrEnd,
                                    size_t uiD,
                                    size_t uiDSkip
                                  )
    {
        if(uiD == uiDSkip)
            uiD++;
        if(uiD < d)
        {
            val_t uiRet = 0;
            vCurrStart[uiD] = 0;
            vCurrEnd[uiD] = vOverlayBottomLeft[uiD];
            uiRet += countToZeroFromAllPrecursers(xDatasetId, vOverlayBottomLeft, vOverlayTopRight, 
                                                  vCurrStart, vCurrEnd, uiD+1, uiDSkip);
            vCurrStart[uiD] = vOverlayBottomLeft[uiD];
            vCurrEnd[uiD] = vOverlayTopRight[uiD];
            uiRet += countToZeroFromAllPrecursers(xDatasetId, vOverlayBottomLeft, vOverlayTopRight, 
                                                  vCurrStart, vCurrEnd, uiD+1, uiDSkip);
            return uiRet;
        }
        else
        {
            vCurrStart[uiDSkip] = vOverlayBottomLeft[uiDSkip];
            vCurrEnd[uiDSkip] = vOverlayTopRight[uiDSkip];
            
            val_t uiRet = count<false>( xDatasetId, vCurrStart, vCurrEnd );
            if constexpr( EXPLAIN_QUERY )
            {
                std::cerr << "\t\tcountToZeroFromAllPrecursers from=(";
                for(size_t uiI = 0; uiI < d; uiI++)
                    std::cerr << vCurrStart[uiI] << ", ";
                std::cerr << ") to=(";
                for(size_t uiI = 0; uiI < d; uiI++)
                    std::cerr << vCurrEnd[uiI] << ", ";
                std::cerr << ") val=" << uiRet << std::endl;
            }
            return uiRet;
        }
    }

    void generateForPointsHelper( class_key_t xDatasetId, std::function<void( offset_t, bool )> fRegisterMeInTree,
                                  size_t uiMaxPointsPerOverlay, size_t uiFrom, size_t uiTo,
                                  std::array<coordinate_t, d> vOverlayBottomLeft,
                                  std::array<coordinate_t, d> vOverlayTopRight,
                                  const std::array<cord_vec_t, d>& vAxisCoordinates
                                   )
    {
        if( uiFrom == uiTo )
            return;
        bool bAllPointsSame = true;
        vPoints.forRange(
            [ & ]( const point_t& rP1 ) {
                vPoints.forRange(
                    [ & ]( const point_t& rP2 ) {
                        for( size_t uiA = 0; uiA < b; uiA++ )
                            if( rP1.vPos[ uiA ] != rP2.vPos[ uiA ] )
                                bAllPointsSame = false;
                        return bAllPointsSame;
                    },
                    uiFrom, uiTo );
                return false;
            },
            uiFrom, uiTo );
        if( uiTo - uiFrom <= uiMaxPointsPerOverlay || bAllPointsSame )
        {
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "generating overlay" << std::endl;
            // recursion termination
            std::array<size_t, d> vBegins{ };
            std::array<size_t, d> vSizes{ };

            // for each dimension
            for( size_t uiI = 0; uiI < d; uiI++ )
            {
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\tgenerating axis for dimension=" << uiI << std::endl;
                // 1) find axis coordinates
                size_t uiAxisCoordinatesStart = 0;
                while(uiAxisCoordinatesStart < vAxisCoordinates[uiI].size() && 
                      vAxisCoordinates[uiI][uiAxisCoordinatesStart] < vOverlayBottomLeft[uiI])
                    ++uiAxisCoordinatesStart;
                size_t uiAxisCoordinatesEnd = uiAxisCoordinatesStart;
                // @todo: continue here -> somehow the end coordinate is messed up...
                while(uiAxisCoordinatesEnd < vAxisCoordinates[uiI].size() && 
                      vAxisCoordinates[uiI][uiAxisCoordinatesEnd] <= vOverlayTopRight[uiI])
                    ++uiAxisCoordinatesEnd;

                // 2) fill a temp vec with the values
                auto vTmp = axis_vec_generator( );
                vPoints.sortByDim( uiI, uiFrom, uiTo );
                size_t uiPointsIt = uiFrom;
                
                pos_t vPos = vOverlayTopRight;
                for( size_t uiA = 0; uiA < b; uiA++ )
                    if(vPos[uiA] < std::numeric_limits<coordinate_t>::max())
                        ++vPos[uiA];
                vPos[ uiI ] = vOverlayBottomLeft[ uiI ];
                if constexpr( EXPLAIN_QUERY )
                {
                    std::cerr << "\toverlay init val pos=(";
                    for( size_t uiA = 0; uiA < b; uiA++ )
                        std::cerr << vPos[uiA] << ", ";
                    std::cerr << ")" << std::endl;
                }
    
                val_t uiInitVal = countToZero<false>( xDatasetId, vPos );

                for(size_t uiJ = uiAxisCoordinatesStart; uiJ < uiAxisCoordinatesEnd; uiJ++)
                {
                    // set point iterator to correct position
                    while(uiPointsIt < uiTo && vPoints.vData[ uiPointsIt ].vPos[ uiI ] <= vAxisCoordinates[uiI][uiJ])
                        ++uiPointsIt;

                    val_t uiVal = uiInitVal;
                    
                    if constexpr( EXPLAIN_QUERY )
                        std::cerr << "\toverlay axis pos=" << vAxisCoordinates[uiI][uiJ] << " init val=" 
                                    << uiVal << std::endl;

                    std::array<coordinate_t, d> vStart{};
                    std::array<coordinate_t, d> vEnd{};
                    pos_t vPos = vOverlayTopRight;
                    vPos[uiI] = vAxisCoordinates[uiI][uiJ];
                    for( size_t uiA = 0; uiA < b; uiA++ )
                        if(vPos[uiA] < std::numeric_limits<coordinate_t>::max())
                            ++vPos[uiA];
                    if constexpr( EXPLAIN_QUERY )
                    {
                        std::cerr << "\toverlay axis pos=" << vAxisCoordinates[uiI][uiJ] << " precurser pos=(";
                        for( size_t uiA = 0; uiA < b; uiA++ )
                            std::cerr << vPos[uiA] << ", ";
                        std::cerr << ")" << std::endl;
                    }

                    uiVal += countToZeroFromAllPrecursers(xDatasetId, vOverlayBottomLeft, vPos,
                                                          vStart, vEnd, 0, uiI);

                    if constexpr( EXPLAIN_QUERY )
                        std::cerr << "\toverlay axis pos=" << vAxisCoordinates[uiI][uiJ] << " init val + precursers =" 
                                    << uiVal << std::endl;

                    uiVal += uiPointsIt - uiFrom;

                    if constexpr( EXPLAIN_QUERY )
                        std::cerr << "\toverlay axis pos=" << vAxisCoordinates[uiI][uiJ] 
                                  << " init val + precursers + points =" 
                                    << uiVal << std::endl;

                    if( ( vTmp.size() == 0 || vTmp.back().second < uiVal ) && uiVal > 0)
                        vTmp.emplace_back( vAxisCoordinates[uiI][uiJ], uiVal );
                }

                assert(vTmp.size() > 0);

                if constexpr( EXPLAIN_QUERY )
                    for( size_t uiX = 0; uiX < vTmp.size( ); uiX++ )
                        std::cerr << "\tfinal: overlay axis pos=" << vTmp[ uiX ].first << " val=" << vTmp[ uiX ].second
                                  << std::endl;

                // 3) save pointers to the constructed axes
                vBegins[ uiI ] = vEntries.size( );
                vSizes[ uiI ] = vTmp.size( );
                vEntries.incSize( vTmp.size( ) );

                // 4) construct eytzinger representation of sorted compressed axis
                constructEytzinger( vTmp, vBegins[ uiI ], 0, 1, vTmp.size( ) );
            } // for

            // std::cout << " insert ";
            // for( size_t uiI = 0; uiI < d; uiI++ )
            //    std::cout << vOverlayPos[ uiI ] << ", ";
            // std::cout << std::endl;
            // insert overlay into correct position in tree

            offset_t uiMyOffset = vTree.vLeaves.size( );
            vTree.vLeaves.emplace_back( vBegins, vSizes, uiFrom, uiTo );

            fRegisterMeInTree( uiMyOffset, true );
        }
        else
        {
            // recursive call: split points in half

            // 1) find dimension that splits points most evenly in half
            size_t uiBestDim = 0;
            size_t uiBestCount = std::numeric_limits<size_t>::max( );
            for( size_t uiI = 0; uiI < d; uiI++ )
            {
                vPoints.sortByDim( uiI, uiFrom, uiTo );
                std::array<coordinate_t, b> uiAxisSplitPos;
                for( size_t uiA = 0; uiA < b; uiA++ )
                    uiAxisSplitPos[ uiA ] = vPoints.vData[ ( uiA * ( uiTo - uiFrom + 1 ) ) / b + uiFrom ].vPos[ uiI ];
                std::array<size_t, b> uiCountInSplit{ };
                vPoints.forRange(
                    [ & ]( const point_t& rP ) {
                        for( size_t uiA = 1; uiA <= b; uiA++ )
                            if( rP.vPos[ uiI ] >= uiAxisSplitPos[ b - uiA ] )
                            {
                                uiCountInSplit[ b - uiA ] += 1;
                                break;
                            }
                        return true;
                    },
                    uiFrom, uiTo );

                size_t uiMaxCount = 0;
                for( size_t uiA = 0; uiA < b; uiA++ )
                    uiMaxCount = std::max( uiMaxCount, uiCountInSplit[ uiA ] );

                if( uiMaxCount < uiBestCount )
                {
                    uiBestDim = uiI;
                    uiBestCount = uiMaxCount;
                }
            }

            // 2) seperate points in that dimension (by sorting and looking for the split pos)
            vPoints.sortByDim( uiBestDim, uiFrom, uiTo );
            std::array<coordinate_t, b> uiAxisSplitPos;
            for( size_t uiA = 0; uiA < b; uiA++ )
                uiAxisSplitPos[ uiA ] = vPoints.vData[ ( uiA * ( uiTo - uiFrom + 1 ) ) / b + uiFrom ].vPos[ uiBestDim ];
            std::array<size_t, b> uiCountInSplit{ };
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "splitting in dimension " << uiBestDim << " start "
                          << vPoints.vData[ uiFrom ].vPos[ uiBestDim ] << " end "
                          << vPoints.vData[ uiTo - 1 ].vPos[ uiBestDim ] << " with max count " << uiBestCount
                          << std::endl;
            vPoints.forRange(
                [ & ]( const point_t& rP ) {
                    for( size_t uiA = 1; uiA <= b; uiA++ )
                        if( rP.vPos[ uiBestDim ] >= uiAxisSplitPos[ b - uiA ] )
                        {
                            if constexpr( EXPLAIN_QUERY )
                            {
                                std::cerr << "point (";
                                for( size_t uiI = 0; uiI < d; uiI++ )
                                    std::cerr << rP.vPos[ uiI ] << ", ";
                                std::cerr << ")" << std::endl;
                            }
                            uiCountInSplit[ b - uiA ] += 1;
                            break;
                        }
                    return true;
                },
                uiFrom, uiTo );

            if constexpr( EXPLAIN_QUERY )
                for( size_t uiA = 0; uiA < b; uiA++ )
                    std::cerr << "\tat coordinate " << uiAxisSplitPos[ uiA ] << " containing " << uiCountInSplit[ uiA ]
                              << " points" << std::endl;

            // 3) insert new kd-tree node

            offset_t uiMyOffset = vTree.vTree.size( );
            vTree.vTree.emplace_back( );
            vTree.vTree[ uiMyOffset ].uiSplitDimension = uiBestDim;
            for( size_t uiA = 0; uiA < b; uiA++ )
            {
                std::get<0>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = std::numeric_limits<coordinate_t>::max( );
                std::get<1>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = std::numeric_limits<offset_t>::max( );
            }

            bool bRegistered = false;

            for( size_t uiA = 0; uiA < b; uiA++ )
            {

                vOverlayBottomLeft[ uiBestDim ] = uiAxisSplitPos[ uiA ];
                if( uiA + 1 < b )
                    vOverlayTopRight[ uiBestDim ] = uiAxisSplitPos[ uiA + 1 ];

                // std::cout << "split dim " << uiBestDim << " at " << uiSplitPos << " into " << 100 * fBestRatio << "%"
                //          << std::endl;


                // 4) recursive calls
                // here it is important to go from left to right so that vvPrefixSumVecs is filled in the correct order
                //


                generateForPointsHelper(
                    xDatasetId,
                    [ & ]( size_t uiOffsetNode, bool bIsLeaf ) {
                        std::get<0>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = uiAxisSplitPos[ uiA ];
                        std::get<1>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = uiOffsetNode;
                        std::get<2>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = bIsLeaf;

                        if(!bRegistered)
                            fRegisterMeInTree( uiMyOffset, false );
                        bRegistered = true;
                    },
                    uiMaxPointsPerOverlay, uiFrom, uiFrom + uiCountInSplit[ uiA ], vOverlayBottomLeft,
                    vOverlayTopRight, vAxisCoordinates );

                uiFrom += uiCountInSplit[ uiA ];
            }
            if(!bRegistered)
                fRegisterMeInTree( uiMyOffset, false );
            assert( uiFrom == uiTo );
        }
    }

    template<bool SILENT>
    std::pair<pos_t, const overlay_meta_t*> getOverlay( const class_key_t& xDatasetId, const pos_t& vPos ) const
    {
        for( size_t uiI = 0; uiI < vTree.vRoots.size( ); uiI++ )
        {
            if( std::get<0>( vTree.vRoots[ uiI ] ) == xDatasetId )
            {
                pos_t vBottomLeft{ };
                size_t uiCurr = std::get<1>( vTree.vRoots[ uiI ] );
                bool bIsLeaf = std::get<2>( vTree.vRoots[ uiI ] );
                while( !bIsLeaf )
                {
                    assert( uiCurr < vTree.vTree.size( ) );
                    bool bSet = false;
                    for( size_t uiA = 1; uiA <= b; uiA++ )
                        if( std::get<0>( vTree.vTree[ uiCurr ].vChildren[ b - uiA ] ) <=
                            vPos[ vTree.vTree[ uiCurr ].uiSplitDimension ] && 
                            std::get<1>( vTree.vTree[ uiCurr ].vChildren[ b - uiA ] ) != std::numeric_limits<offset_t>::max( ))
                        {
                            vBottomLeft[ vTree.vTree[ uiCurr ].uiSplitDimension ] =
                                std::get<0>( vTree.vTree[ uiCurr ].vChildren[ b - uiA ] );
                            bIsLeaf = std::get<2>( vTree.vTree[ uiCurr ].vChildren[ b - uiA ] );
                            uiCurr = std::get<1>( vTree.vTree[ uiCurr ].vChildren[ b - uiA ] );
                            bSet = true;
                            break;
                        }
                    if( !bSet )
                        return std::make_pair( vBottomLeft, nullptr );
                }
                assert( uiCurr < vTree.vLeaves.size( ) );
                if constexpr( EXPLAIN_QUERY && !SILENT )
                    std::cerr << "\t\toverlay index: " << uiCurr << std::endl;
                return std::make_pair( vBottomLeft, &vTree.vLeaves[ uiCurr ] );
            }
        }
        throw std::runtime_error( "getOverlay: dataset Id not found" );
    }

    void iterateOverlaysIn( std::function<void( const overlay_meta_t& )> fDo, size_t uiCurr, bool bIsLeaf,
                            const pos_t& vFrom, const pos_t& vTo ) const
    {
        if( bIsLeaf )
            fDo( vTree.vLeaves[ uiCurr ] );
        else
            for( size_t uiA = 0; uiA < b; )
                if( std::get<1>( vTree.vTree[ uiCurr ].vChildren[ uiA ] ) != std::numeric_limits<offset_t>::max( ) )
                {
                    // find next valid pos
                    size_t uiB = uiA + 1;
                    for( ; uiB < b && std::get<1>( vTree.vTree[ uiCurr ].vChildren[ uiB ] ) ==
                                          std::numeric_limits<offset_t>::max( );
                         uiB++ )
                        ;

                    bool bLeftBorderOk = vTo[ vTree.vTree[ uiCurr ].uiSplitDimension ] >
                                         std::get<0>( vTree.vTree[ uiCurr ].vChildren[ uiA ] );
                    bool bRightBorderOk = uiB == b || vFrom[ vTree.vTree[ uiCurr ].uiSplitDimension ] <=
                                                          std::get<0>( vTree.vTree[ uiCurr ].vChildren[ uiB ] );

                    if( bLeftBorderOk && bRightBorderOk )
                        iterateOverlaysIn( fDo, std::get<1>( vTree.vTree[ uiCurr ].vChildren[ uiA ] ),
                                           std::get<2>( vTree.vTree[ uiCurr ].vChildren[ uiA ] ), vFrom, vTo );
                    uiA = uiB;
                }
                else
                    ++uiA;
    }

    void iterateOverlaysIn( std::function<void( const overlay_meta_t& )> fDo, const class_key_t& xDatasetId,
                            const pos_t& vFrom, const pos_t& vTo ) const
    {
        for( size_t uiI = 0; uiI < d; uiI++ )
            if( vFrom[ uiI ] >= vTo[ uiI ] )
                return;
        for( size_t uiI = 0; uiI < vTree.vRoots.size( ); uiI++ )
            if( std::get<0>( vTree.vRoots[ uiI ] ) == xDatasetId )
                iterateOverlaysIn( fDo, std::get<1>( vTree.vRoots[ uiI ] ), std::get<2>( vTree.vRoots[ uiI ] ), vFrom,
                                   vTo );
    }


    template<bool SILENT>
    val_t countToZero( const overlay_meta_t& rO, pos_t vPos ) const
    {
        val_t uiRet = 0;
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            if constexpr( EXPLAIN_QUERY && !SILENT )
                std::cerr << "\t\toverlay value: dimension=" << uiI << " position=" << vPos[ uiI ] << std::endl;
            val_t uiVal = vEntries.get( vPos[ uiI ], rO.vEntryBegins[ uiI ], rO.vSizes[ uiI ] );
            if constexpr( EXPLAIN_QUERY && !SILENT )
                std::cerr << "\t\toverlay value=" << uiVal << std::endl;
            uiRet += uiVal;
        }

        val_t uiVal = 0;
        vPoints.forRange(
            [ & ]( const point_t& rP ) {
                bool bToTopRight = true;
                for( size_t uiI = 0; uiI < d; uiI++ )
                    bToTopRight = bToTopRight && ( rP.vPos[ uiI ] > vPos[ uiI ] );
                if( bToTopRight )
                    uiVal += 1;
                return true;
            },
            rO.uiPointsBegin, rO.uiPointsEnd );
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\toverlay number of points: " << ( d - 1 ) << " * " << uiVal
                      << " of total: " << rO.uiPointsEnd - rO.uiPointsBegin << std::endl;
        uiRet += ( d - 1 ) * uiVal;

        uiVal = vEntries.getLast( rO.vEntryBegins[ 0 ], rO.vSizes[ 0 ] );
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\toverlay top right corner: value= -" << ( d - 1 ) << " * " << uiVal << std::endl;
        assert( uiRet >= ( d - 1 ) * uiVal ); // can never have a negative return number
        uiRet -= ( d - 1 ) * uiVal;
        return uiRet;
    }

    template<bool SILENT>
    val_t countToZero( class_key_t xDatasetId, pos_t vPos ) const
    {
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\tcountToZero" << std::endl;
        if( vTree.vRoots.size( ) == 0 )
            return 0;
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            if( vPos[ uiI ] == 0 )
            {
                if constexpr( EXPLAIN_QUERY && !SILENT )
                    std::cerr << "\t\tdimension: " << uiI << " is zero. result must be zero as well." << std::endl;
                return 0;
            }
            else
                --vPos[ uiI ];
        }
        std::pair<pos_t, const overlay_meta_t*> rO = getOverlay<SILENT>( xDatasetId, vPos );
        if( rO.second == nullptr )
        {
            if constexpr( EXPLAIN_QUERY && !SILENT )
                std::cerr << "\t\tpoint is to the bottom left of most bottom left overlay. result must be zero."
                          << std::endl;
            return 0;
        }
        if constexpr( EXPLAIN_QUERY && !SILENT )
        {
            std::cerr << "\t\toverlay position: (";
            for( size_t uiI = 0; uiI < d; uiI++ )
                std::cerr << rO.first[ uiI ] << ", ";
            std::cerr << ")" << std::endl;
        }
        return countToZero<SILENT>( *rO.second, vPos );
    }

    template<bool SILENT>
    val_t countHelper( const class_key_t& xDatasetId, const pos_t& vFrom, const pos_t& vTo, pos_t& vCurr,
                       size_t uiDCurr, size_t uiNumStart ) const
    {
        if( uiDCurr == d )
        {
            if constexpr( EXPLAIN_QUERY && !SILENT )
            {
                std::cerr << "\tcounting from (";
                for( size_t uiI = 0; uiI < d; uiI++ )
                    std::cerr << "0, ";
                std::cerr << ") to (";
                for( size_t uiI = 0; uiI < d; uiI++ )
                    std::cerr << vCurr[ uiI ] << ", ";
                std::cerr << ") ";
                if( uiNumStart % 2 == 0 )
                    std::cerr << "positive";
                else
                    std::cerr << "negative";
                std::cerr << std::endl;
            }
            auto uiVal = countToZero<SILENT>( xDatasetId, vCurr );
            if constexpr( EXPLAIN_QUERY && !SILENT )
                std::cerr << "\tresult = " << uiVal << std::endl;
            return uiVal * ( uiNumStart % 2 == 0 ? 1 : -1 );
        }
        else
        {
            vCurr[ uiDCurr ] = vFrom[ uiDCurr ];
            val_t uiRet = countHelper<SILENT>( xDatasetId, vFrom, vTo, vCurr, uiDCurr + 1, uiNumStart + 1 );
            vCurr[ uiDCurr ] = vTo[ uiDCurr ];
            uiRet += countHelper<SILENT>( xDatasetId, vFrom, vTo, vCurr, uiDCurr + 1, uiNumStart );
            return uiRet;
        }
    }

  public:
    Tree( std::string sPrefix ) : vTree( sPrefix ), vEntries( sPrefix ), vPoints( sPrefix ), vDesc( sPrefix )
    {}

    void clear( )
    {
        vTree.clear( );
        vEntries.clear( );
        vPoints.clear( );
        vDesc.clear( );
    }

    size_t numPoints( ) const
    {
        return vPoints.size( );
    }

    size_t addPoint( pos_t vPos, std::string sDesc )
    {
        vPoints.add( vPos, vDesc.add( sDesc ) );
    }

    void generateForPoints( class_key_t xDatasetId, size_t uiMaxPointsPerOverlay, size_t uiFrom, size_t uiTo )
    {
        std::array<coordinate_t, d> vOverlayBottomLeft{ };
        std::array<coordinate_t, d> vOverlayTopRight{ };
        std::array<cord_vec_t, d> vAxisCoordinates{ };
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            vOverlayTopRight[ uiI ] = std::numeric_limits<coordinate_t>::max( );

            vAxisCoordinates[ uiI ] = cord_vec_generator( );
            vPoints.sortByDim(uiI, uiFrom, uiTo);
            vPoints.forRange(
                [ & ]( const point_t& rP2 ) {
                    if( vAxisCoordinates[ uiI ].size() == 0 || vAxisCoordinates[ uiI ].back() < rP2.vPos[ uiI ] )
                        vAxisCoordinates[ uiI ].push_back( rP2.vPos[ uiI ] );
                    return true;
                },
                uiFrom, uiTo );
        }
        generateForPointsHelper(
            xDatasetId, //
            [ & ]( size_t uiOffsetNode, bool bIsLeaf ) {
                vTree.vRoots.emplace_back( xDatasetId, uiOffsetNode, bIsLeaf );
            },
            uiMaxPointsPerOverlay, //
            uiFrom, //
            uiTo, //
            vOverlayBottomLeft, //
            vOverlayTopRight, //
            vAxisCoordinates );
    }

    template<bool SILENT>
    val_t count( class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        pos_t vCurr{ };
        if constexpr( EXPLAIN_QUERY && !SILENT )
        {
            std::cerr << "counting query for datsetId=" << xDatasetId << " bottomLeftCorner=(";
            for( size_t uiI = 0; uiI < d; uiI++ )
                std::cerr << vFrom[ uiI ] << ", ";
            std::cerr << ") topRightCorner=(";
            for( size_t uiI = 0; uiI < d; uiI++ )
                std::cerr << vTo[ uiI ] << ", ";
            std::cerr << ")" << std::endl;
        }
        return countHelper<SILENT>( xDatasetId, vFrom, vTo, vCurr, 0, 0 );
    }

    void iterate( std::function<void( pos_t, std::string )> fDo, class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        iterateOverlaysIn(
            [ & ]( const overlay_meta_t& xO ) {
                vPoints.forRange(
                    [ & ]( const point_t& xPoint ) {
                        for( size_t uiI = 0; uiI < d; uiI++ )
                        {
                            if( xPoint.vPos[ uiI ] >= vTo[ uiI ] )
                                return true;
                            if( xPoint.vPos[ uiI ] < vFrom[ uiI ] )
                                return true;
                        }
                        fDo( xPoint.vPos, vDesc.get( xPoint.uiDescOffset ) );
                        return true;
                    },
                    xO.uiPointsBegin, xO.uiPointsEnd );
            },
            xDatasetId, vFrom, vTo );
    }

    std::vector<std::pair<pos_t, std::string>> get( class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        std::vector<std::pair<pos_t, std::string>> vRet;
        iterate( [ & ]( pos_t vPos, std::string sDesc ) { vRet.emplace_back( vPos, sDesc ); }, xDatasetId, vFrom, vTo );
        return vRet;
    }

    std::string print( ) const
    {
        return "Tree:\n" + vTree.print( ) + "\nEntries:\n" + vEntries.print( ) + "Points:\n" + vPoints.print( ) +
               "Desc:\n" + vDesc.print( );
    }
};

} // namespace kdpstree

#if WITH_PYTHON
template <typename type_defs> void exportTree( pybind11::module& m, std::string sName )
{
    pybind11::class_<kdpstree::Tree<type_defs>>( m, sName.c_str( ) )
        .def( pybind11::init<std::string>( ) ) // constructor
        .def( "add_point", &kdpstree::Tree<type_defs>::addPoint )
        .def( "num_points", &kdpstree::Tree<type_defs>::numPoints )
        .def( "generate_for_points", &kdpstree::Tree<type_defs>::generateForPoints )
        .def( "count", &kdpstree::Tree<type_defs>::template count<false>, "" )
        .def( "get", &kdpstree::Tree<type_defs>::get, "" )
        .def( "__str__", &kdpstree::Tree<type_defs>::print )
        .def( "clear", &kdpstree::Tree<type_defs>::clear )

        ;
}
#endif