#pragma once

#include "kdpstree/overlay_kd_tree.h"
#include "kdpstree/overlay_meta.h"
#include "kdpstree/point.h"
#include "kdpstree/points.h"
#include "kdpstree/type_defs.h"
#include "kdpstree/util.h"
#include <iostream>
#include <limits>


#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif

namespace kdpstree
{

template<typename T, size_t N>
std::array<T, N>& operator+=(std::array<T, N>& vA, const std::array<T, N>& vB)
{
    for(size_t uiI = 1; uiI < N; uiI++)
        vA[uiI] += vB[UiI];
    return vA;
}

template<typename T, size_t N>
std::array<T, N>& operator-=(std::array<T, N>& vA, const std::array<T, N>& vB)
{
    for(size_t uiI = 1; uiI < N; uiI++)
        vA[uiI] -= vB[UiI];
    return vA;
}

template<typename T, size_t N>
bool operator>(std::array<T, N>& vA, size_t uiAs)
{
    for(size_t uiI = 0; uiI < N; uiI++)
        if(!vA[uiI] < uiAs)
            return false;
    return true;
}

template <typename type_defs> class Tree
{
    EXTRACT_TYPE_DEFS; // macro call

    using overlay_meta_t = OverlayMeta<type_defs>;
    using overlay_kd_tree_t = OverlayKdTree<type_defs>;
    using overlay_entries_t = OverlayEntries<type_defs>;
    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;
    using desc_t = Desc<type_defs>;

    vec_generator_t<coordinate_t> cord_vec_generator = vec_generator_t<coordinate_t>( );
    using cord_vec_t = typeof( cord_vec_generator( ) );

    vec_generator_t<data_t> prefix_sum_vec_generator = vec_generator_t<data_t>( );
    using prefix_sum_vec_t = typeof( prefix_sum_vec_generator( "" ) );

    vec_generator_t<point_t> point_vec_generator = vec_generator_t<point_t>( );
    using point_vec_t = typeof( point_vec_generator( ) );


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

    void generateForPointsHelper( class_key_t xDatasetId, std::function<void( offset_t, bool )> fRegisterMeInTree,
                                  size_t uiMaxPointsPerOverlay, size_t uiFrom, size_t uiTo,
                                  std::array<coordinate_t, d> vOverlayBottomLeft,
                                  std::array<coordinate_t, d> vOverlayTopRight,
                                  const std::array<cord_vec_t, d>& vAxisCoordinates,
                                  std::array<size_t, d> vAxisFrom,
                                  std::array<size_t, d> vAxisTo,
                                  point_vec_t& vPointsEnd,
                                  size_t uiFromEnd,
                                  size_t uiToEnd
                                   )
    {
        if( uiFrom == uiTo )
            return;
        bool bAllPointsSame = true;
        vPoints.forRange(
            [ & ]( const point_t& rP1 ) {
                vPoints.forRange(
                    [ & ]( const point_t& rP2 ) {
                        for( size_t uiA = 0; uiA < d; uiA++ )
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
            // @todo should be own function
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "generating overlay: " << vTree.vLeaves.size( ) << " from: " << vOverlayBottomLeft 
                          << " to: " << vOverlayTopRight << std::endl;
            // recursion termination
            std::array<size_t, d> vBegins{ };
            std::array<size_t, d> vSizes{ };

            // for each dimension
            for( size_t uiI = 0; uiI < d; uiI++ )
            {
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\tgenerating axis for dimension=" << uiI << std::endl;
                // 2) fill a temp vec with the values
                auto vTmp = axis_vec_generator( );
                if constexpr( EXPLAIN_QUERY )
                {
                    std::cerr << "\tuiAxisCoordinatesStart/End=" << uiAxisCoordinatesStart << " "
                              << uiAxisCoordinatesEnd 
                              << " from/to pos " 
                              << vAxisCoordinates[uiI][vAxisFrom[uiI]] << " " 
                              << vAxisCoordinates[uiI][vAxisTo[uiI]-1] 
                              << " overlay start/end " 
                              << vOverlayBottomLeft[uiI] << " " 
                              << vOverlayTopRight[uiI]
                              << std::endl;
                              
                    vPoints.forRange(
                        [ & ]( const point_t& rP ) {
                            std::cerr << "\t\tpoint " << rP << std::endl;
                            return true;
                        },
                        uiFrom, uiTo );
                    vPointsEnd.forRange(
                        [ & ]( const point_t& rP ) {
                            std::cerr << "\t\tpoint end " << rP << std::endl;
                            return true;
                        },
                        uiFromEnd, uiToEnd );
                }

                pos_t vPos = vOverlayTopRight;
                for( size_t uiA = 0; uiA < d; uiA++ )
                    if(vPos[uiA] < std::numeric_limits<coordinate_t>::max())
                        ++vPos[uiA];
                vPos[ uiI ] = vOverlayBottomLeft[ uiI ];
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\toverlay init val pos=" << vPos << std::endl;
                data_t uiInitVal = countToZero<true>( xDatasetId, vPos );
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\toverlay axis init val=" 
                                << uiInitVal << std::endl;

                vPos = vOverlayBottomLeft;
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\toverlay bottom left corner pos=" << vPos << std::endl;
                uiInitVal -= countToZero<true>( xDatasetId, vPos );
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\toverlay axis init val - bottom left corner=" 
                                << uiInitVal << std::endl;

                vPoints.sortByLayerAndDim( uiI, uiFrom, uiTo );
                vPointsEnd.sortEndByLayerAndDim( uiI, uiFromEnd, uiToEnd );

                // @todo should not be lambda function
                auto fForPos = [&](coordinate_t uiCoord){
                    data_t uiVal = uiInitVal;

                    std::array<coordinate_t, d> vStart{};
                    std::array<coordinate_t, d> vEnd{};
                    pos_t vPos = vOverlayBottomLeft;
                    vPos[uiI] = uiCoord + 1;
                    if constexpr( EXPLAIN_QUERY )
                        std::cerr << "\toverlay axis pos=" << uiCoord << " precurser pos=" << vPos << std::endl;

                    uiVal += countToZero<true>( xDatasetId, vPos );

                    if constexpr( EXPLAIN_QUERY )
                        std::cerr << "\toverlay axis pos=" << uiCoord 
                                  << " init val - bottom left corner + precursers =" 
                                    << uiVal << std::endl;

                    // set point iterator to correct position
                    size_t uiPointsIt = uiFrom;
                    size_t uiPointsEndIt = uiFromEnd;
                    for(size_t uiA = 0; uiA < LAYERS; uiA++)
                    {
                        // @todo could be bin search
                        while(uiPointsIt < uiTo && vPoints.vData[ uiPointsIt ].uiLayer < uiA 
                              || (vPoints.vData[ uiPointsIt ].uiLayer == uiA 
                                  && vPoints.vData[ uiPointsIt ].vFrom[ uiI ] <= uiCoord))
                            ++uiPointsIt;
                        while(uiPointsEndIt < uiToEnd && vPointsEnd.vData[ uiPointsEndIt ].uiLayer < uiA 
                              || (vPointsEnd.vData[ uiPointsEndIt ].uiLayer == uiA 
                                  && vPointsEnd.vData[ uiPointsEndIt ].vTo[ uiI ] <= uiCoord))
                            ++uiPointsEndIt;

                        uiVal[uiA][0] += uiPointsIt - uiFrom;
                        uiVal[uiA][1] += uiPointsEndIt - uiFromEnd;

                        if constexpr( EXPLAIN_QUERY )
                            std::cerr << "\toverlay axis pos=" << uiCoord
                                    << " init val - bottom left corner + precursers + points ="
                                    << uiVal << std::endl;
                    }

                    if( ( vTmp.size() == 0 || vTmp.back().second < uiVal ) && uiVal > 0)
                        vTmp.emplace_back( uiCoord, uiVal );
                };

                fForPos( vOverlayBottomLeft[uiI] );
                for(size_t uiJ = vAxisFrom[uiI]; uiJ < vAxisTo[uiI]; uiJ++)
                    fForPos( vAxisCoordinates[uiI][uiJ] );

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
            // @todo should be own function

            // recursive call: split points in half

            // 1) find dimension that splits points most evenly in half
            size_t uiBestDim = 0;
            size_t uiBestCount = std::numeric_limits<size_t>::max( );
            for( size_t uiI = 0; uiI < d; uiI++ )
            {
                vPoints.sortByDim( uiI, uiFrom, uiTo );
                std::array<coordinate_t, b> uiAxisSplitPos;
                for( size_t uiA = 0; uiA < b; uiA++ )
                    uiAxisSplitPos[ uiA ] = vPoints.vData[ ( uiA * ( uiTo - uiFrom ) ) / b + uiFrom ].vPos[ uiI ];
                for( size_t uiA = 1; uiA < b; uiA++ )
                    while(uiAxisSplitPos[ uiA - 1 ] >= uiAxisSplitPos[ uiA ] && 
                          uiAxisSplitPos[ uiA ] < vOverlayTopRight[uiI])
                        ++uiAxisSplitPos[ uiA ];
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

            // 2.1) seperate points in that dimension (by sorting and looking for the split pos)
            vPoints.sortByDim( uiBestDim, uiFrom, uiTo );
            vPointsEnd.sortEndByDim( uiBestDim, uiFromEnd, uiToEnd );
            std::array<coordinate_t, b> uiAxisSplitPos;
            for( size_t uiA = 0; uiA < b; uiA++ )
                uiAxisSplitPos[ uiA ] = vPoints.vData[ ( uiA * ( uiTo - uiFrom ) ) / b + uiFrom ].vPos[ uiBestDim ];
            for( size_t uiA = 1; uiA < b; uiA++ )
                while(uiAxisSplitPos[ uiA - 1 ] >= uiAxisSplitPos[ uiA ] && 
                      uiAxisSplitPos[ uiA ] < vOverlayTopRight[uiBestDim])
                    ++uiAxisSplitPos[ uiA ];
            uiAxisSplitPos[ 0 ] = vOverlayBottomLeft[ uiBestDim ];
            std::array<size_t, b> uiCountInSplit{ };
            std::array<size_t, b> uiCountInSplitEnd{ };
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "splitting in dimension " << uiBestDim << " start "
                          << vPoints.vData[ uiFrom ].vPos[ uiBestDim ] << " end "
                          << vPoints.vData[ uiTo - 1 ].vPos[ uiBestDim ] << " with max count " << uiBestCount
                          << std::endl;
            // @todo could be a bin search
            vPoints.forRange(
                [ & ]( const point_t& rP ) {
                    for( size_t uiA = 1; uiA <= b; uiA++ )
                        if( rP.vPos[ uiBestDim ] >= uiAxisSplitPos[ b - uiA ] )
                        {
                            if constexpr( EXPLAIN_QUERY )
                                std::cerr << "point " << rP << std::endl;
                            uiCountInSplit[ b - uiA ] += 1;
                            break;
                        }
                    return true;
                },
                uiFrom, uiTo );
            vPointsEnd.forRange(
                [ & ]( const point_t& rP ) {
                    for( size_t uiA = 1; uiA <= b; uiA++ )
                        if( rP.vPos[ uiBestDim ] >= uiAxisSplitPos[ b - uiA ] )
                        {
                            if constexpr( EXPLAIN_QUERY )
                                std::cerr << "end point " << rP << std::endl;
                            uiCountInSplitEnd[ b - uiA ] += 1;
                            break;
                        }
                    return true;
                },
                uiFromEnd, uiToEnd );

            // 2.2) seperate axis coordinates into bins
            std::array<std::array<size_t, d>, b> vAxisFroms;
            std::array<std::array<size_t, d>, b> vAxisTos;
            for( size_t uiA = 0; uiA < b; uiA++ )
                for( size_t uiI = 0; uiI < d; uiI++ )
                {
                    vAxisFroms[uiA][uiI] = uiA == 0 ? vAxisFrom[uiI] : vAxisTos[uiA-1][uiI];
                    vAxisTos[uiA][uiI] = vAxisFroms[uiA][uiI];
                    // @todo could be a bin search
                    while(vAxisTos[uiA][uiI] < vAxisTo[uiI] && 
                          (uiA+1 < b || vAxisCoordinates[uiI][vAxisTos[uiA][uiI]] < uiAxisSplitPos[uiA+1]))
                        ++vAxisTos[uiA][uiI];
                }

            if constexpr( EXPLAIN_QUERY )
                for( size_t uiA = 0; uiA < b; uiA++ )
                    std::cerr << "\tat coordinate " << uiAxisSplitPos[ uiA ] << " containing " << uiCountInSplit[ uiA ]
                              << " points" << std::endl;

            // 3) insert new kd-tree nodes

            offset_t uiMyOffset = vTree.vTree.size( );
            vTree.vTree.emplace_back( );
            vTree.vTree[ uiMyOffset ].uiSplitDimension = uiBestDim;
            for( size_t uiA = 0; uiA < b; uiA++ )
            {
                std::get<0>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = std::numeric_limits<coordinate_t>::max( );
                std::get<1>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = std::numeric_limits<offset_t>::max( );
            }

            bool bRegistered = false;

            coordinate_t uiTopRightPre = vOverlayTopRight[ uiBestDim ];
            for( size_t uiA = 0; uiA < b; uiA++ )
            {

                vOverlayBottomLeft[ uiBestDim ] = uiAxisSplitPos[ uiA ];
                size_t uiB = uiA + 1;
                while(uiB < b && uiCountInSplit[ uiB ] == 0 )
                    uiB++;
                if( uiB < b )
                    vOverlayTopRight[ uiBestDim ] = uiAxisSplitPos[ uiB ];
                else
                    vOverlayTopRight[ uiBestDim ] = uiTopRightPre;

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
                    vOverlayTopRight, vAxisCoordinates, vAxisFroms[uiA], vAxisFroms[uiA], vPointsEnd, uiFromEnd,
                    uiFromEnd + uiCountInSplitEnd[uiA] );

                uiFrom += uiCountInSplit[ uiA ];
                uiFromEnd += uiCountInSplitEnd[uiA]
            }
            if(!bRegistered)
                fRegisterMeInTree( uiMyOffset, false );
            assert( uiFrom == uiTo );
            assert( uiFromEnd == uiToEnd );
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
    val_t countToZero( const overlay_meta_t& rO, pos_t vPos, std::array<bool, 2> vInclusive ) const
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
                for(size_t uiI = 0; uiI < d; uiI++)
                    bToTopRight = bToTopRight && ( rP.vPos[ uiI ] > vPos[ uiI ] );
                if( bToTopRight )
                    uiVal += 1;
                return true;
            },
            rO.uiPointsBegin, rO.uiPointsEnd );
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\toverlay number of points: " << uiVal
                      << " of total: " << rO.uiPointsEnd - rO.uiPointsBegin << std::endl;
        uiRet += uiVal;

        uiVal = vEntries.getLast( rO.vEntryBegins[ 0 ], rO.vSizes[ 0 ] );
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\toverlay top right corner: value= -" << uiVal << std::endl;
        assert( uiRet >= uiVal ); // can never have a negative return number
        uiRet -= uiVal;
        return uiRet;
    }

    template<bool SILENT>
    val_t countToZero( class_key_t xDatasetId, pos_t vPos, std::array<bool, 2> vInclusive ) const
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
                       size_t uiDCurr, size_t uiNumStart, std::array<bool, 2> vInclusive ) const
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
        // @todo start with empty and build on the go
        std::array<cord_vec_t, d> vAxisCoordinates{ };
        std::array<size_t, d> vFrom{ };
        std::array<size_t, d> vTo{ };
        point_vec_t vPointsEnd = point_vec_generator();
        vPoints.forRange(
            [ & ]( const point_t& rP ) {
                vPointsEnd.push_back(rP);
                return true;
            },
            uiFrom, uiTo );

        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            vOverlayTopRight[ uiI ] = std::numeric_limits<coordinate_t>::max( );

            cord_vec_t vTmpFrom = cord_vec_generator( );
            vPoints.sortByDim(uiI, uiFrom, uiTo);
            vPoints.forRange(
                [ & ]( const point_t& rP2 ) {
                    if( vTmpFrom.size() == 0 || vTmpFrom.back() < rP2.vFrom[ uiI ] )
                        vTmpFrom.push_back( rP2.vFrom[ uiI ] );
                    return true;
                },
                uiFrom, uiTo );
            cord_vec_t vTmpTo = cord_vec_generator( );
            vPoints.sortEndByDim(uiI, uiFrom, uiTo);
            vPoints.forRange(
                [ & ]( const point_t& rP2 ) {
                    if( vTmpTo.size() == 0 || vTmpTo.back() < rP2.vTo[ uiI ] )
                        vTmpTo.push_back( rP2.vTo[ uiI ] );
                    return true;
                },
                uiFrom, uiTo );
            vAxisCoordinates[ uiI ] = cord_vec_generator( );
            size_t uiF = 0; uiT = 0;
            while(uiF < vTmpFrom.size() || uiT < vTmpFrom.size() )
            {
                if(uiF < vTmpFrom.size() && (uiT == vTmpFrom.size() || vTmpFrom[uiF] < vTmpTo[uiT]))
                    vAxisCoordinates[ uiI ].push_back(vTmpFrom[uiF++]);
                else if(uiT < vTmpFrom.size() && (uiF == vTmpFrom.size() || vTmpTo[uiT] < vTmpFrom[uiF]))
                    vAxisCoordinates[ uiI ].push_back(vTmpTo[uiT++]);
                else if(uiT < vTmpFrom.size() && uiF < vTmpFrom.size() && vTmpTo[uiT] == vTmpFrom[uiF]))
                {
                    vAxisCoordinates[ uiI ].push_back(vTmpFrom[uiF++]);
                    uiT++
                }
                else
                    assert(false); // should never reach this point
            }
            vTo[uiI] = vCoordinates[uiI].size();
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
            vAxisCoordinates, //
            vFrom, //
            vTo, //
            vPointsEnd, //
            0, //
            vPointsEnd.size() );
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

    void iterate( std::function<void( pos_t, pos_t, std::string )> fDo, class_key_t xDatasetId, pos_t vFrom, 
                  pos_t vTo ) const
    {
        iterateOverlaysIn(
            [ & ]( const overlay_meta_t& xO ) {
                vPoints.forRange(
                    [ & ]( const point_t& xPoint ) {
                        for( size_t uiI = 0; uiI < d; uiI++ )
                        {
                            if( xPoint.vFrom[ uiI ] >= vTo[ uiI ] )
                                return true;
                            if( xPoint.vTo[ uiI ] < vFrom[ uiI ] )
                                return true;
                        }
                        fDo( xPoint.vFrom, xPoint.vTo, vDesc.get( xPoint.uiDescOffset ) );
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