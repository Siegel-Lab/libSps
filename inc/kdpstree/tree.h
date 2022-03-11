#pragma once

#include "kdpstree/overlay_kd_tree.h"
#include "kdpstree/overlay_meta.h"
#include "kdpstree/point.h"
#include "kdpstree/points.h"
#include "kdpstree/type_defs.h"
#include "kdpstree/util.h"
#include <iostream>
#include <limits>
#include <sstream>


#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif

namespace kdpstree
{

template<typename T, size_t N>
std::array<T, N>& operator+=(std::array<T, N>& vA, const std::array<T, N>& vB)
{
    for(size_t uiI = 0; uiI < N; uiI++)
        vA[uiI] += vB[uiI];
    return vA;
}

template<typename T, size_t N>
std::array<T, N>& operator-=(std::array<T, N>& vA, const std::array<T, N>& vB)
{
    for(size_t uiI = 0; uiI < N; uiI++)
        vA[uiI] -= vB[uiI];
    return vA;
}


template<typename T, size_t N>
bool oneSmaller(std::array<T, N>& vA, std::array<T, N>& vB)
{
    for(size_t uiI = 0; uiI < N; uiI++)
        if(vA[uiI] < vB[uiI])
            return true;
    return false;
}

template<typename T, size_t N>
bool oneLarger(std::array<T, N>& vA, T uiAs)
{
    for(size_t uiI = 0; uiI < N; uiI++)
        if(vA[uiI] > uiAs)
            return true;
    return false;
}

template<typename T, size_t N>
std::array<T, N>& operator*=(std::array<T, N>& vA, T uiB)
{
    for(size_t uiI = 0; uiI < N; uiI++)
        vA[uiI] *= uiB;
    return vA;
}

template <typename type_defs> class Tree;

template <typename type_defs>
std::ostream& operator<<(std::ostream& os, const Tree<type_defs>& rTree);

template <typename type_defs> class Tree
{
    EXTRACT_TYPE_DEFS; // macro call

    using overlay_meta_t = OverlayMeta<type_defs>;
    using overlay_kd_tree_t = OverlayKdTree<type_defs>;
    using overlay_entries_t = OverlayEntries<type_defs>;
    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;
    using desc_t = Desc<type_defs>;

    
    EXTRACT_TMP_VEC_GENERATOR(cord, coordinate_t); // macro call

    using axis_t = std::pair<coordinate_t, data_t>;
    EXTRACT_TMP_VEC_GENERATOR(axis, axis_t); // macro call

    friend std::ostream& operator<< <>(std::ostream& os, const Tree& rTree);

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
                                  std::array<size_t, d> vAxisTo
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
                auto vTmp = axis_vec_generator.vec( );
                if constexpr( EXPLAIN_QUERY )
                {
                    std::cerr << "\tuiAxisCoordinatesStart/End=" << vAxisFrom[uiI] << " "
                              << vAxisTo[uiI] 
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
                }

                pos_t vPos = vOverlayTopRight;
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
                                << uiInitVal << std::endl << "\tuiFrom: " << uiFrom << " uiTo: " << uiTo << std::endl;

                vPoints.sortByLayerAndDim( uiI, uiFrom, uiTo );

                if constexpr( EXPLAIN_QUERY )
                    vPoints.forRange(
                        [ & ]( const point_t& rP ) {
                            std::cerr << "\t\tsorted point " << rP << std::endl;
                            return true;
                        },
                        uiFrom, uiTo );

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
                    size_t uiPointsStartIt = uiFrom;
                    for(size_t uiA = 0; uiA < LAYERS; uiA++)
                    {
                        // @todo could be bin search
                        while(uiPointsStartIt < uiTo && vPoints.vData[ uiPointsStartIt ].uiLayer < uiA)
                            ++uiPointsStartIt;
                        size_t uiPointsEndIt = uiPointsStartIt;
                        while(uiPointsEndIt < uiTo && vPoints.vData[ uiPointsEndIt ].uiLayer == uiA 
                                  && vPoints.vData[ uiPointsEndIt ].vPos[ uiI ] <= uiCoord)
                            ++uiPointsEndIt;

                        uiVal[uiA] += uiPointsEndIt - uiPointsStartIt;

                    }
                    if constexpr( EXPLAIN_QUERY )
                        std::cerr << "\toverlay axis pos=" << uiCoord
                                << " init val - bottom left corner + precursers + points ="
                                << uiVal << std::endl;

                    if( ( vTmp.size() == 0 || oneSmaller(vTmp.back().second, uiVal) ) && oneLarger(uiVal, (val_t)0))
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
            vTree.vLeaves.push_back( overlay_meta_t( vBegins, vSizes, uiFrom, uiTo ) );

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
            std::array<coordinate_t, b> uiAxisSplitPos;
            for( size_t uiA = 0; uiA < b; uiA++ )
                uiAxisSplitPos[ uiA ] = vPoints.vData[ ( uiA * ( uiTo - uiFrom ) ) / b + uiFrom ].vPos[ uiBestDim ];
            for( size_t uiA = 1; uiA < b; uiA++ )
                while(uiAxisSplitPos[ uiA - 1 ] >= uiAxisSplitPos[ uiA ] && 
                      uiAxisSplitPos[ uiA ] < vOverlayTopRight[uiBestDim])
                    ++uiAxisSplitPos[ uiA ];
            uiAxisSplitPos[ 0 ] = vOverlayBottomLeft[ uiBestDim ];
            std::array<size_t, b> uiCountInSplit{ };
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

            // 2.2) seperate axis coordinates into bins
            std::array<std::array<size_t, d>, b> vAxisFroms;
            std::array<std::array<size_t, d>, b> vAxisTos;
            size_t uiA = 0;
            while( uiA < b && uiCountInSplit[uiA] == 0 )
                ++uiA;
            size_t uiLastA = 0;
            while( uiA < b )
            {
                size_t uiNextA = uiA + 1;
                while(uiNextA < b && uiCountInSplit[uiNextA] == 0)
                    ++uiNextA;

                for( size_t uiI = 0; uiI < d; uiI++ )
                {
                    vAxisFroms[uiA][uiI] = vAxisFrom[uiI];
                    vAxisTos[uiA][uiI] = vAxisTo[uiI];
                }
                vAxisFroms[uiA][uiBestDim] = uiA == 0 ? vAxisFrom[uiBestDim] : vAxisTos[uiLastA][uiBestDim];
                vAxisTos[uiA][uiBestDim] = vAxisFroms[uiA][uiBestDim];
                // @todo could be a bin search
                while(vAxisTos[uiA][uiBestDim] < vAxisTo[uiBestDim] && 
                        (uiNextA == b || vAxisCoordinates[uiBestDim][vAxisTos[uiA][uiBestDim]] < uiAxisSplitPos[uiNextA]))
                    ++vAxisTos[uiA][uiBestDim];

                if constexpr( EXPLAIN_QUERY )
                {
                    std::cerr << "\taxis from " << vAxisFroms[uiA] << " to " << vAxisTos[uiA] << std::endl;
                    std::cerr << "\tcoords: " << vAxisCoordinates[uiBestDim][vAxisFroms[uiA][uiBestDim]] << " " 
                              << vAxisCoordinates[uiBestDim][vAxisTos[uiA][uiBestDim]-1] << std::endl;
                    std::cerr << "\tat coordinate " << uiAxisSplitPos[ uiA ] << " containing " << uiCountInSplit[ uiA ]
                              << " points" << std::endl;
                }
                uiLastA = uiA;
                uiA = uiNextA;
            }

            // 3) insert new kd-tree nodes

            offset_t uiMyOffset = vTree.vTree.size( );
            vTree.vTree.push_back( typename OverlayKdTree<type_defs>::OverlayKdTreeNode(uiBestDim) );
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
                    vOverlayTopRight, vAxisCoordinates, vAxisFroms[uiA], vAxisTos[uiA] );

                uiFrom += uiCountInSplit[ uiA ];
            }
            if(!bRegistered)
                fRegisterMeInTree( uiMyOffset, false );
            assert( uiFrom == uiTo );
            
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "split in dimension " << uiBestDim << " done. (" << uiMyOffset << ": " 
                          << (size_t)vTree.vTree[uiMyOffset].uiSplitDimension << ")" << std::endl;
        }
    }

    template<bool SILENT>
    std::pair<pos_t, const overlay_meta_t*> getOverlay( const class_key_t& xDatasetId, const pos_t& vPos ) const
    {
        auto cIter = vTree.vRoots.begin( );
        while(cIter != vTree.vRoots.end())
        {
            if( std::get<0>( *cIter ) == xDatasetId )
            {
                pos_t vBottomLeft{ };
                size_t uiCurr = std::get<1>( *cIter );
                bool bIsLeaf = std::get<2>( *cIter );
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
            ++cIter;
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
    data_t countToZero( const overlay_meta_t& rO, pos_t vPos ) const
    {
        data_t uiRet{};
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            if constexpr( EXPLAIN_QUERY && !SILENT )
                std::cerr << "\t\toverlay value: dimension=" << uiI << " position=" << vPos[ uiI ] << std::endl;
            data_t uiVal = vEntries.get( vPos[ uiI ], rO.vEntryBegins[ uiI ], rO.vSizes[ uiI ] );
            if constexpr( EXPLAIN_QUERY && !SILENT )
                std::cerr << "\t\toverlay value=" << uiVal << std::endl;
            uiRet += uiVal;
        }

        data_t uiVal{};
        vPoints.forRange(
            [ & ]( const point_t& rP ) {
                bool bToTopRight = true;
                for(size_t uiI = 0; uiI < d; uiI++)
                    bToTopRight = bToTopRight && ( rP.vPos[ uiI ] > vPos[ uiI ] );
                if( bToTopRight )
                    uiVal[rP.uiLayer] += 1;
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
    data_t countToZero( class_key_t xDatasetId, pos_t vPos ) const
    {
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\tcountToZero" << std::endl;
        if( vTree.vRoots.size( ) == 0 )
            return data_t{}; // return zero
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            if( vPos[ uiI ] == 0 )
            {
                if constexpr( EXPLAIN_QUERY && !SILENT )
                    std::cerr << "\t\tdimension: " << uiI << " is zero. result must be zero as well." << std::endl;
                return data_t{}; // return zero
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
            return data_t{}; // return zero
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
    data_t countHelper( const class_key_t& xDatasetId, const pos_t& vFrom, const pos_t& vTo, pos_t& vCurr,
                       size_t uiDCurr, size_t uiNumStart ) const
    {
        if( uiDCurr == d )
        {
            if constexpr( EXPLAIN_QUERY && !SILENT )
            {
                std::cerr << "\tcounting from (";
                for( size_t uiI = 0; uiI < d; uiI++ )
                    std::cerr << "0, ";
                std::cerr << ") to " << vCurr;
                if( uiNumStart % 2 == 0 )
                    std::cerr << " positive";
                else
                    std::cerr << " negative";
                std::cerr << std::endl;
            }
            auto uiVal = countToZero<SILENT>( xDatasetId, vCurr );
            if constexpr( EXPLAIN_QUERY && !SILENT )
                std::cerr << "\tresult = " << uiVal << std::endl;
            uiVal *= (val_t)( uiNumStart % 2 == 0 ? 1 : -1 );
            return uiVal;
        }
        else
        {
            vCurr[ uiDCurr ] = vFrom[ uiDCurr ];
            data_t uiRet = countHelper<SILENT>( xDatasetId, vFrom, vTo, vCurr, uiDCurr + 1, uiNumStart + 1 );
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

    size_t addPoint( pos_t vPos, layers_t uiLayer, std::string sDesc )
    {
        vPoints.add( vPos, vDesc.add( sDesc ), uiLayer );
    }

    void generateForPoints( class_key_t xDatasetId, size_t uiMaxPointsPerOverlay, size_t uiFrom, size_t uiTo )
    {
        std::array<coordinate_t, d> vOverlayBottomLeft{ };
        std::array<coordinate_t, d> vOverlayTopRight{ };
        std::array<cord_vec_t, d> vAxisCoordinates{ };
        std::array<size_t, d> vFrom{ };
        std::array<size_t, d> vTo{ };

        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            vOverlayTopRight[ uiI ] = std::numeric_limits<coordinate_t>::max( );

            vPoints.sortByDim(uiI, uiFrom, uiTo);
            vAxisCoordinates[ uiI ] = cord_vec_generator.vec( );
            vPoints.forRange(
                [ & ]( const point_t& rP2 ) {
                    if( vAxisCoordinates[ uiI ].size() == 0 || vAxisCoordinates[ uiI ].back() < rP2.vPos[ uiI ] )
                        vAxisCoordinates[ uiI ].push_back( rP2.vPos[ uiI ] );
                    return true;
                },
                uiFrom, uiTo );
            vTo[uiI] = vAxisCoordinates[uiI].size();
        }
        generateForPointsHelper(
            xDatasetId, //
            [ & ]( size_t uiOffsetNode, bool bIsLeaf ) {
                vTree.vRoots.push_back( std::make_tuple( xDatasetId, uiOffsetNode, bIsLeaf ) );
            },
            uiMaxPointsPerOverlay, //
            uiFrom, //
            uiTo, //
            vOverlayBottomLeft, //
            vOverlayTopRight, //
            vAxisCoordinates, //
            vFrom, //
            vTo );
    }

    template<bool SILENT>
    data_t count( class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        pos_t vCurr{ };
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "counting query for datsetId=" << xDatasetId << " bottomLeftCorner=" << vFrom
                      << " topRightCorner=" << vTo << std::endl;
        return countHelper<SILENT>( xDatasetId, vFrom, vTo, vCurr, 0, 0 );
    }

    void iterate( std::function<void( pos_t, layers_t, std::string )> fDo, class_key_t xDatasetId, pos_t vFrom, 
                  pos_t vTo ) const
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
                        fDo( xPoint.vPos, xPoint.uiLayer, vDesc.get( xPoint.uiDescOffset ) );
                        return true;
                    },
                    xO.uiPointsBegin, xO.uiPointsEnd );
            },
            xDatasetId, vFrom, vTo );
    }

    std::array<std::vector<std::pair<pos_t, std::string>>, LAYERS> get( class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        std::array<std::vector<std::pair<pos_t, std::string>>, LAYERS> vRet;
        iterate( [ & ]( pos_t vPos, layers_t uiLayer, std::string sDesc ) { vRet[uiLayer].emplace_back( vPos, sDesc ); }, xDatasetId, vFrom, vTo );
        return vRet;
    }

    std::string print()
    {
        std::stringstream ss;
        ss << *this << std::endl;
        return ss.str();
    }
};


template <typename type_defs>
std::ostream& operator<<(std::ostream& os, const Tree<type_defs>& rTree)
{
    os << "Tree:" << std::endl << rTree.vTree 
       << "Entries:" << std::endl << rTree.vEntries 
       << "Points:" << std::endl << rTree.vPoints
       <<"Desc:" << std::endl << rTree.vDesc << std::endl;
    return os;
}

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