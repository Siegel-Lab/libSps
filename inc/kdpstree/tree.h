#pragma once

#include "kdpstree/util.h"

#include "kdpstree/overlay_kd_tree.h"
#include "kdpstree/overlay_meta.h"
#include "kdpstree/point.h"
#include "kdpstree/points.h"
#include "kdpstree/type_defs.h"
#include <iostream>
#include <limits>
#include <sstream>


#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
//#include "kdpstree/util/pybind11.h"
#endif

namespace kdpstree
{

template <typename T, size_t N> std::array<T, N>& operator+=( std::array<T, N>& vA, const std::array<T, N>& vB )
{
    for( size_t uiI = 0; uiI < N; uiI++ )
        vA[ uiI ] += vB[ uiI ];
    return vA;
}

template <typename T, size_t N> std::array<T, N>& operator-=( std::array<T, N>& vA, const std::array<T, N>& vB )
{
    for( size_t uiI = 0; uiI < N; uiI++ )
        vA[ uiI ] -= vB[ uiI ];
    return vA;
}


template <typename T, size_t N> bool oneSmaller( std::array<T, N>& vA, std::array<T, N>& vB )
{
    for( size_t uiI = 0; uiI < N; uiI++ )
        if( vA[ uiI ] < vB[ uiI ] )
            return true;
    return false;
}

template <typename T, size_t N> bool oneLarger( std::array<T, N>& vA, T uiAs )
{
    for( size_t uiI = 0; uiI < N; uiI++ )
        if( vA[ uiI ] > uiAs )
            return true;
    return false;
}

template <typename T, size_t N> std::array<T, N>& operator*=( std::array<T, N>& vA, T uiB )
{
    for( size_t uiI = 0; uiI < N; uiI++ )
        vA[ uiI ] *= uiB;
    return vA;
}

template <typename type_defs> class Tree;

} // namespace kdpstree

namespace std
{

template <typename type_defs> ostream& operator<<( ostream& os, const typename kdpstree::Tree<type_defs>& rTree );

}

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


    EXTRACT_TMP_VEC_GENERATOR( cord, coordinate_t ); // macro call

    using axis_t = std::pair<coordinate_t, data_t>;
    EXTRACT_TMP_VEC_GENERATOR( axis, axis_t ); // macro call

    friend std::ostream& std::operator<<<>( std::ostream& os, const Tree& rTree );

    overlay_kd_tree_t vTree;
    overlay_entries_t vEntries;
    points_t vPoints;
    desc_t vDesc;

    /**
     * @brief  see https://algorithmica.org/en/eytzinger
     * @todo this should be a static B-tree, where
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
                                  const std::array<cord_vec_t, d>& vAxisCoordinates, std::array<size_t, d> vAxisFrom,
                                  std::array<size_t, d> vAxisTo, std::optional<progress_stream_t>& xProg, 
                                  Profiler& rProfiler )
    {
        if( uiFrom == uiTo )
            return;
        rProfiler.step("overlay - check points same");
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
                rProfiler.step("overlay - init");
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\tgenerating axis for dimension=" << uiI << std::endl;
                // 2) fill a temp vec with the values
                auto vTmp = axis_vec_generator.vec( );
                if constexpr( EXPLAIN_QUERY )
                {
                    std::cerr << "\tuiAxisCoordinatesStart/End=" << vAxisFrom[ uiI ] << " " << vAxisTo[ uiI ]
                              << " from/to pos " << vAxisCoordinates[ uiI ][ vAxisFrom[ uiI ] ] << " "
                              << vAxisCoordinates[ uiI ][ vAxisTo[ uiI ] - 1 ] << " overlay start/end "
                              << vOverlayBottomLeft[ uiI ] << " " << vOverlayTopRight[ uiI ] << std::endl;

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
                    std::cerr << "\toverlay axis init val=" << uiInitVal << std::endl;

                vPos = vOverlayBottomLeft;
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\toverlay bottom left corner pos=" << vPos << std::endl;
                uiInitVal -= countToZero<true>( xDatasetId, vPos );
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\toverlay axis init val - bottom left corner=" << uiInitVal << std::endl
                              << "\tuiFrom: " << uiFrom << " uiTo: " << uiTo << std::endl;

                vPoints.sortByLayerAndDim( uiI, uiFrom, uiTo );

                if constexpr( EXPLAIN_QUERY )
                    vPoints.forRange(
                        [ & ]( const point_t& rP ) {
                            std::cerr << "\t\tsorted point " << rP << std::endl;
                            return true;
                        },
                        uiFrom, uiTo );

                std::array<offset_t, LAYERS> vPointsStartIt;
                std::array<offset_t, LAYERS> vPointsCurrIt;
                for( size_t uiA = 0; uiA < LAYERS; uiA++ )
                {
                    vPointsStartIt[ uiA ] =
                        vPoints.lowerBound( uiA > 0 ? vPointsStartIt[ uiA - 1 ] : uiFrom, uiTo,
                                            [ & ]( const point_t& rP ) { return rP.uiLayer < uiA; } );
                    vPointsCurrIt[ uiA ] = vPointsStartIt[ uiA ];
                }

                const overlay_meta_t* pOverlayCache = nullptr;
                std::vector<std::pair<coordinate_t, const data_t>> vCache;
                size_t uiCacheIt = 0;
                coordinate_t uiOverlayCacheTop = 0;

                const auto& rvPData = vPoints.vData;
                coordinate_t uiCoord = vOverlayBottomLeft[ uiI ];
                while(uiCoord < vOverlayTopRight[ uiI ])
                {
                    rProfiler.step("overlay - tmp vec - cache");
                    while(uiCacheIt < vCache.size() && vCache[uiCacheIt].first < uiCoord)
                        ++uiCacheIt;

                    // setup cache if necessary
                    if(vCache.size() == 0 || uiOverlayCacheTop <= uiCoord)
                    {
                        vCache.clear();
                        pos_t vPos = vOverlayBottomLeft;
                        vPos[ uiI ] = uiCoord + 1;
                        bool bOk = true;
                        for( size_t uiI = 0; uiI < d && bOk; uiI++ )
                        {
                            if(vPos[ uiI ] == 0)
                                bOk = false;
                            else
                                --vPos[ uiI ];
                        }
                        if (bOk)
                        {
                            if constexpr( EXPLAIN_QUERY )
                                std::cerr << "\toverlay axis pos=" << uiCoord << " precurser pos= " << vPos 
                                          << std::endl;
                            auto rO = getOverlay<true>( xDatasetId, vPos );
                            pOverlayCache = std::get<2>(rO);
                            uiOverlayCacheTop = std::get<1>(rO)[ uiI ];
                            vEntries.forRange(
                                [&](coordinate_t vPos, const data_t& vData){
                                    vCache.emplace_back(vPos, vData);
                                },
                                std::get<2>(rO)->vEntryBegins[ uiI ], 
                                std::get<2>(rO)->vSizes[ uiI ],
                                uiCoord,
                                vOverlayTopRight[ uiI ]
                            );
                            uiCacheIt = 0;
                            assert(vCache[uiCacheIt].first >= uiCoord);
                        }
                        else
                        {
                            // if there are no more overlays and no more points we can cancel the loop as well
                            bool bOK = true;
                            for( size_t uiA = 0; uiA < LAYERS && bOK; uiA++ )
                                if(vPointsCurrIt[ uiA ] < (uiA + 1 < LAYERS ? vPointsStartIt[ uiA + 1 ] : uiTo))
                                    bOK = false;
                            if(bOK)
                                break;
                            else
                            {
                                coordinate_t uiNewCoord = std::numeric_limits<coordinate_t>::max();
                                for( size_t uiA = 0; uiA < LAYERS; uiA++ )
                                    if(vPointsCurrIt[ uiA ] < (uiA + 1 < LAYERS ? vPointsStartIt[ uiA + 1 ] : uiTo) && 
                                        rvPData[ vPointsCurrIt[ uiA ] ].vPos[ uiI ] < uiNewCoord)
                                        uiNewCoord = rvPData[ vPointsCurrIt[ uiA ] ].vPos[ uiI ];
                                if(uiNewCoord < std::numeric_limits<coordinate_t>::max())
                                uiCoord = std::max(uiCoord, uiNewCoord);
                            }
                        }
                    }

                    rProfiler.step("overlay - tmp vec - points");
                    for( size_t uiA = 0; uiA < LAYERS; uiA++ )
                    {
                        offset_t uiPointsEndIt = vPointsCurrIt[ uiA ];
                        while(uiPointsEndIt < (uiA + 1 < LAYERS ? vPointsStartIt[ uiA + 1 ] : uiTo) &&
                              rvPData[uiPointsEndIt].uiLayer <= uiA && 
                              rvPData[uiPointsEndIt].vPos[ uiI ] <= uiCoord)
                              ++uiPointsEndIt;
                        assert( vPoints.vData[ vPointsCurrIt[ uiA ] ].uiLayer == uiA ||
                                vPointsCurrIt[ uiA ] == uiPointsEndIt );

                        if constexpr( EXPLAIN_QUERY )
                            std::cerr << "\t\tpos=" << uiCoord << " layer=" << uiA << " from=" << vPointsStartIt[ uiA ]
                                      << " to=" << uiPointsEndIt << std::endl;

                        vPointsCurrIt[ uiA ] = uiPointsEndIt;
                    }

                    coordinate_t uiNewCoord = std::numeric_limits<coordinate_t>::max();
                    if(uiCacheIt < vCache.size())
                        uiNewCoord = vCache[uiCacheIt].first;
                    for( size_t uiA = 0; uiA < LAYERS; uiA++ )
                        if(vPointsCurrIt[ uiA ] > vPointsStartIt[ uiA ] && 
                           rvPData[ vPointsCurrIt[ uiA ]-1 ].vPos[ uiI ] < uiNewCoord)
                            uiNewCoord = rvPData[ vPointsCurrIt[ uiA ]-1 ].vPos[ uiI ];
                    uiCoord = std::max(uiCoord, uiNewCoord);

                    // compute value
                    rProfiler.step("overlay - tmp vec - comp");
                    data_t uiVal = uiInitVal;
                    if(uiCacheIt < vCache.size() && vCache[uiCacheIt].first == uiCoord)
                        uiVal += vCache[uiCacheIt].second;

                    if constexpr( EXPLAIN_QUERY )
                        std::cerr << "\toverlay axis pos=" << uiCoord
                                  << " init val - bottom left corner + precursers = " << uiVal << std::endl;

#ifndef NDEBUG
                    {
                        pos_t vPos = vOverlayBottomLeft;
                        vPos[ uiI ] = uiCoord + 1;
                        bool bOk = true;
                        for( size_t uiI = 0; uiI < d && bOk; uiI++ )
                        {
                            if(vPos[ uiI ] == 0)
                                bOk = false;
                            else
                                --vPos[ uiI ];
                        }
                        if(bOk)
                        {
                            data_t uiX = countToZero<true>( *pOverlayCache, vPos );
                            uiX += uiInitVal;
                            assert(uiX == uiVal);
                        }
                    }
                    {
                        pos_t vPos = vOverlayBottomLeft;
                        vPos[ uiI ] = uiCoord + 1;
                        data_t uiX = countToZero<true>( xDatasetId, vPos );
                        uiX += uiInitVal;
                        assert(uiX == uiVal);
                    }
#endif

                    for( size_t uiA = 0; uiA < LAYERS; uiA++ )
                        uiVal[ uiA ] += vPointsCurrIt[ uiA ] - vPointsStartIt[ uiA ];
                        
                    if constexpr( EXPLAIN_QUERY )
                        std::cerr << "\toverlay axis pos=" << uiCoord
                                  << " init val - bottom left corner + precursers + points = " << uiVal << std::endl;

                    if( ( vTmp.size( ) == 0 || oneSmaller( vTmp.back( ).second, uiVal ) ) &&
                                oneLarger( uiVal, (val_t)0 ) )
                        vTmp.emplace_back( uiCoord, uiVal );

                    if(uiCoord < std::numeric_limits<coordinate_t>::max())
                        ++uiCoord;
                }

                for( size_t uiA = 0; uiA < LAYERS - 1; uiA++ )
                    assert( vPointsCurrIt[ uiA ] == vPointsStartIt[ uiA + 1 ] );
                assert( vTmp.size( ) > 0 );

                if constexpr( EXPLAIN_QUERY )
                    for( size_t uiX = 0; uiX < vTmp.size( ); uiX++ )
                        std::cerr << "\tfinal: overlay axis pos=" << vTmp[ uiX ].first << " val=" << vTmp[ uiX ].second
                                  << std::endl;

                // 3) save pointers to the constructed axes
                rProfiler.step("overlay - eytzinger");
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

            if( xProg && xProg->printAgain( ) )
            {
                rProfiler.print("overlay");
                xProg << CLRLN << "Generated overlays for " << uiTo << " out of " << vPoints.size( )
                      << " points. That's " << ( 100.0 * uiTo ) / vPoints.size( ) << "%...";
            }
        }
        else
        {
            rProfiler.step("overlay - split");
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
                    while( uiAxisSplitPos[ uiA - 1 ] >= uiAxisSplitPos[ uiA ] &&
                           uiAxisSplitPos[ uiA ] < vOverlayTopRight[ uiI ] )
                        ++uiAxisSplitPos[ uiA ];
                std::array<size_t, b> uiCountInSplit{ };

                offset_t uiCurrPos = uiFrom;
                for( size_t uiA = 0; uiA + 1 < b; uiA++ )
                {
                    offset_t uiX = vPoints.lowerBound( uiCurrPos, uiTo, [ & ]( const point_t& rP ) {
                        return rP.vPos[ uiI ] < uiAxisSplitPos[ uiA + 1 ];
                    } );
                    uiCountInSplit[ uiA ] = uiX - uiCurrPos;
                    uiCurrPos = uiX;
                }
                uiCountInSplit[ b - 1 ] = uiTo - uiCurrPos;

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
                while( uiAxisSplitPos[ uiA - 1 ] >= uiAxisSplitPos[ uiA ] &&
                       uiAxisSplitPos[ uiA ] < vOverlayTopRight[ uiBestDim ] )
                    ++uiAxisSplitPos[ uiA ];
            uiAxisSplitPos[ 0 ] = vOverlayBottomLeft[ uiBestDim ];
            std::array<size_t, b> uiCountInSplit{ };
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "splitting in dimension " << uiBestDim << " start "
                          << vPoints.vData[ uiFrom ].vPos[ uiBestDim ] << " end "
                          << vPoints.vData[ uiTo - 1 ].vPos[ uiBestDim ] << " with max count " << uiBestCount
                          << std::endl;
            offset_t uiCurrPos = uiFrom;
            for( size_t uiA = 0; uiA + 1 < b; uiA++ )
            {
                offset_t uiX = vPoints.lowerBound( uiCurrPos, uiTo, [ & ]( const point_t& rP ) {
                    return rP.vPos[ uiBestDim ] < uiAxisSplitPos[ uiA + 1 ];
                } );
                uiCountInSplit[ uiA ] = uiX - uiCurrPos;
                uiCurrPos = uiX;
            }
            uiCountInSplit[ b - 1 ] = uiTo - uiCurrPos;
            if constexpr( EXPLAIN_QUERY )
                vPoints.forRange(
                    [ & ]( const point_t& rP ) {
                        std::cerr << "point " << rP << std::endl;
                        return true;
                    },
                    uiFrom, uiTo );

            // 2.2 move all empty bins to end of array
            size_t uiNextNonEmpty = 0;
            for( size_t uiA = 0; uiA < b; uiA++ )
            {
                while( uiNextNonEmpty < b && uiCountInSplit[ uiNextNonEmpty ] == 0 )
                    ++uiNextNonEmpty;

                if( uiNextNonEmpty == b )
                {
                    uiAxisSplitPos[ uiA ] = vOverlayTopRight[ uiBestDim ];
                    uiCountInSplit[ uiA ] = 0;
                }
                else
                {
                    uiAxisSplitPos[ uiA ] = uiAxisSplitPos[ uiNextNonEmpty ];
                    uiCountInSplit[ uiA ] = uiCountInSplit[ uiNextNonEmpty ];
                }
                if( uiNextNonEmpty < b )
                    ++uiNextNonEmpty;
            }


            // 2.3) seperate axis coordinates into bins
            std::array<std::array<size_t, d>, b> vAxisFroms;
            std::array<std::array<size_t, d>, b> vAxisTos;
            size_t uiA = 0;
            while( uiA < b && uiCountInSplit[ uiA ] == 0 )
                ++uiA;
            size_t uiLastA = 0;
            while( uiA < b )
            {
                size_t uiNextA = uiA + 1;
                while( uiNextA < b && uiCountInSplit[ uiNextA ] == 0 )
                    ++uiNextA;

                for( size_t uiI = 0; uiI < d; uiI++ )
                {
                    vAxisFroms[ uiA ][ uiI ] = vAxisFrom[ uiI ];
                    vAxisTos[ uiA ][ uiI ] = vAxisTo[ uiI ];
                }
                vAxisFroms[ uiA ][ uiBestDim ] = uiA == 0 ? vAxisFrom[ uiBestDim ] : vAxisTos[ uiLastA ][ uiBestDim ];

                if( uiNextA < b )
                {
                    auto xEnd = vAxisTo[ uiBestDim ] < vAxisCoordinates[ uiBestDim ].size( )
                                    ? vAxisCoordinates[ uiBestDim ].begin( ) + vAxisTo[ uiBestDim ]
                                    : vAxisCoordinates[ uiBestDim ].end( );
                    auto xRet =
                        std::lower_bound( vAxisCoordinates[ uiBestDim ].begin( ) + vAxisFroms[ uiA ][ uiBestDim ],
                                          xEnd,
                                          uiAxisSplitPos[ uiNextA ] );
                    vAxisTos[ uiA ][ uiBestDim ] = xRet != xEnd ? xRet - vAxisCoordinates[ uiBestDim ].begin( )
                                                                : vAxisCoordinates[ uiBestDim ].size( );
                }

                if constexpr( EXPLAIN_QUERY )
                {
                    std::cerr << "\taxis from " << vAxisFroms[ uiA ] << " to " << vAxisTos[ uiA ] << std::endl;
                    std::cerr << "\tcoords: " << vAxisCoordinates[ uiBestDim ][ vAxisFroms[ uiA ][ uiBestDim ] ] << " "
                              << vAxisCoordinates[ uiBestDim ][ vAxisTos[ uiA ][ uiBestDim ] - 1 ] << std::endl;
                    std::cerr << "\tat coordinate " << uiAxisSplitPos[ uiA ] << " containing " << uiCountInSplit[ uiA ]
                              << " points" << std::endl;
                }
                uiLastA = uiA;
                uiA = uiNextA;
            }

            // 3) insert new kd-tree nodes

            offset_t uiMyOffset = vTree.vTree.size( );
            vTree.vTree.push_back( typename OverlayKdTree<type_defs>::OverlayKdTreeNode( uiBestDim ) );
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
                while( uiB < b && uiCountInSplit[ uiB ] == 0 )
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

                std::get<0>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = uiAxisSplitPos[ uiA ];

                generateForPointsHelper(
                    xDatasetId,
                    [ & ]( size_t uiOffsetNode, bool bIsLeaf ) {
                        std::get<1>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = uiOffsetNode;
                        std::get<2>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = bIsLeaf;

                        if( !bRegistered )
                            fRegisterMeInTree( uiMyOffset, false );
                        bRegistered = true;
                    },
                    uiMaxPointsPerOverlay, uiFrom, uiFrom + uiCountInSplit[ uiA ], vOverlayBottomLeft, vOverlayTopRight,
                    vAxisCoordinates, vAxisFroms[ uiA ], vAxisTos[ uiA ], xProg, rProfiler );

                uiFrom += uiCountInSplit[ uiA ];
            }
            if( !bRegistered )
                fRegisterMeInTree( uiMyOffset, false );
            assert( uiFrom == uiTo );

            if constexpr( EXPLAIN_QUERY )
                std::cerr << "split in dimension " << uiBestDim << " done. (" << uiMyOffset << ": "
                          << (size_t)vTree.vTree[ uiMyOffset ].uiSplitDimension << ")" << std::endl;
        }
    }

    struct CompOverlay
    {
        bool operator( )( coordinate_t uiPos, const std::tuple<coordinate_t, offset_t, bool>& rEle )
        {
            return uiPos < std::get<0>( rEle );
        }
    };

    template <bool SILENT>
    std::tuple<pos_t, pos_t, const overlay_meta_t*> getOverlay( const class_key_t& xDatasetId, const pos_t& vPos ) const
    {
        auto cIter =
            std::lower_bound( vTree.vRoots.cbegin( ), vTree.vRoots.cend( ), xDatasetId,
                              [ & ]( const std::tuple<class_key_t, offset_t, bool>& rA,
                                     const class_key_t& xDatasetId ) { return std::get<0>( rA ) < xDatasetId; } );

        if( cIter == vTree.vRoots.cend( ) || std::get<0>( *cIter ) != xDatasetId )
            throw std::runtime_error( "getOverlay: dataset Id not found" );

        pos_t vBottomLeft{ };
        pos_t vTopRight{ };
        for( size_t uiI = 0; uiI < d; uiI++ )
            vTopRight[uiI] = std::numeric_limits<coordinate_t>::max();
        size_t uiCurr = std::get<1>( *cIter );
        bool bIsLeaf = std::get<2>( *cIter );
        while( !bIsLeaf )
        {
            assert( uiCurr < vTree.vTree.size( ) );
            auto xIt = std::upper_bound( vTree.vTree[ uiCurr ].vChildren.begin( ),
                                         vTree.vTree[ uiCurr ].vChildren.end( ),
                                         vPos[ vTree.vTree[ uiCurr ].uiSplitDimension ],
                                         CompOverlay( ) );
            if( xIt == vTree.vTree[ uiCurr ].vChildren.begin( ) )
            {
                vTopRight[vTree.vTree[ uiCurr ].uiSplitDimension] = std::get<0>( *(xIt+1) );
                return std::make_tuple( vBottomLeft, vTopRight, nullptr );
            }
            else
            {
                if( xIt != vTree.vTree[ uiCurr ].vChildren.end( ) )
                    vTopRight[vTree.vTree[ uiCurr ].uiSplitDimension] = std::get<0>( *xIt );
                --xIt;
                vBottomLeft[vTree.vTree[ uiCurr ].uiSplitDimension] = std::get<0>( *xIt );
                if( std::get<1>( *xIt ) == std::numeric_limits<offset_t>::max( ) )
                    return std::make_tuple( vBottomLeft, vTopRight, nullptr );
                bIsLeaf = std::get<2>( *xIt );
                uiCurr = std::get<1>( *xIt );
            }
        }
        assert( uiCurr < vTree.vLeaves.size( ) );
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\toverlay index: " << uiCurr << std::endl;
        return std::make_tuple( vBottomLeft, vTopRight, &vTree.vLeaves[ uiCurr ] );
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

        auto cIter =
            std::lower_bound( vTree.vRoots.cbegin( ), vTree.vRoots.cend( ), xDatasetId,
                              [ & ]( const std::tuple<class_key_t, offset_t, bool>& rA,
                                     const class_key_t& xDatasetId ) { return std::get<0>( rA ) < xDatasetId; } );

        if( cIter == vTree.vRoots.cend( ) || std::get<0>( *cIter ) != xDatasetId )
            throw std::runtime_error( "getOverlay: dataset Id not found" );

        iterateOverlaysIn( fDo, std::get<1>( *cIter ), std::get<2>( *cIter ), vFrom, vTo );
    }


    template <bool SILENT> data_t countToZero( const overlay_meta_t& rO, pos_t vPos ) const
    {
        data_t uiRet{ };
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            if constexpr( EXPLAIN_QUERY && !SILENT )
                std::cerr << "\t\toverlay value: dimension=" << uiI << " position=" << vPos[ uiI ] << std::endl;
            data_t uiVal = vEntries.get( vPos[ uiI ], rO.vEntryBegins[ uiI ], rO.vSizes[ uiI ] );
            if constexpr( EXPLAIN_QUERY && !SILENT )
                std::cerr << "\t\toverlay value=" << uiVal << std::endl;
            uiRet += uiVal;
        }

        data_t uiVal{ };
        vPoints.forRange(
            [ & ]( const point_t& rP ) {
                bool bToTopRight = true;
                for( size_t uiI = 0; uiI < d; uiI++ )
                    bToTopRight = bToTopRight && ( rP.vPos[ uiI ] > vPos[ uiI ] );
                if( bToTopRight )
                    uiVal[ rP.uiLayer ] += 1;
                return true;
            },
            rO.uiPointsBegin, rO.uiPointsEnd );
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\toverlay number of points: " << uiVal << " of total: " << rO.uiPointsEnd - rO.uiPointsBegin
                      << std::endl;
        uiRet += uiVal;

        uiVal = vEntries.getLast( rO.vEntryBegins[ 0 ], rO.vSizes[ 0 ] );
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\toverlay top right corner: value= -" << uiVal << std::endl;
        assert( uiRet >= uiVal ); // can never have a negative return number
        uiRet -= uiVal;
        return uiRet;
    }

    template <bool SILENT> data_t countToZero( class_key_t xDatasetId, pos_t vPos ) const
    {
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\tcountToZero" << std::endl;
        if( vTree.vRoots.size( ) == 0 )
            return data_t{ }; // return zero
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            if( vPos[ uiI ] == 0 )
            {
                if constexpr( EXPLAIN_QUERY && !SILENT )
                    std::cerr << "\t\tdimension: " << uiI << " is zero. result must be zero as well." << std::endl;
                return data_t{ }; // return zero
            }
            else
                --vPos[ uiI ];
        }
        std::tuple<pos_t, pos_t, const overlay_meta_t*> rO = getOverlay<SILENT>( xDatasetId, vPos );
        if( std::get<2>(rO) == nullptr )
        {
            if constexpr( EXPLAIN_QUERY && !SILENT )
                std::cerr << "\t\tpoint is to the bottom left of most bottom left overlay. result must be zero."
                          << std::endl;
            return data_t{ }; // return zero
        }
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\toverlay: " << rO << std::endl;
        return countToZero<SILENT>( *std::get<2>(rO), vPos );
    }

    template <bool SILENT>
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
    Tree( std::string sPrefix, bool bWrite = false ) : 
        vTree( sPrefix, bWrite ), vEntries( sPrefix, bWrite ), vPoints( sPrefix, bWrite ), vDesc( sPrefix, bWrite )
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

    void addPoint( class_key_t, pos_t vPos, layers_t uiLayer, std::string sDesc )
    {
        assert(uiLayer < LAYERS);
        vPoints.add( vPos, vDesc.add( sDesc ), uiLayer );
    }

    void generate( class_key_t xDatasetId, size_t uiMaxPointsPerOverlay, size_t uiFrom, size_t uiTo,
                   std::optional<progress_stream_t> xProg = { } )
    {
        std::array<coordinate_t, d> vOverlayBottomLeft{ };
        std::array<coordinate_t, d> vOverlayTopRight{ };
        std::array<cord_vec_t, d> vAxisCoordinates{ };
        std::array<size_t, d> vFrom{ };
        std::array<size_t, d> vTo{ };
        
        Profiler xProfiler("generate axis coords");

        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            xProg << CLRLN << "Generating axis coordinates " << uiI << " of " << d << "...";
            vOverlayTopRight[ uiI ] = std::numeric_limits<coordinate_t>::max( );

            vPoints.sortByDim( uiI, uiFrom, uiTo );
            vAxisCoordinates[ uiI ] = cord_vec_generator.vec( );
            vPoints.forRange(
                [ & ]( const point_t& rP2 ) {
                    if( vAxisCoordinates[ uiI ].size( ) == 0 || vAxisCoordinates[ uiI ].back( ) < rP2.vPos[ uiI ] )
                        vAxisCoordinates[ uiI ].push_back( rP2.vPos[ uiI ] );
                    return true;
                },
                uiFrom, uiTo );
            vTo[ uiI ] = vAxisCoordinates[ uiI ].size( );
        }

        xProfiler.print("overlays");
        xProg << CLRLN << "Generating overlays...";
        generateForPointsHelper(
            xDatasetId, //
            [ & ]( size_t uiOffsetNode, bool bIsLeaf ) {
                assert( vTree.vRoots.size( ) == 0 || std::get<0>( vTree.vRoots.back( ) ) < xDatasetId );
                vTree.vRoots.push_back( std::make_tuple( xDatasetId, uiOffsetNode, bIsLeaf ) );
            },
            uiMaxPointsPerOverlay, //
            uiFrom, //
            uiTo, //
            vOverlayBottomLeft, //
            vOverlayTopRight, //
            vAxisCoordinates, //
            vFrom, //
            vTo, //
            xProg,
            xProfiler );
        xProfiler.print("");
        xProg << CLRLN << "Done.\n";
    }

    template <bool SILENT> data_t count( class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        pos_t vCurr{ };
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "counting query for datsetId=" << xDatasetId << " bottomLeftCorner=" << vFrom
                      << " topRightCorner=" << vTo << std::endl;
        return countHelper<SILENT>( xDatasetId, vFrom, vTo, vCurr, 0, 0 );
    }

    void iterate( std::function<void( pos_t, layers_t, std::string )> fDo, const class_key_t& xDatasetId, pos_t vFrom,
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

    using vec_return_get_t = std::vector<std::pair<pos_t, std::string>>;
    using return_get_t = std::array<vec_return_get_t, LAYERS>;

    return_get_t get( class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        return_get_t vRet {};
        //for( size_t uiI = 0; uiI < LAYERS; uiI++ )
        //    vRet[ uiI ] = std::make_shared<vec_return_get_t>( );

        iterate(
            [ & ]( pos_t vPos, layers_t uiLayer, std::string sDesc ) { assert(uiLayer < vRet.size()); vRet[ uiLayer ].emplace_back( vPos, sDesc ); },
            xDatasetId, vFrom, vTo );
        return vRet;
    }

    std::string print( )
    {
        std::stringstream ss;
        ss << *this << std::endl;
        return ss.str( );
    }

    size_t size( ) const
    {
        return vPoints.size( );
    }
};


} // namespace kdpstree

namespace std
{
template <typename type_defs> ostream& operator<<( ostream& os, const typename kdpstree::Tree<type_defs>& rTree )
{
    os << "Tree:" << std::endl
       << rTree.vTree << "Entries:" << std::endl
       << rTree.vEntries << "Points:" << std::endl
       << rTree.vPoints << "Desc:" << std::endl
       << rTree.vDesc << std::endl;
    return os;
}
} // namespace std


#if WITH_PYTHON
template <typename type_defs> void exportStream( pybind11::module& m, std::string sName )
{
    pybind11::class_<typename type_defs::progress_stream_t>( m, sName.c_str( ) );
}
template <typename type_defs> void exportTree( pybind11::module& m, std::string sName )
{
    //pybind11::bind_vector<typename kdpstree::Tree<type_defs>::vec_return_get_t,
    //                      std::shared_ptr<typename kdpstree::Tree<type_defs>::vec_return_get_t>>(
    //    m, ( "__" + sName + "_vec_return_get" ).c_str( ) );
    pybind11::class_<kdpstree::Tree<type_defs>>( m, sName.c_str( ) )
        .def( pybind11::init<std::string, bool>( ), pybind11::arg( "path" ), pybind11::arg( "write_mode" ) = false ) // constructor
        .def( "add_point", &kdpstree::Tree<type_defs>::addPoint )
        .def( "num_points", &kdpstree::Tree<type_defs>::numPoints )
        .def( "generate", &kdpstree::Tree<type_defs>::generate, pybind11::arg( "xDatasetId" ),
              pybind11::arg( "uiMaxPointsPerOverlay" ), pybind11::arg( "uiFrom" ), pybind11::arg( "uiTo" ),
              pybind11::arg( "xProg" ) = typename type_defs::progress_stream_t( ) )
        // pybind11::arg( "xProg" ) = std::optional<typename type_defs::progress_stream_t>( ) )
        .def( "count", &kdpstree::Tree<type_defs>::template count<false>, "" )
        .def( "get", &kdpstree::Tree<type_defs>::get, "" )
        .def( "__str__", &kdpstree::Tree<type_defs>::print )
        .def( "__len__", &kdpstree::Tree<type_defs>::size )
        .def( "clear", &kdpstree::Tree<type_defs>::clear )

        ;
}
#endif