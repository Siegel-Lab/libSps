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

    std::tuple<offset_t, bool> generateForPointsHelper( class_key_t xDatasetId, 
                                                        size_t uiMaxPointsPerOverlay, size_t uiFrom, size_t uiTo,
                                                        std::array<coordinate_t, d> vOverlayBottomLeft,
                                                        std::array<coordinate_t, d> vOverlayTopRight )
    {
        if( uiFrom == uiTo )
            return std::make_tuple( std::numeric_limits<offset_t>::max( ), true );
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
            // recursion termination
            std::array<size_t, d> vBegins{};
            std::array<size_t, d> vSizes{};

            // for each dimension
            for( size_t uiI = 0; uiI < d; uiI++ )
            {
                // std::cout << "dim " << uiI << std::endl;
                // 1) construct compressed axis of overlay
                vPoints.sortByDim( uiI, uiFrom, uiTo );
                auto vTmp = axis_vec_generator( );
                val_t uiContSum = 0;
                size_t uiX = uiFrom;
                forRange( [&](coordinate_t uiP, const val_t& uiVal){
                    while(uiX < uiTo && vPoints.vData[uiX].vPos[uiI] <= uiP)
                    {
                        uiContSum ++;
                        if(vTmp.size() == 0 || vTmp.back().first < vPoints.vData[uiX].vPos[uiI])
                            vTmp.emplace_back(vPoints.vData[uiX].vPos[uiI], uiContSum + uiVal);
                        else
                            vTmp.back().second = uiContSum + uiVal;
                        uiX++;
                    }
                    if(vTmp.size() == 0 || vTmp.back().first < uiP)
                        vTmp.emplace_back(uiP, uiContSum + uiVal);
                    else
                        vTmp.back().second = uiContSum + uiVal;
                }, xDatasetId, uiI, vOverlayBottomLeft[uiI], vOverlayTopRight[uiI] );
                assert(uiX == uiTo);

                // 3) save pointers to the constructed axes
                vBegins[ uiI ] = vEntries.size( );
                vSizes[ uiI ] = vTmp.size( );
                vEntries.incSize( vTmp.size( ) );

                if constexpr( EXPLAIN_QUERY )
                {
                    std::cerr << uiI << " inc by " << vTmp.size( ) << std::endl;
                    for( size_t uiK = 0; uiK < vTmp.size( ); uiK++ )
                    {
                        std::cerr << "vTmp " << vTmp[ uiK ].first << ": " << vTmp[ uiK ].second << std::endl;
                    }
                }

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

            return std::make_tuple( uiMyOffset, true );
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
            {
                std::cerr << "splitting in dimension " << uiBestDim << " start "
                          << vPoints.vData[uiFrom].vPos[ uiBestDim ] << " end "
                          << vPoints.vData[ uiTo-1 ].vPos[ uiBestDim ] << " with max count "
                          << uiBestCount << std::endl;
                for( size_t uiA = 0; uiA < b; uiA++ )
                    std::cerr << "\tat coordinate " << uiAxisSplitPos[ uiA ]
                              << " containing " << uiCountInSplit[ uiA ] << " points" << std::endl;
            }

            // 3) insert new kd-tree node
            offset_t uiMyOffset = vTree.vTree.size( );
            vTree.vTree.emplace_back( );
            vTree.vTree[ uiMyOffset ].uiSplitDimension = uiBestDim;
            for( size_t uiA = 0; uiA < b; uiA++ )
            {
                std::get<0>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = uiAxisSplitPos[ uiA ];
                std::get<1>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = std::numeric_limits<offset_t>::max( );
            }

            for( size_t uiA = 0; uiA < b; uiA++ )
            {
                vOverlayBottomLeft[uiBestDim] = uiAxisSplitPos[ uiA ];
                if(uiA + 1 < b)
                    vOverlayTopRight[uiBestDim] = uiAxisSplitPos[ uiA + 1 ];

                // std::cout << "split dim " << uiBestDim << " at " << uiSplitPos << " into " << 100 * fBestRatio << "%"
                //          << std::endl;


                // 4) recursive calls
                // here it is important to go from left to right so that vvPrefixSumVecs is filled in the correct order
                //

                auto xRecRet =
                    generateForPointsHelper( xDatasetId, uiMaxPointsPerOverlay, uiFrom, uiFrom + uiCountInSplit[ uiA ],
                                             vOverlayBottomLeft, vOverlayTopRight );
                

                // set offset and is leaf node for outgoing pointer of tree
                std::get<0>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = std::get<0>( xRecRet );
                std::get<1>( vTree.vTree[ uiMyOffset ].vChildren[ uiA ] ) = std::get<1>( xRecRet );

                uiFrom += uiCountInSplit[ uiA ];
            }
            assert( uiFrom == uiTo );

            return std::make_tuple( uiMyOffset, false );
        }
    }

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
                            vPos[ vTree.vTree[ uiCurr ].uiSplitDimension ] )
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
                if constexpr( EXPLAIN_QUERY )
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
            for( size_t uiA = 0; uiA < b; uiA++ )
                if( std::get<1>( vTree.vTree[ uiCurr ].vChildren[ uiA ] ) != std::numeric_limits<offset_t>::max( ) &&
                    ( uiA + 1 == b || vFrom[ vTree.vTree[ uiCurr ].uiSplitDimension ] <
                                          std::get<0>( vTree.vTree[ uiCurr ].vChildren[ uiA + 1 ] ) ) &&
                    vTo[ vTree.vTree[ uiCurr ].uiSplitDimension ] >
                        std::get<0>( vTree.vTree[ uiCurr ].vChildren[ uiA ] ) )
                    iterateOverlaysIn( fDo, std::get<1>( vTree.vTree[ uiCurr ].vChildren[ uiA ] ),
                                       std::get<2>( vTree.vTree[ uiCurr ].vChildren[ uiA ] ), vFrom, vTo );
    }

    void iterateOverlaysIn( std::function<void( const overlay_meta_t& )> fDo, const class_key_t& xDatasetId,
                            const pos_t& vFrom, const pos_t& vTo ) const
    {
        for( size_t uiI = 0; uiI < vTree.vRoots.size( ); uiI++ )
            if( std::get<0>( vTree.vRoots[ uiI ] ) == xDatasetId )
                iterateOverlaysIn( fDo, std::get<1>( vTree.vRoots[ uiI ] ), std::get<2>( vTree.vRoots[ uiI ] ), vFrom,
                                   vTo );
    }

    void forRange(std::function<void( coordinate_t, const val_t& )> fDo, size_t uiCurr, bool bIsLeaf, size_t uiD, 
                  coordinate_t uiFrom, coordinate_t uiTo)
    {
        if( bIsLeaf )
            vEntries.forRange(fDo, vTree.vLeaves[ uiCurr ].vEntryBegins[uiD], 
                                             vTree.vLeaves[ uiCurr ].vSizes[uiD], uiFrom, uiTo);
        else if(vTree.vTree[ uiCurr ].uiSplitDimension == uiD)
        {
            for( size_t uiA = 0; uiA < b; uiA++ )
                if( std::get<1>( vTree.vTree[ uiCurr ].vChildren[ uiA ] ) != std::numeric_limits<offset_t>::max( ) &&
                    ( uiA + 1 == b || uiFrom < std::get<0>( vTree.vTree[ uiCurr ].vChildren[ uiA + 1 ] ) ) &&
                    uiTo > std::get<0>( vTree.vTree[ uiCurr ].vChildren[ uiA ] ) )
                    forRange( fDo, std::get<1>( vTree.vTree[ uiCurr ].vChildren[ uiA ] ),
                              std::get<2>( vTree.vTree[ uiCurr ].vChildren[ uiA ] ), uiD, uiFrom, uiTo );
        }
        else
            for( size_t uiA = 1; uiA <= b; uiA++ )
                if( std::get<1>( vTree.vTree[ uiCurr ].vChildren[ b-uiA ] ) != std::numeric_limits<offset_t>::max( )  )
                {
                    forRange( fDo, std::get<1>( vTree.vTree[ uiCurr ].vChildren[ b-uiA ] ),
                              std::get<2>( vTree.vTree[ uiCurr ].vChildren[ b-uiA ] ), uiD, uiFrom, uiTo );
                    return;
                }
    }

    void forRange(std::function<void( coordinate_t, const val_t& )> fDo, const class_key_t& xDatasetId, size_t uiD, 
                  coordinate_t uiFrom, coordinate_t uiTo)
    {
        for( size_t uiI = 0; uiI < vTree.vRoots.size( ); uiI++ )
            if( std::get<0>( vTree.vRoots[ uiI ] ) == xDatasetId )
                forRange( fDo, std::get<1>( vTree.vRoots[ uiI ] ), std::get<2>( vTree.vRoots[ uiI ] ), uiD, uiFrom,
                          uiTo);
    }

    val_t countToZero( class_key_t xDatasetId, pos_t vPos ) const
    {
        if( vTree.vRoots.size( ) == 0 )
            return 0;
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            if( vPos[ uiI ] == 0 )
            {
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\t\tdimension: " << uiI << " is zero. result must be zero as well." << std::endl;
                return 0;
            }
            else
                --vPos[ uiI ];
        }
        std::pair<pos_t, const overlay_meta_t*> rO = getOverlay( xDatasetId, vPos );
        if( rO.second == nullptr )
        {
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "\t\tpoint is to the bottom left of most bottom left overlay. result must be zero."
                          << std::endl;
            return 0;
        }
        if constexpr( EXPLAIN_QUERY )
        {
            std::cerr << "\t\toverlay position: (";
            for( size_t uiI = 0; uiI < d; uiI++ )
                std::cerr << rO.first[ uiI ] << ", ";
            std::cerr << ")" << std::endl;
        }

        val_t uiRet = 0;
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "\t\toverlay value: dimension=" << uiI << " position=" << vPos[ uiI ] << std::endl;
            val_t uiVal = vEntries.get( vPos[ uiI ], rO.second->vEntryBegins[ uiI ], rO.second->vSizes[ uiI ] );
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "\t\toverlay value=" << uiVal << std::endl;
            uiRet += uiVal;
        }

        val_t uiVal = rO.second->uiBottomLeftContSum;
        if constexpr( EXPLAIN_QUERY )
            std::cerr << "\t\toverlay bottom left corner: value= " << -( d - 1 ) << " * " << uiVal << std::endl;
        uiRet -= ( d - 1 ) * uiVal;
        uiVal = 0;
        vPoints.forRange(
            [ & ]( const point_t& rP ) {
                for( size_t uiI = 0; uiI < d; uiI++ )
                    if( rP.vPos[ uiI ] > vPos[ uiI ] )
                        return true;
                uiVal += 1;
                return true;
            },
            rO.second->uiPointsBegin, rO.second->uiPointsEnd );
        if constexpr( EXPLAIN_QUERY )
            std::cerr << "\t\toverlay number of points: " << uiVal
                      << " of total: " << rO.second->uiPointsEnd - rO.second->uiPointsBegin << std::endl;
        uiRet += uiVal;
        return uiRet;
    }

    val_t countHelper( const class_key_t& xDatasetId, const pos_t& vFrom, const pos_t& vTo, pos_t& vCurr,
                       size_t uiDCurr, size_t uiNumStart ) const
    {
        if( uiDCurr == d )
        {
            if constexpr( EXPLAIN_QUERY )
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
            auto uiVal = countToZero( xDatasetId, vCurr );
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "\tresult = " << uiVal << std::endl;
            return uiVal * ( uiNumStart % 2 == 0 ? 1 : -1 );
        }
        else
        {
            vCurr[ uiDCurr ] = vFrom[ uiDCurr ];
            val_t uiRet = countHelper( xDatasetId, vFrom, vTo, vCurr, uiDCurr + 1, uiNumStart + 1 );
            vCurr[ uiDCurr ] = vTo[ uiDCurr ];
            uiRet += countHelper( xDatasetId, vFrom, vTo, vCurr, uiDCurr + 1, uiNumStart );
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
        std::array<coordinate_t, d> vOverlayBottomLeft {};
        std::array<coordinate_t, d> vOverlayTopRight {};
        for( size_t uiI = 0; uiI < d; uiI++ )
            vOverlayTopRight[ uiI ] = std::numeric_limits<coordinate_t>::max( );
        auto xRet = generateForPointsHelper( xDatasetId, //
                                             uiMaxPointsPerOverlay, //
                                             uiFrom, //
                                             uiTo, //
                                             vOverlayBottomLeft, //
                                             vOverlayTopRight );
        // for( size_t uiI = 0; uiI < d; uiI++ )
        //    assert( std::get<0>( xRet )[ uiI ] == ( uiTo - uiFrom ) );

        vTree.vRoots.emplace_back( xDatasetId, std::get<0>( xRet ), std::get<1>( xRet ) );
    }

    val_t count( class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        pos_t vCurr{ };
        if constexpr( EXPLAIN_QUERY )
        {
            std::cerr << "counting query for datsetId=" << xDatasetId << " bottomLeftCorner=(";
            for( size_t uiI = 0; uiI < d; uiI++ )
                std::cerr << vFrom[ uiI ] << ", ";
            std::cerr << ") topRightCorner=(";
            for( size_t uiI = 0; uiI < d; uiI++ )
                std::cerr << vTo[ uiI ] << ", ";
            std::cerr << ")" << std::endl;
        }
        return countHelper( xDatasetId, vFrom, vTo, vCurr, 0, 0 );
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
        .def( "count", &kdpstree::Tree<type_defs>::count, "" )
        .def( "get", &kdpstree::Tree<type_defs>::get, "" )
        .def( "__str__", &kdpstree::Tree<type_defs>::print )
        .def( "clear", &kdpstree::Tree<type_defs>::clear )

        ;
}
#endif