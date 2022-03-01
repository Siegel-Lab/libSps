#pragma once

#include "kdpstree/overlay_meta.h"
#include "kdpstree/point.h"
#include "kdpstree/points.h"
#include "kdpstree/type_defs.h"


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
    using overlay_entries_t = OverlayEntries<type_defs>;
    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;
    using desc_t = Desc<type_defs>;

    map_generator_t<overlay_key_t, overlay_meta_t> overlay_meta_map_generator =
        map_generator_t<overlay_key_t, overlay_meta_t>( );
    using overlay_meta_map_t = typeof( overlay_meta_map_generator( "" ) );

    vec_generator_t<coordinate_t> axis_vec_generator = vec_generator_t<coordinate_t>( );
    using axis_vec_t = typeof( axis_vec_generator( "" ) );

    vec_generator_t<val_t> prefix_sum_vec_generator = vec_generator_t<val_t>( );
    using prefix_sum_vec_t = typeof( prefix_sum_vec_generator( "" ) );


    overlay_entries_t vEntries;
    overlay_meta_map_t vMeta;
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
            vEntries.vData[ uiK - 1 + uiJ ].first = vAxisVecs[ uiI++ ];
            uiI = constructEytzinger( vAxisVecs, uiJ, uiI, uiK * 2 + 1, uiN );
        }
        return uiI;
    }

    template <typename axis_vecs_t, typename prefix_sums_vecs_t>
    std::array<val_t, d>
    generateForPointsHelper( class_key_t& rDatasetId, size_t uiMaxPointsPerOverlay, size_t uiFrom, size_t uiTo,
                             const axis_vecs_t& vvAxisVecs, std::array<size_t, d> vvAxisVecsIntervalStart,
                             std::array<size_t, d> vvAxisVecsIntervalEnd, prefix_sums_vecs_t& vvPrefixSumVecs,
                             std::array<val_t, d> vInitialPrefixSum )
    {
        if( uiTo - uiFrom <= uiMaxPointsPerOverlay )
        {
            // recursion termination
            std::array<size_t, d> vBegins;
            std::array<size_t, d> vSizes;
            std::array<coordinate_t, d> vOverlayPos;

            // for each dimension
            for( size_t uiI = 0; uiI < d; uiI++ )
            {
                // 1) construct compressed axis of overlay
                auto vTmp = axis_vec_generator( );
                size_t uiX = vvAxisVecsIntervalStart[ uiI ];
                vPoints.forRange(
                    [ & ]( const point_t& rP ) {
                        while( uiX < vvAxisVecsIntervalEnd[ uiI ] && vvAxisVecs[ uiI ][ uiX ] <= rP.vPos[ uiI ] )
                        {
                            if( vTmp.size( ) == 0 || vvAxisVecs[ uiI ][ uiX ] != vTmp.back( ) )
                                vTmp.push_back( vvAxisVecs[ uiI ][ uiX ] );
                            ++uiX;
                        }
                        if( vTmp.size( ) == 0 || rP.vPos[ uiI ] != vTmp.back( ) )
                            vTmp.push_back( rP.vPos[ uiI ] );
                    },
                    uiFrom, uiTo );

                while( uiX < vvAxisVecsIntervalEnd[ uiI ] )
                {
                    if( vTmp.size( ) == 0 || vvAxisVecs[ uiI ][ uiX ] != vTmp.back( ) )
                        vTmp.push_back( vvAxisVecs[ uiI ][ uiX ] );
                    ++uiX;
                }

                // 2) set bottom left position of overlay
                vOverlayPos[ uiI ] = vvAxisVecs[ uiI ][ vvAxisVecsIntervalStart[ uiI ] ];

                // 3) save pointers to the constructed axes
                vBegins[ uiI ] = vEntries.size( );
                vSizes[ uiI ] = vTmp.size( );
                vEntries.incSize( vTmp.size( ) );

                // 4) construct eytzinger representation of sorted compressed axis
                constructEytzinger( vTmp, vBegins[ uiI ], 0, 1, vTmp.size( ) );
            } // for

            // insert overaly into correct position in tree
            vMeta.insert( std::pair<overlay_key_t, overlay_meta_t>(
                overlay_key_t( rDatasetId, vOverlayPos ),
                // create overlay
                overlay_meta_t( vBegins, vSizes, uiFrom, uiTo, vPoints, vEntries, vvAxisVecs, vvAxisVecsIntervalStart,
                                vvAxisVecsIntervalEnd, vvPrefixSumVecs, vInitialPrefixSum ) ) );
            return vInitialPrefixSum;
        }
        else
        {
            // recursive call: split points in half

            // 1) find dimension that splits points most evenly in half
            size_t uiBestDim = 0;
            double fBestRatio = 1;
            for( size_t uiI = 0; uiI < d; uiI++ )
            {
                size_t uiAxisSplitIndex = ( vvAxisVecsIntervalEnd[ uiI ] + vvAxisVecsIntervalStart[ uiI ] ) / 2;
                coordinate_t uiSplitPos = vvAxisVecs[ uiI ][ uiAxisSplitIndex ];
                size_t uiNumBefore = 0;
                vPoints.forRange(
                    [ & ]( const point_t& rP ) {
                        if( rP.vPos[ uiI ] <= uiSplitPos )
                            uiNumBefore += 1;
                    },
                    uiFrom, uiTo );
                double fRatio = std::abs( 0.5 - ( uiNumBefore / (double)( uiTo - uiFrom ) ) );
                if( fRatio < fBestRatio )
                {
                    uiBestDim = uiI;
                    fBestRatio = fRatio;
                }
            }

            // 2) seperate points in that dimension (by sorting and looking for the split pos)
            vPoints.sortByDim( uiBestDim, uiFrom, uiTo );
            size_t uiAxisSplitIndex = ( vvAxisVecsIntervalEnd[ uiBestDim ] + vvAxisVecsIntervalStart[ uiBestDim ] ) / 2;
            coordinate_t uiSplitPos = vvAxisVecs[ uiBestDim ][ uiAxisSplitIndex ];
            size_t uiPointSplitIndex = uiFrom;
            vPoints.forRange(
                [ & ]( const point_t& rP ) {
                    if( rP.vPos[ uiBestDim ] <= uiSplitPos )
                        uiPointSplitIndex += 1;
                },
                uiFrom, uiTo );
            std::array<size_t, d> vvAxisVecsIntervalEndNew = vvAxisVecsIntervalEnd;
            vvAxisVecsIntervalEndNew[ uiBestDim ] = uiAxisSplitIndex;
            std::array<size_t, d> vvAxisVecsIntervalStartNew = vvAxisVecsIntervalStart;
            vvAxisVecsIntervalStartNew[ uiBestDim ] = uiAxisSplitIndex;


            // 3) recursive calls
            // here it is important to go from left to right so that vvPrefixSumVecs is filled in the correct order
            //

            std::array<val_t, d> vFinalPrefixSumLeft = generateForPointsHelper(
                rDatasetId, uiMaxPointsPerOverlay, uiFrom, uiPointSplitIndex, vvAxisVecs, vvAxisVecsIntervalStart,
                vvAxisVecsIntervalEndNew, vvPrefixSumVecs, vInitialPrefixSum );
            vInitialPrefixSum[ uiBestDim ] = vFinalPrefixSumLeft[ uiBestDim ];
            std::array<val_t, d> vFinalPrefixSumRight = generateForPointsHelper(
                rDatasetId, uiMaxPointsPerOverlay, uiPointSplitIndex, uiTo, vvAxisVecs, vvAxisVecsIntervalStartNew,
                vvAxisVecsIntervalEnd, vvPrefixSumVecs, vInitialPrefixSum );

            std::array<val_t, d> vFinalPrefixSum;
            for( size_t uiI = 0; uiI < d; uiI++ )
            {
                if( uiI == uiBestDim )
                    vFinalPrefixSum[ uiI ] = vFinalPrefixSumRight[ uiI ];
                else
                    vFinalPrefixSum[ uiI ] =
                        vFinalPrefixSumLeft[ uiI ] + vFinalPrefixSumRight[ uiI ] - vInitialPrefixSum[ uiI ];
            }
            return vFinalPrefixSum;
        }
    }

    typename overlay_meta_map_t::const_iterator getOverlay( const class_key_t& xDatasetId, const pos_t& vPos ) const
    {
        for( size_t uiI = 0; uiI < d; uiI++ )
            assert( vPos[ uiI ] > 0 );
        // get the first element that is smaller than vPos
        // since the first point in vMeta must be at 0 in each dimension this element must exist and the -1 is valid in
        // all cases
        return --vMeta.lower_bound( overlay_key_t( xDatasetId, vPos ) );
    }

    bool someDimLarger( const pos_t& x, const pos_t& vThan ) const
    {
        for( size_t uiI = 0; uiI < d; uiI++ )
            if( x[ uiI ] > vThan[ uiI ] )
                return true;
        return false;
    }
    bool allDimLarger( const pos_t& x, const pos_t& vThan ) const
    {
        for( size_t uiI = 0; uiI < d; uiI++ )
            if( x[ uiI ] <= vThan[ uiI ] )
                return false;
        return true;
    }

    void iterateOverlaysIn( std::function<void( typename overlay_meta_map_t::const_iterator )> fDo,
                            const class_key_t& xDatasetId, const pos_t& vFrom, const pos_t& vTo ) const
    {
        const pos_t vCurr = vTo;
        do
        {
            typename overlay_meta_map_t::const_iterator xIt = getOverlay( xDatasetId, vCurr ) + 1;
            do
            {
                if( xIt != vMeta.begin( ) )
                    xIt--;
                fDo( xIt );
            } while( xIt != vMeta.begin( ) && allDimLarger( xIt->first->second, vFrom ) );

            if( !someDimLarger( xIt->first->second, vFrom ) )
                break;

            for( size_t uiI = 0; uiI < d; uiI++ )
            {
                if( xIt->first->second[ uiI ] <= vFrom[ uiI ] )
                    vCurr[ uiI ] = vTo[ uiI ];
                else
                    vCurr[ uiI ] = xIt->first->second[ uiI ];
            }
        } while( true );
    }

    val_t countToZero( class_key_t xDatasetId, pos_t vPos ) const
    {
        if( vMeta.size( ) == 0 )
            return 0;
        typename overlay_meta_map_t::const_iterator xIt = getOverlay( xDatasetId, vPos );

        val_t uiRet = 0;
        for( size_t uiI = 0; uiI < d; uiI++ )
            uiRet += vEntries.get( vPos[ uiI ] - 1, xIt->second.vEntryBegins[ uiI ], xIt->second.vSizes[ uiI ] );
        uiRet -=
            ( d - 1 ) * vEntries.get( xIt->first.second[ 0 ], xIt->second.vEntryBegins[ 0 ], xIt->second.vSizes[ 0 ] );
        vPoints.forRange(
            [ & ]( const point_t& rP ) {
                for( size_t uiI = 0; uiI < d; uiI++ )
                    if( rP.vPos[ uiI ] >= vPos[ uiI ] )
                        return;
                uiRet += 1;
            },
            xIt->second.uiPointsBegin, xIt->second.uiPointsEnd );
        return uiRet;
    }

    val_t countHelper( const class_key_t& xDatasetId, const pos_t& vFrom, const pos_t& vTo, pos_t& vCurr,
                       size_t uiDCurr, size_t uiNumStart ) const
    {
        if( uiDCurr == d )
            return countToZero( xDatasetId, vCurr ) * ( uiNumStart % 2 == 0 ? 1 : -1 );
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
    Tree( std::string sPrefix )
        : vMeta( overlay_meta_map_generator( sPrefix + ".overlay_pos" ) ),
          vEntries( sPrefix ),
          vPoints( sPrefix ),
          vDesc( sPrefix )
    {}

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
        std::array<axis_vec_t, d> vvAxisVecs;
        std::array<prefix_sum_vec_t, d> vvPrefixSumVecs;
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            vvAxisVecs[ uiI ] = axis_vec_generator( );
            vvPrefixSumVecs[ uiI ] = prefix_sum_vec_generator( );
        }

        std::array<size_t, d> vvAxisVecsIntervalStart;
        std::array<size_t, d> vvAxisVecsIntervalEnd;
        std::array<val_t, d> vInitialPrefixSum;
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            vInitialPrefixSum[ uiI ] = 0;
            vPoints.sortByDim( uiI, uiFrom, uiTo );
            vPoints.forRange(
                [ & ]( const point_t& xPoint ) {
                    if( vvAxisVecs[ uiI ].size( ) == 0 || vvAxisVecs[ uiI ].back( ) != xPoint.vPos[ uiI ] )
                    {
                        vvAxisVecs[ uiI ].push_back( xPoint.vPos[ uiI ] );
                        vvPrefixSumVecs[ uiI ].push_back( 0 );
                    }
                },
                0, vPoints.size( ) );
            vvAxisVecsIntervalStart[ uiI ] = 0;
            vvAxisVecsIntervalEnd[ uiI ] = vvAxisVecs[ uiI ].size( );
        }
        std::array<val_t, d> vFinalPrefixSum = generateForPointsHelper( xDatasetId, //
                                                                        uiMaxPointsPerOverlay, //
                                                                        uiFrom, //
                                                                        uiTo, //
                                                                        vvAxisVecs, //
                                                                        vvAxisVecsIntervalStart, //
                                                                        vvAxisVecsIntervalEnd, //
                                                                        vvPrefixSumVecs, //
                                                                        vInitialPrefixSum );
        for( size_t uiI = 0; uiI < d; uiI++ )
            assert( vFinalPrefixSum[ uiI ] == ( uiTo - uiFrom ) );
    }

    val_t count( class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        pos_t vCurr{ };
        return countHelper( xDatasetId, vFrom, vTo, vCurr, 0, 0 );
    }

    void iterate( std::function<void( pos_t, std::string )> fDo, class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        iterateOverlaysIn(
            [ & ]( typename overlay_meta_map_t::const_iterator xIt ) {
                vPoints.forRange(
                    [ & ]( point_t& xPoint ) {
                        for( size_t uiI = 0; uiI < d; uiI++ )
                        {
                            if( xPoint.vPos[ uiI ] >= vTo[ uiI ] )
                                return;
                            if( xPoint.vPos[ uiI ] < vFrom[ uiI ] )
                                return;
                        }
                        fDo( xPoint.vPos, vDesc.get( xPoint.uiDescOffset ) );
                    },
                    xIt->second.uiPointsBegin, xIt->second.uiPointsEnd );
            },
            xDatasetId, vFrom, vTo );
    }

    std::vector<std::pair<pos_t, std::string>> get( class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        std::vector<std::pair<pos_t, std::string>> vRet;
        iterate( [ & ]( pos_t vPos, std::string sDesc ) { vRet.emplace_back( vPos, sDesc ); }, xDatasetId, vFrom, vTo );
        return vRet;
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
        //.def( "__str__", &kdpstree::Tree<type_defs>::print )
        //.def( "str_raw", &kdpstree::Tree<type_defs>::printRaw )

        ;
}
#endif