#pragma once

#include "bps/overlay_meta.h"
#include "bps/points.h"
#include "bps/type_defs.h"


namespace bps
{

template <typename type_defs> class Overlays
{
    using val_t = typename type_defs::val_t;
    using coordinate_t = typename type_defs::coordinate_t;
    static const coordinate_t d = type_defs::d;
    using pos_t = typename type_defs::pos_t;
    using point = Point<type_defs>;
    using map_generator = typename type_defs::map_generator;
    using overlay_key_t = typename type_defs::overlay_key_t;
    using overlay_meta_t = OverlayMeta<type_defs>;
    auto overlay_meta_map_generator = map_generator<overlay_key_t, overlay_meta>( );
    auto axis_vec_generator = vec_generator<coordinate_t>( );
    auto prefix_sum_vec_generator = vec_generator<val_t>( );

    using overlay_meta_map = typeof( overlay_meta_map_generator( "" ) );
    using overlay_entries = OverlayEntries<type_defs>;
    using points = Points<type_defs>;

    overlay_meta_map vMeta;
    overlay_entries vEntries;

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
    size_t constructEytzinger( axis_vec_t& vAxisVec, size_t uiJ, size_t uiI, size_t uiK, size_t uiN )
    {
        if( uiK <= n )
        {
            uiI = constructEytzinger( vAxisVecs, uiJ, uiI, uiK * 2, uiN );
            vEntries.vData[ uiK - 1 + uiJ ].first = vAxisVec[ uiI++ ];
            uiI = constructEytzinger( vAxisVecs, uiJ, uiI, uiK * 2 + 1, uiN );
        }
        return uiI;
    }

    template <typename axis_vecs_t, typename prefix_sums_vecs_t>
    std::array<val_t, d>
    generateForPoints( class_key_t& rDatasetId, size_t uiMaxPointsPerOverlay, points& vPoints, size_t uiFrom,
                       size_t uiTo, const axis_vecs_t& vvAxisVecs, std::array<size_t, d> vvAxisVecsIntervalStart,
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
                    [ & ]( const point& rP ) {
                        while( uiX < vvAxisVecsIntervalEnd[ uiI ] vvAxisVecs[ uiI ][ uiX ] <= rP.vPos[ uiI ] )
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
            vMeta.push_back( std::pair<overlay_key_t, overlay_meta_t>(
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
                    [ & ]( const point& rP ) {
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
            size_t uiAxisSplitIndex = ( vvAxisVecsIntervalEnd[ uiI ] + vvAxisVecsIntervalStart[ uiI ] ) / 2;
            coordinate_t uiSplitPos = vvAxisVecs[ uiI ][ uiAxisSplitIndex ];
            size_t uiPointSplitIndex = uiFrom;
            vPoints.forRange(
                [ & ]( const point& rP ) {
                    if( rP.vPos[ uiI ] <= uiSplitPos )
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

            std::array<val_t, d> vFinalPrefixSumLeft = generateForPoints(
                rDatasetId, uiMaxPointsPerOverlay, vPoints, uiFrom, uiPointSplitIndex, vvAxisVecs,
                vvAxisVecsIntervalStart, vvAxisVecsIntervalEndNew, vvPrefixSumVecs, vInitialPrefixSum );
            vInitialPrefixSum[ uiBestDim ] = vFinalPrefixSumLeft[ uiBestDim ];
            std::array<val_t, d> vFinalPrefixSumRight = generateForPoints(
                rDatasetId, uiMaxPointsPerOverlay, vPoints, uiPointSplitIndex, uiTo, vvAxisVecs,
                vvAxisVecsIntervalStartNew, vvAxisVecsIntervalEnd, vvPrefixSumVecs, vInitialPrefixSum );

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

  public:
    Overlays( std::string sPrefix )
        : vMeta( overlay_meta_map_generator( sPrefix + ".overlay_pos" ) ), vEntries( sPrefix )
    {}


    void generateForPoints( class_key_t xDatasetId, size_t uiMaxPointsPerOverlay, points& vPoints, size_t uiFrom,
                            size_t uiTo )
    {
        std::array<typeof( axis_vec_generator( ) ), d> vvAxisVecs = axis_vec_generator( );
        std::array<typeof( prefix_sum_vec_generator( ) ), d> vvPrefixSumVecs = prefix_sum_vec_generator( );
        std::array<size_t, d> vvAxisVecsIntervalStart;
        std::array<size_t, d> vvAxisVecsIntervalEnd;
        std::array<val_t, d> vInitialPrefixSum;
        for( size_t uiI = 0; uiI < d; uiI++ )
        {
            vInitialPrefixSum[ uiI ] = 0;
            vPoints.sortByDim( uiI, 0, vPoints.size( ) );
            vPoints.forRange(
                [ & ]( point& xPoint ) {
                    if( vvAxisVecs[ uiI ].size( ) == 0 || vvAxisVecs[ uiI ].back( ).vPos[ uiI ] != xPoint.vPos[ uiI ] )
                    {
                        vvAxisVecs[ uiI ].push_back( xPoint.vPos[ uiI ] );
                        vvPrefixSumVecs[ uiI ].push_back( 0 );
                    }
                },
                0, vPoints.size( ) );
            vvAxisVecsIntervalStart[ uiI ] = 0;
            vvAxisVecsIntervalEnd[ uiI ] = vvAxisVecs[ uiI ].size( );
        }
        std::array<val_t, d> vFinalPrefixSum =
            generateForPoints( xDatasetId, uiMaxPointsPerOverlay, vPoints, uiFrom, uiTo, vvAxisVecs,
                               vvAxisVecsIntervalStart, vvAxisVecsIntervalEnd, vvPrefixSumVecs, vInitialPrefixSum );
        for( size_t uiI = 0; uiI < d; uiI++ )
            assert( vFinalPrefixSum[ uiI ] == ( uiTo - uiFrom ) );
    }
};

} // namespace bps