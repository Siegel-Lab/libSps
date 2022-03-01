#pragma once

#include "kdpstree/overlay_entries.h"
#include "kdpstree/point.h"
#include "kdpstree/points.h"
#include "kdpstree/type_defs.h"


namespace kdpstree
{

template <typename type_defs> class OverlayMeta
{
    EXTRACT_TYPE_DEFS; // macro call

    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;
    using overlay_entries_t = OverlayEntries<type_defs>;

    vec_generator_t<val_t> prefix_sum_vec_generator = vec_generator_t<val_t>( );
    using prefix_sum_vec_t = typeof( prefix_sum_vec_generator( ) );

    using axis_vec_t = std::array<prefix_sum_vec_t, d>;

  public:
    const std::array<size_t, d> vEntryBegins;
    const std::array<size_t, d> vSizes;
    const size_t uiPointsBegin;
    const size_t uiPointsEnd;

    OverlayMeta( std::array<size_t, d> vEntryBegins, //
                 std::array<size_t, d> vSizes, //
                 size_t uiPointsBegin, //
                 size_t uiPointsEnd, //
                 const points_t& vPoints, //
                 overlay_entries_t& vEntries, //
                 const axis_vec_t& vvAxisVec, //
                 std::array<size_t, d> vvAxisVecsIntervalStart, //
                 std::array<size_t, d> vvAxisVecsIntervalEnd, //
                 std::array<prefix_sum_vec_t, d>& vvPrefixSumVec, //
                 std::array<val_t, d>& vInitialPrefixSum )
        : vEntryBegins( vEntryBegins ), //
          vSizes( vSizes ), //
          uiPointsBegin( uiPointsBegin ), //
          uiPointsEnd( uiPointsEnd )
    {
        // this will get overridded in a sec but we can conviniently use it as intermediate memory
        vPoints.forRange(
            [ & ]( const point_t& xPoint ) {
                for( size_t uiI = 0; uiI < d; uiI++ )
                {
                    assert( vEntries.has( xPoint.vPos[ uiI ], vEntryBegins[ uiI ], vSizes[ uiI ] ) );
                    vEntries.variableGet( xPoint.vPos[ uiI ], vEntryBegins[ uiI ], vSizes[ uiI ] ) += 1;
                }
            },
            uiPointsBegin, uiPointsEnd );

        for( size_t uiI = 0; uiI < d; uiI++ )
            for( size_t uiX = vvAxisVecsIntervalStart[ uiI ]; uiX < vvAxisVecsIntervalEnd[ uiI ]; uiX++ )
            {
                if( vEntries.has( vvAxisVec[ uiI ][ uiX ], vEntryBegins[ uiI ], vSizes[ uiI ] ) )
                {
                    vInitialPrefixSum[ uiI ] +=
                        vEntries.get( vvAxisVec[ uiI ][ uiX ], vEntryBegins[ uiI ], vSizes[ uiI ] );
                    vEntries.variableGet( vvAxisVec[ uiI ][ uiX ], vEntryBegins[ uiI ], vSizes[ uiI ] ) =
                        vvPrefixSumVec[ uiI ][ uiX ];
                }
                vvPrefixSumVec[ uiI ][ uiX ] += vInitialPrefixSum[ uiI ];
            }
    }
};

} // namespace kdpstree