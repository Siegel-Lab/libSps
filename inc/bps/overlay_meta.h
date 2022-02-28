#pragma once

#include "bps/type_defs.h"

namespace bps
{

template <typename type_defs> class OverlayMeta
{
    using val_t = typename type_defs::val_t;
    using coordinate_t = typename type_defs::coordinate_t;
    using point = Point<type_defs>;
    using points = Points<type_defs>;
    using overlay_entries = OverlayEntries<type_defs>;
    static const coordinate_t d = type_defs::d;

    using prefix_sum_vec = typeof( vec_generator<val_t>( )( ) );

    using axis_vec = std::array<typeof( vec_generator<coordinate_t>( )( ) ), d>;

    static const coordinate_t d = type_defs::d;
    std::array<size_t, d> vEntryBegins;
    std::array<size_t, d> vSizes;
    size_t uiPointsBegin;
    size_t uiPointsEnd;

  public:
    OverlayMeta( std::array<size_t, d> vEntryBegins, std::array<size_t, d> vSizes, size_t uiPointsBegin,
                 size_t uiPointsEnd, const points& vPoints, overlay_entries& vEntries, const axis_vec& vvAxisVec,
                 std::array<size_t, d> vvAxisVecsIntervalStart, std::array<size_t, d> vvAxisVecsIntervalEnd,
                 std::array<prefix_sum_vec, d>& vvPrefixSumVec, std::array<val_t, d>& vInitialPrefixSum )
        : vEntryBegins( vEntryBegins ), vSizes( vSizes ), uiPointsBegin( uiPointsBegin ), uiPointsEnd( uiPointsEnd )
    {
        // this will get overridded in a sec but we can conviniently use it as intermediate memory
        vPoints.forRange(
            [ & ]( const point& xPoint ) {
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

} // namespace bps