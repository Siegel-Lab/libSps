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

  public:
    const std::array<size_t, d> vEntryBegins;
    const std::array<size_t, d> vSizes;
    const size_t uiPointsBegin;
    const size_t uiPointsEnd;

    OverlayMeta( std::array<size_t, d> vEntryBegins, //
                 std::array<size_t, d>
                     vSizes, //
                 size_t uiPointsBegin, //
                 size_t uiPointsEnd )
        : vEntryBegins( vEntryBegins ), //
          vSizes( vSizes ), //
          uiPointsBegin( uiPointsBegin ), //
          uiPointsEnd( uiPointsEnd )
    {
    }

    std::string print( ) const
    {
        std::string sRet = "p" + std::to_string( uiPointsBegin ) + " - p" + std::to_string( uiPointsEnd );
        for( size_t uiI = 0; uiI < d; uiI++ )
            sRet += ", ( d" + std::to_string( uiI ) + ": e" + std::to_string( vEntryBegins[ uiI ] ) + ", s" +
                    std::to_string( vSizes[ uiI ] ) + ")";
        return sRet;
    }
};

} // namespace kdpstree