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
};


template <typename type_defs>
std::ostream& operator<<(std::ostream& os, const OverlayMeta<type_defs>& rMeta)
{
    os << "p" << rMeta.uiPointsBegin << " - p" << rMeta.uiPointsEnd;
    for( size_t uiI = 0; uiI < type_defs::d; uiI++ )
        os << ", (d" << uiI << ": e" << rMeta.vEntryBegins[ uiI ] << ", s" << rMeta.vSizes[ uiI ] << ")";
    return os;
}


} // namespace kdpstree