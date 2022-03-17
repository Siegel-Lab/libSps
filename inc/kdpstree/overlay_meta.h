#pragma once

#include "kdpstree/overlay_entries.h"
#include "kdpstree/point.h"
#include "kdpstree/points.h"
#include "kdpstree/type_defs.h"


namespace kdpstree
{

template <typename type_defs> class UnalignedOverlayMeta
{
    EXTRACT_TYPE_DEFS; // macro call

    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;
    using overlay_entries_t = OverlayEntries<type_defs>;

  public:
    std::array<offset_t, d> vEntryBegins;
    std::array<offset_t, d> vSizes;
    offset_t uiPointsBegin;
    offset_t uiPointsEnd;

    UnalignedOverlayMeta( std::array<size_t, d> vEntryBegins, //
                          std::array<size_t, d>
                              vSizes, //
                          size_t uiPointsBegin, //
                          size_t uiPointsEnd )
        : vEntryBegins( vEntryBegins ), //
          vSizes( vSizes ), //
          uiPointsBegin( uiPointsBegin ), //
          uiPointsEnd( uiPointsEnd )
    {}

    UnalignedOverlayMeta( )
        : vEntryBegins{ }, //
          vSizes{ }, //
          uiPointsBegin( 0 ), //
          uiPointsEnd( 0 )
    {}
};

template <typename type_defs> class OverlayMeta : public UnalignedOverlayMeta<type_defs>
{
  public:
    using UnalignedOverlayMeta<type_defs>::UnalignedOverlayMeta;

  private:
    static constexpr size_t ALIGN_TO = 64;

    static_assert( sizeof( UnalignedOverlayMeta<type_defs> ) <= ALIGN_TO );
    // make sure sizeof(this) % 4096 == 0
    std::array<char, ALIGN_TO % sizeof( UnalignedOverlayMeta<type_defs> )> __buffer;
}; // struct


} // namespace kdpstree

namespace std
{


template <typename type_defs> std::ostream& operator<<( ostream& os, const kdpstree::OverlayMeta<type_defs>& rMeta )
{
    os << "p" << rMeta.uiPointsBegin << " - p" << rMeta.uiPointsEnd;
    for( size_t uiI = 0; uiI < type_defs::d; uiI++ )
        os << ", (d" << uiI << ": e" << rMeta.vEntryBegins[ uiI ] << ", s" << rMeta.vSizes[ uiI ] << ")";
    return os;
}

} // namespace std