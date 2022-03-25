#pragma once

#include "kdpstree/type_defs.h"

namespace kdpstree
{
template <typename type_defs> class Point
{
    EXTRACT_TYPE_DEFS; // macro call

  public:
    pos_t vPos;
    size_t uiDescOffset;

    Point( pos_t vPos, size_t uiDescOffset )
        : vPos( vPos ), uiDescOffset( uiDescOffset )
    {}

    Point( ) : vPos{ }, uiDescOffset( 0 )
    {}
};


} // namespace kdpstree

namespace std
{

template <typename type_defs>
std::ostream& operator<<( std::ostream& os, const typename kdpstree::Point<type_defs>& xPoint )
{
    os << xPoint.vPos << " d" << xPoint.uiDescOffset;
    return os;
}

} // namespace std