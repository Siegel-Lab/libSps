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
    layers_t uiLayer;

    Point( pos_t vPos, size_t uiDescOffset, layers_t uiLayer )
        : vPos( vPos ), uiDescOffset( uiDescOffset ), uiLayer( uiLayer )
    {}

    Point( ) : vPos{ }, uiDescOffset( 0 ), uiLayer( 0 )
    {}
};


} // namespace kdpstree

namespace std
{

template <typename type_defs>
std::ostream& operator<<( std::ostream& os, const typename kdpstree::Point<type_defs>& xPoint )
{
    os << xPoint.vPos << " l" << (size_t)xPoint.uiLayer << " d" << xPoint.uiDescOffset;
    return os;
}

} // namespace std