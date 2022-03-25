#pragma once

#include "sps/type_defs.h"

namespace sps
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


} // namespace sps

namespace std
{

template <typename type_defs>
std::ostream& operator<<( std::ostream& os, const typename sps::Point<type_defs>& xPoint )
{
    os << xPoint.vPos << " d" << xPoint.uiDescOffset;
    return os;
}

} // namespace std