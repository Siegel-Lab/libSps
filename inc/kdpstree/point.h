#pragma once

#include "kdpstree/type_defs.h"

namespace kdpstree
{
template <typename type_defs> class Point
{
    using pos_t = typename type_defs::pos_t;

  public:
    pos_t vFrom;
    pos_t vTo;
    size_t uiDescOffset;
    size_t uiLayer;

    Point( pos_t vFrom, pos_t vTo, size_t uiDescOffset, size_t uiLayer ) : vFrom( vFrom ), vTo( vTo ), 
                                                                           uiDescOffset( uiDescOffset ), 
                                                                           uiLayer(uiLayer)
    {}
};


std::ostream& operator<<(std::ostream& os, const Point& xPoint)
{
    os << xPoint.vFrom << " - " << xPoint.vTo << " l" << xPoint.uiLayer << " d" << xPoint.uiDescOffset;
    return os;
}

} // namespace kdpstree