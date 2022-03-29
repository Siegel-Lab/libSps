#pragma once

#include "sps/type_defs.h"

namespace sps
{
template <typename type_defs> class Point
{
    EXTRACT_TYPE_DEFS; // macro call
    
    using desc_t = Desc<type_defs>;

  public:
    pos_t vPos;
    size_t uiDescOffset;

    Point( pos_t vPos, size_t uiDescOffset )
        : vPos( vPos ), uiDescOffset( uiDescOffset )
    {}

    Point( ) : vPos{ }, uiDescOffset( 0 )
    {}

    friend std::ostream& operator<<( std::ostream& os, const Point& xPoint )
    {
        os << xPoint.vPos << " d" << xPoint.uiDescOffset;
        return os;
    }

    std::ostream& stream( std::ostream& os, const desc_t& vDesc ) const
    {
        os << vPos << ": " << vDesc.get(uiDescOffset);
        return os;
    }
};


} // namespace sps
