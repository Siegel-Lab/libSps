#pragma once

#include "kdpstree/type_defs.h"

namespace kdpstree
{
template <typename type_defs> class Point
{
    using pos_t = typename type_defs::pos_t;

  public:
    pos_t vPos;
    size_t uiDescOffset;

    Point( pos_t vPos, size_t uiDescOffset ) : vPos( vPos ), uiDescOffset( uiDescOffset )
    {}
};
} // namespace kdpstree