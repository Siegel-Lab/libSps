#pragma once

#include "cstree/type_defs.h"

namespace cstree
{
template <typename type_defs> class DataPoint
{
  public:
    typename type_defs::point_t vPos;
    typename type_defs::points_vec_offset_t uiDesc;

    DataPoint( typename type_defs::point_t vPos, typename type_defs::points_vec_offset_t uiDesc )
        : vPos( vPos ), uiDesc( uiDesc )
    {}
};
} // namespace cstree