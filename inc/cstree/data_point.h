#pragma once

#include "cstree/type_defs.h"

namespace cstree
{
template <typename type_defs> class DataPoint
{
  private:
    type_defs::point_t vPos;
    type_defs::points_vec_offset_t uiDesc;
}
} // namespace cstree