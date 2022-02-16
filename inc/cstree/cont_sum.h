#pragma once

#include "cstree/subtree_info.h"
#include "cstree/type_defs.h"


namespace cstree
{
template <typename type_defs> class ContSum
{
  private:
    type_defs::cont_sum_val_t uiVal;
    SubtreeInfo<type_defs> xSubtree;

  public:
    ContSum( type_defs::cont_sum_val_t uiVal ) : uiVal( uiVal ), xSubtree( )
    {}
}
} // namespace cstree