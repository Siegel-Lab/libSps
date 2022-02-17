#pragma once

#include "cstree/subtree_info.h"
#include "cstree/type_defs.h"


namespace cstree
{
template <typename type_defs> class ContSum
{
  public:
    typename type_defs::cont_sum_val_t uiVal;
    SubtreeInfo<type_defs> xSubtree;

    ContSum( typename type_defs::cont_sum_val_t uiVal ) : uiVal( uiVal ), xSubtree( )
    {}

    ContSum( ) : ContSum( 0 )
    {}
};
} // namespace cstree