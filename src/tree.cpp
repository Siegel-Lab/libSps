#include "kdpstree/default.h"

PYBIND11_MODULE( libKdpsTree, m )
{
    exportTree<InMemTypeDef>( m, "KdpsTree" );
}