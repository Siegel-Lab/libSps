#include "kdpstree/default.h"

PYBIND11_MODULE( libKdpsTree, m )
{
    exportTree<InMemTypeDef<1>>( m, "KdpsTree_1L" );
}