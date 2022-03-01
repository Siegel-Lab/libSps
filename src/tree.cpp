#include "kdpstree/default.h"

PYBIND11_MODULE( libKdpsTree, m )
{
    exportTree<InMemTypeDef<1>>( m, "KdpsTree_1D" );
    exportTree<InMemTypeDef<2>>( m, "KdpsTree_2D" );
    exportTree<InMemTypeDef<3>>( m, "KdpsTree_3D" );
    exportTree<InMemTypeDef<4>>( m, "KdpsTree_4D" );
    exportTree<InMemTypeDef<5>>( m, "KdpsTree_5D" );
}