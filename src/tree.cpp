#include "cstree/default.h"

PYBIND11_MODULE( libcstree, m )
{
    exportTree<InMemTypeDef<1>>( m, "CsTree_1D" );
    exportTree<InMemTypeDef<2>>( m, "CsTree_2D" );
    exportTree<InMemTypeDef<3>>( m, "CsTree_3D" );
    exportTree<InMemTypeDef<4>>( m, "CsTree_4D" );
    exportTree<InMemTypeDef<5>>( m, "CsTree_5D" );
}