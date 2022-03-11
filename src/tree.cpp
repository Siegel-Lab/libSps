#include "kdpstree/default.h"

PYBIND11_MODULE( libKdpsTree, m )
{
    exportTree<OnDiskTypeDef<1>>( m, "KdpsTree1" );
    exportTree<OnDiskTypeDef<2>>( m, "KdpsTree2" );
    exportTree<OnDiskTypeDef<3>>( m, "KdpsTree3" );
    exportTree<OnDiskTypeDef<4>>( m, "KdpsTree4" );
    exportTree<OnDiskTypeDef<5>>( m, "KdpsTree5" );
    exportTree<OnDiskTypeDef<6>>( m, "KdpsTree6" );
    exportTree<OnDiskTypeDef<7>>( m, "KdpsTree7" );
    exportTree<OnDiskTypeDef<8>>( m, "KdpsTree8" );
    exportTree<OnDiskTypeDef<9>>( m, "KdpsTree9" );
    exportTree<OnDiskTypeDef<10>>( m, "KdpsTree10" );

    exportArray<OnDiskTypeDef<1>>( m, "PsArray1" );
    exportArray<OnDiskTypeDef<2>>( m, "PsArray2" );
    exportArray<OnDiskTypeDef<3>>( m, "PsArray3" );
    exportArray<OnDiskTypeDef<4>>( m, "PsArray4" );
    exportArray<OnDiskTypeDef<5>>( m, "PsArray5" );
    exportArray<OnDiskTypeDef<6>>( m, "PsArray6" );
    exportArray<OnDiskTypeDef<7>>( m, "PsArray7" );
    exportArray<OnDiskTypeDef<8>>( m, "PsArray8" );
    exportArray<OnDiskTypeDef<9>>( m, "PsArray9" );
    exportArray<OnDiskTypeDef<10>>( m, "PsArray10" );
}