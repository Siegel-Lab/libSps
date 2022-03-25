#include "kdpstree/default.h"
#include "kdpstree/main.h"


PYBIND11_MODULE( libKdpsTree, m )
{

    exportStream<OnDiskTypeDef<2>>( m, "__ProgressOutStream" );
    exportMain<OnDiskTypeDef<2>>( m, "SparsePrefixSum_2D" );
}