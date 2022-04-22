#include "sps/default.h"
#include "sps/main.h"
#include <stdlib.h>

PYBIND11_MODULE( libSps, m )
{
    // prevent creation of stxxl log files
    if(getenv((char*) "STXXLLOGFILE") == nullptr)
        putenv((char*) "STXXLLOGFILE=/dev/null");
    if(getenv((char*) "STXXLERRLOGFILE") == nullptr)
        putenv((char*) "STXXLERRLOGFILE=/dev/null");

    // export various types
    exportSparseCoord<DiskTypeDef<2, false>>( m, "SparseCoords" );

    exportStream<DiskTypeDef<2, false>>( m, "__ProgressOutStream" );
    
    exportMain<DiskTypeDef<2, false>>( m, "SparsePrefixSum_2D" );
    exportMain<DiskTypeDef<3, false>>( m, "SparsePrefixSum_3D" );
    exportMain<DiskTypeDef<4, false>>( m, "SparsePrefixSum_4D" );
    exportMain<DiskTypeDef<5, false>>( m, "SparsePrefixSum_5D" );

    exportMain<DiskTypeDef<2, true>>( m, "DependantDimSparsePrefixSum_2D" );
    exportMain<DiskTypeDef<3, true>>( m, "DependantDimSparsePrefixSum_3D" );
    exportMain<DiskTypeDef<4, true>>( m, "DependantDimSparsePrefixSum_4D" );
    exportMain<DiskTypeDef<5, true>>( m, "DependantDimSparsePrefixSum_5D" );
}