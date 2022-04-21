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
    exportSparseCoord<OnDiskTypeDef<2, false>>( m, "SparseCoords" );

    exportStream<OnDiskTypeDef<2, false>>( m, "__ProgressOutStream" );
    
    exportMain<OnDiskTypeDef<2, false>>( m, "SparsePrefixSum_2D" );
    exportMain<OnDiskTypeDef<3, false>>( m, "SparsePrefixSum_3D" );
    exportMain<OnDiskTypeDef<4, false>>( m, "SparsePrefixSum_4D" );
    exportMain<OnDiskTypeDef<5, false>>( m, "SparsePrefixSum_5D" );

    exportMain<OnDiskTypeDef<2, true>>( m, "DependantDimSparsePrefixSum_2D" );
    exportMain<OnDiskTypeDef<3, true>>( m, "DependantDimSparsePrefixSum_3D" );
    exportMain<OnDiskTypeDef<4, true>>( m, "DependantDimSparsePrefixSum_4D" );
    exportMain<OnDiskTypeDef<5, true>>( m, "DependantDimSparsePrefixSum_5D" );
}