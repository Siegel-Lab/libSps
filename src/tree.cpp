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
    exportSparseCoord<OnDiskTypeDef<2>>( m, "SparseCoords" );

    exportStream<OnDiskTypeDef<2>>( m, "__ProgressOutStream" );
    exportMain<OnDiskTypeDef<2>>( m, "SparsePrefixSum_2D" );
    exportMain<OnDiskTypeDef<3>>( m, "SparsePrefixSum_3D" );
    exportMain<OnDiskTypeDef<4>>( m, "SparsePrefixSum_4D" );
    exportMain<OnDiskTypeDef<5>>( m, "SparsePrefixSum_5D" );
}