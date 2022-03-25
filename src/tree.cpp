#include "sps/default.h"
#include "sps/main.h"


PYBIND11_MODULE( libSps, m )
{

    exportStream<OnDiskTypeDef<2>>( m, "__ProgressOutStream" );
    exportMain<OnDiskTypeDef<2>>( m, "SparsePrefixSum_2D" );
}