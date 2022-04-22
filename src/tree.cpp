#include "sps/default.h"
#include "sps/main.h"
#include <stdlib.h>

template <size_t D, bool dependant_dim> using TD = CachedTypeDef<D, dependant_dim>;

PYBIND11_MODULE( libSps, m )
{
    // prevent creation of stxxl log files
    if( getenv( (char*)"STXXLLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLLOGFILE=/dev/null" );
    if( getenv( (char*)"STXXLERRLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLERRLOGFILE=/dev/null" );

    // export various types
    exportSparseCoord<TD<2, false>>( m, "SparseCoords" );

    exportStream<TD<2, false>>( m, "__ProgressOutStream" );

    exportMain<TD<2, false>>( m, "SparsePrefixSum_2D" );
    exportMain<TD<3, false>>( m, "SparsePrefixSum_3D" );
    exportMain<TD<4, false>>( m, "SparsePrefixSum_4D" );
    exportMain<TD<5, false>>( m, "SparsePrefixSum_5D" );

    exportMain<TD<2, true>>( m, "DependantDimSparsePrefixSum_2D" );
    exportMain<TD<3, true>>( m, "DependantDimSparsePrefixSum_3D" );
    exportMain<TD<4, true>>( m, "DependantDimSparsePrefixSum_4D" );
    exportMain<TD<5, true>>( m, "DependantDimSparsePrefixSum_5D" );
}