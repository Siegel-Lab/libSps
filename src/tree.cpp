#include "sps/default.h"
#include "sps/main.h"
#include <stdlib.h>


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template<size_t D, bool dependant_dim>
void exportStorage(pybind11::module& m, std::string sPref, std::string sSuff)
{
#if DISK == ON
    exportMain<DiskTypeDef<D, dependant_dim>>( m, ("Disk" + sPref + "PrefixSum" + sSuff).c_str() );
#endif
#if CACHED == ON
    exportMain<CachedTypeDef<D, dependant_dim>>( m, ("Cached" + sPref + "PrefixSum" + sSuff).c_str() );
#endif
#if RAM == ON
    exportMain<InMemTypeDef<D, dependant_dim>>( m, ("Ram" + sPref + "PrefixSum" + sSuff).c_str() );
#endif
}

template<size_t D>
void exportDependant(pybind11::module& m, std::string sSuff)
{
#if W_DEPENDANT_DIM == ON
    exportStorage<D, true>(m, "DependantDim", sSuff);
#endif
#if WO_DEPENDANT_DIM == ON
    exportStorage<D, false>(m, "", sSuff);
#endif
}

void exportDims(pybind11::module& m)
{
#if NUM_DIMENSIONS_A != 0
    exportDependant<NUM_DIMENSIONS_A>(m, "_" + std::to_string(NUM_DIMENSIONS_A) + "D");
#endif
#if NUM_DIMENSIONS_B != 0
    exportDependant<NUM_DIMENSIONS_B>(m, "_" + std::to_string(NUM_DIMENSIONS_B) + "D");
#endif
#if NUM_DIMENSIONS_C != 0
    exportDependant<NUM_DIMENSIONS_C>(m, "_" + std::to_string(NUM_DIMENSIONS_C) + "D");
#endif
#if NUM_DIMENSIONS_D != 0
    exportDependant<NUM_DIMENSIONS_D>(m, "_" + std::to_string(NUM_DIMENSIONS_D) + "D");
#endif
}
#pragma GCC diagnostic pop

PYBIND11_MODULE( libSps, m )
{
    // prevent creation of stxxl log files
    if( getenv( (char*)"STXXLLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLLOGFILE=/dev/null" );
    if( getenv( (char*)"STXXLERRLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLERRLOGFILE=/dev/null" );

    // export various types
    exportSparseCoord<InMemTypeDef<2, false>>( m, "SparseCoords" );
    exportStream<InMemTypeDef<2, false>>( m, "__ProgressOutStream" );

    exportDims(m);
}