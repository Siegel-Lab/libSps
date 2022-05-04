#include "sps/default.h"
#include "sps/main.h"
#include <stdlib.h>
#include <fstream>
#include <unistd.h>

#define STRINGIFY(s) XSTRINGIFY(s)
#define XSTRINGIFY(s) #s

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template<size_t D, bool dependant_dim, size_t orthope>
void exportStorage(pybind11::module& m, std::string sPref, std::string sSuff)
{
#ifdef DISK
    exportMain<DiskTypeDef<D, dependant_dim, orthope>>( m, ("Disk" + sPref + "PrefixSum" + sSuff).c_str() );
#endif
#ifdef CACHED
    exportMain<CachedTypeDef<D, dependant_dim, orthope>>( m, ("Cached" + sPref + "PrefixSum" + sSuff).c_str() );
#endif
#ifdef RAM
    exportMain<InMemTypeDef<D, dependant_dim, orthope>>( m, ("Ram" + sPref + "PrefixSum" + sSuff).c_str() );
#endif
}

template<size_t D, bool dependant_dim>
void exportOrthope(pybind11::module& m, std::string sPref, std::string sSuff)
{
#ifdef W_CUBES
    exportStorage<D + 3, dependant_dim, 3>(m, sPref + "Cubes", sSuff);
#endif
#ifdef W_RECTANGLES
    exportStorage<D + 2, dependant_dim, 2>(m, sPref + "Rectangles", sSuff);
#endif
#ifdef W_INTERVALS
    exportStorage<D + 1, dependant_dim, 1>(m, sPref + "Intervals", sSuff);
#endif
#ifdef W_POINTS
    exportStorage<D, dependant_dim, 0>(m, sPref + "Points", sSuff);
#endif
}

template<size_t D>
void exportDependant(pybind11::module& m, std::string sSuff)
{
#ifdef W_DEPENDANT_DIM
    exportOrthope<D, true>(m, "DependantDim", sSuff);
#endif
#ifdef WO_DEPENDANT_DIM
    exportOrthope<D, false>(m, "", sSuff);
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

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

size_t getTotalSystemMemory()
{
    size_t pages = sysconf(_SC_PHYS_PAGES);
    size_t page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}

template<size_t D, bool dependant_dim, size_t orthope>
std::unique_ptr<AbstractMain> factoryHelper(std::string sStorageType, std::string sPrefix, bool bWrite )
{
#ifdef DISK 
#ifdef CACHED
    if(sStorageType == "PickByFileSize")
    {
        if(bWrite) // cached
            return std::make_unique<sps::Main<CachedTypeDef<D, dependant_dim, orthope>>>(sPrefix, bWrite);
        
        std::ifstream::pos_type uiTotalSize = 0;
        uiTotalSize += filesize((sPrefix + ".coords").c_str());
        uiTotalSize += filesize((sPrefix + ".datasets").c_str());
        uiTotalSize += filesize((sPrefix + ".desc").c_str());
        uiTotalSize += filesize((sPrefix + ".meta").c_str());
        uiTotalSize += filesize((sPrefix + ".overlays").c_str());
        uiTotalSize += filesize((sPrefix + ".points").c_str());
        uiTotalSize += filesize((sPrefix + ".prefix_sums").c_str());

        if( (size_t)uiTotalSize * 2 < getTotalSystemMemory()) // RAM
            return std::make_unique<sps::Main<InMemTypeDef<D, dependant_dim, orthope>>>(sPrefix, bWrite);
        else // CACHED
            return std::make_unique<sps::Main<CachedTypeDef<D, dependant_dim, orthope>>>(sPrefix, bWrite);
    }
#endif
#endif

#ifdef DISK
    if(sStorageType == "Disk")
        return std::make_unique<sps::Main<DiskTypeDef<D, dependant_dim, orthope>>>(sPrefix, bWrite);
#endif
#ifdef CACHED
    if(sStorageType == "Cached")
        return std::make_unique<sps::Main<CachedTypeDef<D, dependant_dim, orthope>>>(sPrefix, bWrite);
#endif
#ifdef RAM
    if(sStorageType == "Ram")
        return std::make_unique<sps::Main<InMemTypeDef<D, dependant_dim, orthope>>>(sPrefix, bWrite);
#endif
    throw std::invalid_argument("libSps has not been compiled with the requested storage type.");
}

template<size_t D, bool dependant_dim>
std::unique_ptr<AbstractMain> factoryHelper(size_t uiOrthtopeDims, std::string sStorageType, std::string sPrefix, bool bWrite )
{
#ifdef W_CUBES
    if(uiOrthtopeDims == 3)
        return factoryHelper<D + 3, true, 3>(sStorageType, sPrefix, bWrite);
#endif
#ifdef W_RECTANGLES
    if(uiOrthtopeDims == 2)
        return factoryHelper<D + 2, true, 2>(sStorageType, sPrefix, bWrite);
#endif
#ifdef W_INTERVALS
    if(uiOrthtopeDims == 1)
        return factoryHelper<D + 1, true, 1>(sStorageType, sPrefix, bWrite);
#endif
#ifdef W_POINTS
    if(uiOrthtopeDims == 0)
        return factoryHelper<D, true, 0>(sStorageType, sPrefix, bWrite);
#endif
    throw std::invalid_argument("libSps has not been compiled with the requested number of orthotope dimensions.");
}

template<size_t D>
std::unique_ptr<AbstractMain> factoryHelper(bool bDependentDimension, size_t uiOrthtopeDims, 
                                      std::string sStorageType, std::string sPrefix, bool bWrite )
{
    #ifdef W_DEPENDANT_DIM
    if(bDependentDimension)
        return factoryHelper<D, true>(uiOrthtopeDims, sStorageType, sPrefix, bWrite);
#endif
#ifdef WO_DEPENDANT_DIM
    if(!bDependentDimension)
        return factoryHelper<D, false>(uiOrthtopeDims, sStorageType, sPrefix, bWrite);
#endif
    throw std::invalid_argument("libSps has not been compiled with the requested dependent dimension configuration.");
}

std::unique_ptr<AbstractMain> factory(std::string sPrefix, size_t uiD, bool bDependentDimension, size_t uiOrthtopeDims, 
                                      std::string sStorageType, bool bWrite )
{
#if NUM_DIMENSIONS_A != 0
    if (NUM_DIMENSIONS_A == uiD)
        return factoryHelper<NUM_DIMENSIONS_A>(bDependentDimension, uiOrthtopeDims, sStorageType, sPrefix, bWrite);
#endif
#if NUM_DIMENSIONS_B != 0
    if (NUM_DIMENSIONS_B == uiD)
        return factoryHelper<NUM_DIMENSIONS_B>(bDependentDimension, uiOrthtopeDims, sStorageType, sPrefix, bWrite);
#endif
#if NUM_DIMENSIONS_C != 0
    if (NUM_DIMENSIONS_C == uiD)
        return factoryHelper<NUM_DIMENSIONS_C>(bDependentDimension, uiOrthtopeDims, sStorageType, sPrefix, bWrite);
#endif
#if NUM_DIMENSIONS_D != 0
    if (NUM_DIMENSIONS_D == uiD)
        return factoryHelper<NUM_DIMENSIONS_D>(bDependentDimension, uiOrthtopeDims, sStorageType, sPrefix, bWrite);
#endif
    throw std::invalid_argument("libSps has not been compiled with the requested number of dimensions.");
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
    pybind11::class_<sps::AbstractMain>( m, "__AbstractMain" );
    exportDims(m);

    // export the factory function
    m.def("make_sps_index", &factory, 
            pybind11::arg( "filepath_prefix" ),
            pybind11::arg( "num_dimensions" ),
            pybind11::arg( "with_dependent_dimension" ),
            pybind11::arg( "num_orthotope_dimensions" ),
            pybind11::arg( "storage_type" ),
            pybind11::arg( "open_in_write_mode" ),
          "A factory for the sparse prefix sum datastructure");
}