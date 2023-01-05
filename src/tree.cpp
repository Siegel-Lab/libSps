#include "sps/abstract_index.h"

#include "sps/build_time.h"
#include "sps/default.h"
#include "sps/index.h"
#include "sps/simple_vector.h"
#include "sps/version.h"
#include <fstream>
#include <string_view>
#include <stdlib.h>
#include <unistd.h>


constexpr bool bBinSearchSparse = false;

template <size_t D, size_t O, template <size_t, size_t, bool> typename storage_t>
std::string exportPrefixSumIndex( pybind11::module& m, std::string sPref )
{
    if constexpr(D != 0)
    {
        std::string sDesc = "for";
        if( O == 0 )
            sDesc += " points";
        else if( O == 1 )
            sDesc += " intervals";
        else if( O == 2 )
            sDesc += " rectangles";
        else if( O == 3 )
            sDesc += " cubes";
        else
            sDesc += " " + std::to_string( D ) + "-orthotopes";
        sDesc += " in " + std::to_string( D - O ) + "-dimensional space";

        return exportIndex<storage_t<D + O, O, bBinSearchSparse>>(
            m, ( sPref + "PrefixSum_" + std::to_string(D) + "D_" + std::to_string(O) + "O" ).c_str( ), sDesc );
    }
    return "";
}

template <size_t D, size_t O, size_t storage_t>
std::string exportPrefixSumIndex( pybind11::module& m )
{
    if constexpr(storage_t == 0)
        return exportPrefixSumIndex<D, O, InMemTypeDef>(m, "Ram");
    else if constexpr(storage_t == 1)
        return exportPrefixSumIndex<D, O, DiskTypeDef>(m, "Disk");
    else
        return exportPrefixSumIndex<D, O, CachedTypeDef>(m, "Cached");
}

template <size_t D> std::string exportStorageSimpleVec( pybind11::module& m, std::string sSuff )
{
    std::string sRet = "";
#ifdef DISK
    sRet += exportSimpleVector<DiskTypeDef<D, D, true>>( m, ( "DiskSimpleVector" + sSuff ).c_str( ) );
#endif
#ifdef CACHED
    sRet += exportSimpleVector<CachedTypeDef<D, D, true>>( m, ( "CachedSimpleVector" + sSuff ).c_str( ) );
#endif
#ifdef RAM
    sRet += exportSimpleVector<InMemTypeDef<D, D, true>>( m, ( "RamSimpleVector" + sSuff ).c_str( ) );
#endif
    return sRet;
}


template <size_t D, size_t O, template <size_t, size_t, bool> typename storage_t>
std::unique_ptr<AbstractIndex> factoryHelper( size_t uiD, size_t uiOrthtopeDims,
                                              std::string sPrefix, bool bWrite, bool bSimpleVec )
{
    if(uiD == D && uiOrthtopeDims == O)
    {
        if( bSimpleVec )
            return std::make_unique<sps::SimpleVector<storage_t<D, D, true>>>( sPrefix, bWrite );
        else
            return std::make_unique<sps::Index<storage_t<D + O, O, bBinSearchSparse>>>( sPrefix, bWrite );
    }
    return nullptr;
}

template <size_t D, size_t O, size_t storage_t>
std::unique_ptr<AbstractIndex> factoryHelper( [[maybe_unused]] size_t uiD, [[maybe_unused]] size_t uiOrthtopeDims,
                                              [[maybe_unused]] std::string sStorageType,
                                              [[maybe_unused]] std::string sPrefix, [[maybe_unused]] bool bWrite, [[maybe_unused]] bool bSimpleVec )
{
    if constexpr(D != 0)
    {
        if constexpr(storage_t == 0)
            if(sStorageType == "Ram")
                return factoryHelper<D, O, InMemTypeDef>(uiD, uiOrthtopeDims, sPrefix, bWrite, bSimpleVec);
        if constexpr(storage_t == 1)
            if(sStorageType == "Disk")
                return factoryHelper<D, O, DiskTypeDef>(uiD, uiOrthtopeDims, sPrefix, bWrite, bSimpleVec);
        if constexpr(storage_t == 2)
            if(sStorageType == "Cached")
                return factoryHelper<D, O, CachedTypeDef>(uiD, uiOrthtopeDims, sPrefix, bWrite, bSimpleVec);
    }
    return nullptr;
}



std::unique_ptr<AbstractIndex> factory( std::string sPrefix, size_t uiD, size_t uiOrthtopeDims,
                                        std::string sStorageType, bool bWrite )
{
    bool bSimpleVec = false;
    std::unique_ptr<AbstractIndex> pRet = nullptr;
    pRet = factoryHelper<DIMENSIONS_A, ORTHOTOPE_A, STORAGE_A>(uiD, uiOrthtopeDims, sStorageType, sPrefix, bWrite,  
                                                               bSimpleVec);
    if(pRet != nullptr)
        return pRet;
    pRet = factoryHelper<DIMENSIONS_B, ORTHOTOPE_B, STORAGE_B>(uiD, uiOrthtopeDims, sStorageType, sPrefix, bWrite,  
                                                               bSimpleVec);
    if(pRet != nullptr)
        return pRet;
    pRet = factoryHelper<DIMENSIONS_C, ORTHOTOPE_C, STORAGE_C>(uiD, uiOrthtopeDims, sStorageType, sPrefix, bWrite,  
                                                               bSimpleVec);
    if(pRet != nullptr)
        return pRet;
    pRet = factoryHelper<DIMENSIONS_D, ORTHOTOPE_D, STORAGE_D>(uiD, uiOrthtopeDims, sStorageType, sPrefix, bWrite,  
                                                               bSimpleVec);
    if(pRet != nullptr)
        return pRet;
    pRet = factoryHelper<DIMENSIONS_E, ORTHOTOPE_E, STORAGE_E>(uiD, uiOrthtopeDims, sStorageType, sPrefix, bWrite,  
                                                               bSimpleVec);
    if(pRet != nullptr)
        return pRet;
    throw std::invalid_argument( "libSps has not been compiled with the requested parameter combination." );
}
#pragma GCC diagnostic pop

PYBIND11_MODULE( libSps, m )
{
    // prevent creation of stxxl log files
    if( getenv( (char*)"STXXLLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLLOGFILE=/dev/null" );
    if( getenv( (char*)"STXXLERRLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLERRLOGFILE=/dev/null" );

    m.attr( "VERSION" ) = SPS_VERSION;
    m.attr( "BUILD_TIME" ) = SPS_BUILD_TIME;

    exportEnum( m );

    // export various types
    pybind11::class_<sps::AbstractIndex>( m, "AbstractIndex",
                                          R"pbdoc(
    Abstract Index class. Not usefull on it's own.
)pbdoc" );

    std::string sIndices;
    sIndices += exportPrefixSumIndex<DIMENSIONS_A, ORTHOTOPE_A, STORAGE_A>(m);
    sIndices += exportPrefixSumIndex<DIMENSIONS_B, ORTHOTOPE_B, STORAGE_B>(m);
    sIndices += exportPrefixSumIndex<DIMENSIONS_C, ORTHOTOPE_C, STORAGE_C>(m);
    sIndices += exportPrefixSumIndex<DIMENSIONS_D, ORTHOTOPE_D, STORAGE_D>(m);
    sIndices += exportPrefixSumIndex<DIMENSIONS_E, ORTHOTOPE_E, STORAGE_E>(m);

    std::string sRaw = R"pbdoc(
Documentation for the Python Module
-----------------------------------
.. currentmodule:: libSps
.. autosummary::
    :toctree: _generate

    make_sps_index
    IntersectionType
    AbstractIndex
)pbdoc";

    m.doc( ) = sRaw + sIndices;

    // export the factory function
    m.def( "make_sps_index", &factory, pybind11::arg( "filepath_prefix" ) = "", pybind11::arg( "num_dimensions" ) = 2,
           pybind11::arg( "num_orthotope_dimensions" ) = 0,
           pybind11::arg( "storage_type" ) = "Ram", pybind11::arg( "open_in_write_mode" ) = true,
           R"pbdoc(
    Factory for the sparse prefix sum indices.
    
    :param filepath_prefix: Prefix path of the index on the filesystem (multiple files with different endings will be created), defaults to "".
    :type filepath_prefix: str

    :param num_dimensions: Number of dimensions for datapoints in the index, defaults to 2.
    :type num_dimensions: int

    :param num_orthotope_dimensions: Number of orthotope dimensions for the index (i.e. in how many dimensions do the data-hyperrectangles have thickness), defaults to 0.
    :type num_orthotope_dimensions: int

    :param storage_type: The way the datastructure interacts with the filesystem, defaults to Ram.
    :type storage_type: str

    :param open_in_write_mode: Open the index in write mode (if this is set to False no changes can be made to the index), defaults to True.
    :type open_in_write_mode: str

    :return: The sparse prefix sum datastructure.
    :rtype: A child of libSps.AbstractIndex.

    
    Acceptable storage_types are:

    * "Disk" to create indices that load all data from a file on startup and store it back to the file on shutdown. Expect them to consume as much RAM as the filesize.
    * "Cached" to create indices that use a cache to load data from and store data to a file dynamically as needed during runtime. Expect this storage type to be slightly slower than the other two options. For large datasets this storage is necessary, as it allows the RAM usage to be independent of the amount of data stored.
    * "Ram" to create indices that store all information in RAM and never interact with the filesystem.
    * "PickByFileSize" to use the Disk option if the index is opened with open_in_write_mode=False and the index file requires less than half of the available RAM, otherwise Cached is used.
)pbdoc" );
}