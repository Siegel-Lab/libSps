#include "sps/abstract_index.h"

#include "sps/default.h"
#include "sps/index.h"
#include "sps/simple_vector.h"
#include "sps/version.h"
#include <fstream>
#include <stdlib.h>
#include <unistd.h>

#define STRINGIFY( s ) XSTRINGIFY( s )
#define XSTRINGIFY( s ) #s

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

template <size_t D> std::string exportStorageSimpleVec( pybind11::module& m, std::string sSuff )
{
    std::string sRet = "";
#ifdef DISK
    sRet += exportSimpleVector<DiskTypeDef<D, false, D>>( m, ( "DiskSimpleVector" + sSuff ).c_str( ) );
#endif
#ifdef CACHED
    sRet += exportSimpleVector<CachedTypeDef<D, false, D>>( m, ( "CachedSimpleVector" + sSuff ).c_str( ) );
#endif
#ifdef RAM
    sRet += exportSimpleVector<InMemTypeDef<D, false, D>>( m, ( "RamSimpleVector" + sSuff ).c_str( ) );
#endif
    return sRet;
}

template <size_t D, bool dependant_dim, size_t orthope>
std::string exportStorage( pybind11::module& m, std::string sPref, std::string sSuff )
{
    std::string sDesc = "for";
    if( orthope == 0 )
        sDesc += " points";
    else if( orthope == 1 )
        sDesc += " intervals";
    else if( orthope == 2 )
        sDesc += " rectangles";
    else if( orthope == 3 )
        sDesc += " cubes";
    else
        sDesc += " " + std::to_string( D ) + "-orthotopes";
    sDesc += " in " + std::to_string( D - orthope ) + "-dimensional space";
    if( dependant_dim )
        sDesc += ", with a dependent dimension";
    std::string sRet = "";
#ifdef DISK
    sRet += exportIndex<DiskTypeDef<D, dependant_dim, orthope>>( m, ( "Disk" + sPref + "PrefixSum" + sSuff ).c_str( ),
                                                                 sDesc );
#endif
#ifdef CACHED
    sRet += exportIndex<CachedTypeDef<D, dependant_dim, orthope>>(
        m, ( "Cached" + sPref + "PrefixSum" + sSuff ).c_str( ), sDesc );
#endif
#ifdef RAM
    sRet += exportIndex<InMemTypeDef<D, dependant_dim, orthope>>( m, ( "Ram" + sPref + "PrefixSum" + sSuff ).c_str( ),
                                                                  sDesc );
#endif
    return sRet;
}

template <size_t D, bool dependant_dim>
std::string exportOrthope( pybind11::module& m, std::string sPref, std::string sSuff )
{
    std::string sRet = "";
#ifdef W_CUBES
    sRet += exportStorage<D + 3, dependant_dim, 3>( m, sPref + "Cubes", sSuff );
#endif
#ifdef W_RECTANGLES
    sRet += exportStorage<D + 2, dependant_dim, 2>( m, sPref + "Rectangles", sSuff );
#endif
#ifdef W_INTERVALS
    sRet += exportStorage<D + 1, dependant_dim, 1>( m, sPref + "Intervals", sSuff );
#endif
#ifdef W_POINTS
    sRet += exportStorage<D, dependant_dim, 0>( m, sPref + "Points", sSuff );
#endif
    return sRet;
}

template <size_t D> std::string exportDependant( pybind11::module& m, std::string sSuff )
{
    std::string sRet = "";
#ifdef W_DEPENDANT_DIM
    sRet += exportOrthope<D, true>( m, "DependantDim", sSuff );
#endif
#ifdef WO_DEPENDANT_DIM
    sRet += exportOrthope<D, false>( m, "", sSuff );
#endif
    return sRet;
}

std::string exportDims( pybind11::module& m )
{
    std::string sRet = "";
#if NUM_DIMENSIONS_A != 0
    sRet += exportDependant<NUM_DIMENSIONS_A>( m, "_" + std::to_string( NUM_DIMENSIONS_A ) + "D" );
    sRet += exportStorageSimpleVec<NUM_DIMENSIONS_A>( m, "_" + std::to_string( NUM_DIMENSIONS_A ) + "D" );
#endif
#if NUM_DIMENSIONS_B != 0
    sRet += exportDependant<NUM_DIMENSIONS_B>( m, "_" + std::to_string( NUM_DIMENSIONS_B ) + "D" );
    sRet += exportStorageSimpleVec<NUM_DIMENSIONS_B>( m, "_" + std::to_string( NUM_DIMENSIONS_B ) + "D" );
#endif
#if NUM_DIMENSIONS_C != 0
    sRet += exportDependant<NUM_DIMENSIONS_C>( m, "_" + std::to_string( NUM_DIMENSIONS_C ) + "D" );
    sRet += exportStorageSimpleVec<NUM_DIMENSIONS_C>( m, "_" + std::to_string( NUM_DIMENSIONS_C ) + "D" );
#endif
#if NUM_DIMENSIONS_D != 0
    sRet += exportDependant<NUM_DIMENSIONS_D>( m, "_" + std::to_string( NUM_DIMENSIONS_D ) + "D" );
    sRet += exportStorageSimpleVec<NUM_DIMENSIONS_D>( m, "_" + std::to_string( NUM_DIMENSIONS_D ) + "D" );
#endif
    return sRet;
}

std::ifstream::pos_type filesize( const char* filename )
{
    std::ifstream in( filename, std::ifstream::ate | std::ifstream::binary );
    return in.tellg( );
}

size_t getTotalSystemMemory( )
{
    size_t pages = sysconf( _SC_PHYS_PAGES );
    size_t page_size = sysconf( _SC_PAGE_SIZE );
    return pages * page_size;
}

template <size_t D, bool dependant_dim, size_t orthope, template <size_t, bool, size_t> typename storage_t>
std::unique_ptr<AbstractIndex> factoryHelperFinal( std::string sPrefix, bool bWrite, bool bSimpleVec )
{
    if( bSimpleVec )
        return std::make_unique<sps::SimpleVector<storage_t<D, false, D>>>( sPrefix, bWrite );
    else
        return std::make_unique<sps::Index<storage_t<D, dependant_dim, orthope>>>( sPrefix, bWrite );
}

template <size_t D, bool dependant_dim, size_t orthope>
std::unique_ptr<AbstractIndex> factoryHelper( std::string sStorageType, std::string sPrefix, bool bWrite,
                                              bool bSimpleVec )
{
#ifdef DISK
#ifdef CACHED
    if( sStorageType == "PickByFileSize" )
    {
        if( bWrite ) // cached
            return factoryHelperFinal<D, dependant_dim, orthope, CachedTypeDef>( sPrefix, bWrite, bSimpleVec );

        std::ifstream::pos_type uiTotalSize = 0;
        uiTotalSize += filesize( ( sPrefix + ".coords" ).c_str( ) );
        uiTotalSize += filesize( ( sPrefix + ".datasets" ).c_str( ) );
        uiTotalSize += filesize( ( sPrefix + ".desc" ).c_str( ) );
        uiTotalSize += filesize( ( sPrefix + ".meta" ).c_str( ) );
        uiTotalSize += filesize( ( sPrefix + ".overlays" ).c_str( ) );
        uiTotalSize += filesize( ( sPrefix + ".points" ).c_str( ) );
        uiTotalSize += filesize( ( sPrefix + ".prefix_sums" ).c_str( ) );

        if( (size_t)uiTotalSize * 2 < getTotalSystemMemory( ) ) // Disk
            return factoryHelperFinal<D, dependant_dim, orthope, DiskTypeDef>( sPrefix, bWrite, bSimpleVec );
        else // CACHED
            return factoryHelperFinal<D, dependant_dim, orthope, CachedTypeDef>( sPrefix, bWrite, bSimpleVec );
    }
#endif
#endif

#ifdef DISK
    if( sStorageType == "Disk" )
        return factoryHelperFinal<D, dependant_dim, orthope, DiskTypeDef>( sPrefix, bWrite, bSimpleVec );
#endif
#ifdef CACHED
    if( sStorageType == "Cached" )
        return factoryHelperFinal<D, dependant_dim, orthope, CachedTypeDef>( sPrefix, bWrite, bSimpleVec );
#endif
#ifdef RAM
    if( sStorageType == "Ram" )
        return factoryHelperFinal<D, dependant_dim, orthope, InMemTypeDef>( sPrefix, bWrite, bSimpleVec );
#endif
    throw std::invalid_argument( "libSps has not been compiled with the requested storage type." );
}

template <size_t D, bool dependant_dim>
std::unique_ptr<AbstractIndex> factoryHelper( size_t uiOrthtopeDims, std::string sStorageType, std::string sPrefix,
                                              bool bWrite, bool bSimpleVec )
{
#ifdef W_CUBES
    if( uiOrthtopeDims == 3 )
        return factoryHelper<D + 3, dependant_dim, 3>( sStorageType, sPrefix, bWrite, bSimpleVec );
#endif
#ifdef W_RECTANGLES
    if( uiOrthtopeDims == 2 )
        return factoryHelper<D + 2, dependant_dim, 2>( sStorageType, sPrefix, bWrite, bSimpleVec );
#endif
#ifdef W_INTERVALS
    if( uiOrthtopeDims == 1 )
        return factoryHelper<D + 1, dependant_dim, 1>( sStorageType, sPrefix, bWrite, bSimpleVec );
#endif
#ifdef W_POINTS
    if( uiOrthtopeDims == 0 )
        return factoryHelper<D, dependant_dim, 0>( sStorageType, sPrefix, bWrite, bSimpleVec );
#endif
    throw std::invalid_argument( "libSps has not been compiled with the requested number of orthotope dimensions." );
}

template <size_t D>
std::unique_ptr<AbstractIndex> factoryHelper( bool bDependentDimension, size_t uiOrthtopeDims, std::string sStorageType,
                                              std::string sPrefix, bool bWrite, bool bSimpleVec )
{
#ifdef W_DEPENDANT_DIM
    if( bDependentDimension )
        return factoryHelper<D, true>( uiOrthtopeDims, sStorageType, sPrefix, bWrite, bSimpleVec );
#endif
#ifdef WO_DEPENDANT_DIM
    if( !bDependentDimension )
        return factoryHelper<D, false>( uiOrthtopeDims, sStorageType, sPrefix, bWrite, bSimpleVec );
#endif
    throw std::invalid_argument( "libSps has not been compiled with the requested dependent dimension configuration." );
}

std::unique_ptr<AbstractIndex> factory( std::string sPrefix, size_t uiD, bool bDependentDimension,
                                        size_t uiOrthtopeDims, std::string sStorageType, bool bWrite )
{
    bool bSimpleVec = false;
#if NUM_DIMENSIONS_A != 0
    if( NUM_DIMENSIONS_A == uiD )
        return factoryHelper<NUM_DIMENSIONS_A>( bDependentDimension, uiOrthtopeDims, sStorageType, sPrefix, bWrite,
                                                bSimpleVec );
#endif
#if NUM_DIMENSIONS_B != 0
    if( NUM_DIMENSIONS_B == uiD )
        return factoryHelper<NUM_DIMENSIONS_B>( bDependentDimension, uiOrthtopeDims, sStorageType, sPrefix, bWrite,
                                                bSimpleVec );
#endif
#if NUM_DIMENSIONS_C != 0
    if( NUM_DIMENSIONS_C == uiD )
        return factoryHelper<NUM_DIMENSIONS_C>( bDependentDimension, uiOrthtopeDims, sStorageType, sPrefix, bWrite,
                                                bSimpleVec );
#endif
#if NUM_DIMENSIONS_D != 0
    if( NUM_DIMENSIONS_D == uiD )
        return factoryHelper<NUM_DIMENSIONS_D>( bDependentDimension, uiOrthtopeDims, sStorageType, sPrefix, bWrite,
                                                bSimpleVec );
#endif
    throw std::invalid_argument( "libSps has not been compiled with the requested number of dimensions." );
}
#pragma GCC diagnostic pop

PYBIND11_MODULE( libSps, m )
{
    // prevent creation of stxxl log files
    if( getenv( (char*)"STXXLLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLLOGFILE=/dev/null" );
    if( getenv( (char*)"STXXLERRLOGFILE" ) == nullptr )
        putenv( (char*)"STXXLERRLOGFILE=/dev/null" );

    m.attr( "VERSION" ) = VERSION;

    // export various types
    pybind11::class_<sps::AbstractIndex>( m, "AbstractIndex",
                                          R"pbdoc(
    Abstract Index class. Not usefull on it's own.
)pbdoc" );

    std::string sIndices = exportDims( m );
    std::string sRaw = R"pbdoc(
Documentation for the Python Module
-----------------------------------
.. currentmodule:: libSps
.. autosummary::
    :toctree: _generate

    make_sps_index
    AbstractIndex
)pbdoc";

    m.doc( ) = sRaw + sIndices;

    // export the factory function
    m.def( "make_sps_index", &factory, pybind11::arg( "filepath_prefix" ) = "", pybind11::arg( "num_dimensions" ) = 2,
           pybind11::arg( "with_dependent_dimension" ) = false, pybind11::arg( "num_orthotope_dimensions" ) = 0,
           pybind11::arg( "storage_type" ) = "Ram", pybind11::arg( "open_in_write_mode" ) = true,
           R"pbdoc(
    Factory for the sparse prefix sum indices.
    
    :param filepath_prefix: Prefix path of the index on the filesystem (multiple files with different endings will be created), defaults to "".
    :type filepath_prefix: str

    :param num_dimensions: Number of dimensions for datapoints in the index, defaults to 2.
    :type num_dimensions: int

    :param with_dependent_dimension: Whether the overlay grid in dimension 1 is dependent on the grid in dimension 0, defaults to False.
    :type with_dependent_dimension: bool

    :param num_orthotope_dimensions: Number of orthotope dimensions (set this to 1, 2, 3, ... to add intervals, rectangles, cubes, ... instead of points), defaults to 0.
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