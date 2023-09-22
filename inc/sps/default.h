#pragma once

#include "sps/type_defs.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>

#ifdef WITH_STXXL
#include <stxxl/io>
#include <stxxl/sort>
#include <stxxl/vector>
#endif

using namespace sps;

/**
 * @brief Catches all print output.
 *
 * Purpose is to be able to direct printed output into a file, to std::cout or somewhere else.
 * StdOutProgressStream directs the prints to std::cout.
 *
 */
struct StdOutProgressStream
{
    std::chrono::time_point<std::chrono::high_resolution_clock> xLastPrint{ };
    Verbosity xVerb;
    Verbosity xCurr;

    /**
     * @brief Whether it is time to print again.
     *
     * Intended for progress prints, that would otherwise spam the output.
     *
     * @return true If enough time e.g. 1sec has passed since the last call to printAgain.
     * @return false otherwise.
     */
    bool printAgain( )
    {
        auto xCurr = std::chrono::high_resolution_clock::now( );
        std::chrono::duration<double, std::milli> xTime = xCurr - xLastPrint;
        if( xTime.count( ) <= 1000 )
            return false;
        xLastPrint = xCurr;
        return true;
    }

    /**
     * @brief Should the given string be printed based on the current verbosity.
     *
     * @return true if the current verbosity is smaller than the maximal verbosity for this stream.
     * @return false otherwise.
     */
    bool active( ) const
    {
        return xVerb > xCurr;
    }

    /**
     * @brief Print some string-like.
     *
     * @tparam T some printable type.
     * @param sStr message to print.
     * @return StdOutProgressStream& this stream.
     */
    template <typename T> StdOutProgressStream& operator<<( const T& sStr )
    {
        if( active( ) )
            std::cout << sStr << std::flush;
        return *this;
    }

    /**
     * @brief change the current verbosity of this stream.
     *
     * @param xNew new verbosity.
     * @return StdOutProgressStream& this stream.
     */
    StdOutProgressStream& operator<<( const Verbosity& xNew )
    {
        this->xCurr = xNew;
        return *this;
    }

    /**
     * @brief Construct a new Std Out Progress Stream object.
     *
     * @param uiVerb the maximal verbosity for this stream.
     */
    StdOutProgressStream( size_t uiVerb )
        : xLastPrint( std::chrono::high_resolution_clock::now( ) ), xVerb( uiVerb ), xCurr( 0 )
    {}
};

/**
 * @brief Sorter for a std::vector
 *
 * @tparam it_t type of std::vector like stl container iterator
 * @tparam cmp_t type of comparison function object (see std::sort)
 */
template <typename it_t, typename cmp_t> struct RamVectorSorter
{
    /**
     * @brief sorting function.
     *
     * @param xBegin iterator to the first element that shall be sorted.
     * @param xEnd iterator to one past the last element that shall be sorted.
     * @param xComp comparison function object (see std::sort)
     */
    void operator( )( it_t xBegin, it_t xEnd, cmp_t xComp ) const
    {
        std::sort( xBegin, xEnd, xComp );
    }
};

template <typename val_t> struct RamVecGenerator
{
    using file_t = size_t;
    using vec_t = std::vector<val_t>;
    static const bool THREADSAVE = true;

    template <typename it_t, typename cmp_t> using sorter_t = RamVectorSorter<it_t, cmp_t>;

    vec_t vec( file_t& )
    {
        return vec_t( );
    }

    void reserve( size_t uiSize, vec_t& rVec )
    {
        if( uiSize > rVec.capacity( ) )
            rVec.reserve( std::max( uiSize, rVec.capacity( ) + rVec.capacity( ) / 5 ) );
    }

    void try_shrink_to_fit( vec_t& /*rVec*/ )
    {
        // intentionally ignore shrink to fit for avoiding endless reallocations
        // rVec.shrink_to_fit( );
    }

    file_t file( std::string, bool )
    {
        return 0;
    }
};


using default_coordinate_t = uint64_t;
using default_val_t = uint32_t;
using default_class_key_t = uint32_t;
static const size_t DEFAULT_ORTHOTOPE_POWER = 0;
static const bool EXPLAIN = false;

template <size_t D, size_t orthope, size_t orthotope_power>
using PowerInMemTypeDef = TypeDefs<default_coordinate_t, //
                                   default_val_t, //
                                   D, //
                                   default_class_key_t, //
                                   RamVecGenerator, //
                                   RamVecGenerator, //
                                   RamVecGenerator, //
                                   RamVecGenerator, //
                                   RamVecGenerator, //
                                   RamVecGenerator, //
                                   orthope, //
                                   orthotope_power, //
                                   EXPLAIN, //
                                   StdOutProgressStream>;

template <size_t D, size_t orthope, typename val_t>
using ValInMemTypeDef = TypeDefs<default_coordinate_t, //
                                 val_t, //
                                 D, //
                                 default_class_key_t, //
                                 RamVecGenerator, //
                                 RamVecGenerator, //
                                 RamVecGenerator, //
                                 RamVecGenerator, //
                                 RamVecGenerator, //
                                 RamVecGenerator, //
                                 orthope, //
                                 DEFAULT_ORTHOTOPE_POWER, //
                                 EXPLAIN, //
                                 StdOutProgressStream>;
/**
 * @brief Type definitions for a RAM Index
 *
 * Index stores all information in RAM and never interacts with the filesystem.
 *
 * @tparam D number of dimensions
 * @tparam orthope number of orthope dimensions
 */
template <size_t D, size_t orthope> using InMemTypeDef = ValInMemTypeDef<D, orthope, default_val_t>;


#ifdef WITH_STXXL
template <typename val_t, size_t PageSize, size_t CachePages, size_t BlockSize>
struct CacheVec : public stxxl::VECTOR_GENERATOR<val_t, PageSize, CachePages, BlockSize>::result
{
    using vec_t = typename stxxl::VECTOR_GENERATOR<val_t, PageSize, CachePages, BlockSize>::result;

  public:
    using vec_t::vec_t;

    CacheVec( )
    {}
};

template <typename it_t, typename cmp_t> struct CachedVectorSorter
{
    void operator( )( it_t xBegin, it_t xEnd, cmp_t xComp ) const
    {
        stxxl::sort( xBegin, xEnd, xComp, 1024 * 1024 * 1024 );
    }
};

template <typename val_t> struct CachedVecGenerator
{
    using file_t = stxxl::syscall_file;
    using vec_t = CacheVec<val_t, 1, 4096 * 16, 4096>;

    static const bool THREADSAVE = false;

    template <typename it_t, typename cmp_t> using sorter_t = CachedVectorSorter<it_t, cmp_t>;

    vec_t vec( file_t& rFile )
    {
        return vec_t( &rFile );
    }

    void try_shrink_to_fit( vec_t& /*rVec*/ )
    {
        // intentionally ignore shrink to fit for avoiding endless reallocations
        // rVec.resize( rVec.size( ), true );
    }

    void reserve( size_t uiSize, vec_t& rVec )
    {
        rVec.reserve( uiSize );
    }

    file_t file( std::string sPath, bool bOpenInWriteMode )
    {
        return file_t( sPath, ( bOpenInWriteMode ? stxxl::file::open_mode::RDWR | stxxl::file::open_mode::CREAT
                                                 : stxxl::file::open_mode::RDONLY | stxxl::file::open_mode::NO_LOCK ) |
                                  stxxl::file::open_mode::DIRECT );
    }
};
#endif


template <typename val_t> struct DiskVec : public std::vector<val_t>
{
    std::pair<std::string, bool>* fileInfo;

  public:
    DiskVec( std::pair<std::string, bool>* fileInfo ) : std::vector<val_t>( ), fileInfo( fileInfo )
    {
        auto ifstream =
            std::ifstream( fileInfo->first, std::ios_base::in | std::ios_base::out | std::ios_base::binary );

        ifstream.unsetf( std::ios::skipws );
        std::streampos fileSize;

        ifstream.seekg( 0, std::ios::end );
        fileSize = ifstream.tellg( );
        if( fileSize > 0 )
        {
            ifstream.seekg( 0, std::ios::beg );

            assert( fileSize % sizeof( val_t ) == 0 );

            this->resize( fileSize / sizeof( val_t ), val_t( ) );
            ifstream.read( (char*)this->data( ), fileSize );
        }
        ifstream.close( );
    }

    DiskVec( )
    {}

    ~DiskVec( )
    {
        if( fileInfo->second )
        {
            auto ofstream =
                std::ofstream( fileInfo->first, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc );

            ofstream.write( (char*)this->data( ), sizeof( val_t ) * this->size( ) );
            ofstream.flush( );
            ofstream.close( );
        }
    }

    friend std::ostream& operator<<( std::ostream& os, const DiskVec& rIt )
    {
        return os << (std::vector<val_t>)rIt;
    }
};

/**
 * @brief Generator class for a std::vector like
 *
 * @tparam val_t content type for the vector
 */
template <typename val_t> struct DiskVecGenerator
{

    /// @brief type of the file information
    using file_t = std::pair<std::string, bool>;

    /// @brief type of the actual vector
    using vec_t = DiskVec<val_t>;

    /// @brief this vector is threadsave (as long as it's size is not changed)
    static const bool THREADSAVE = false; // true;

    /**
     * @brief defines a sorter object for the generated vector
     */
    template <typename it_t, typename cmp_t> using sorter_t = RamVectorSorter<it_t, cmp_t>;


    /**
     * @brief generate a vector for a given file.
     *
     * @param rFile the file.
     * @return vec_t the generated vector.
     */
    vec_t vec( file_t& rFile )
    {
        return vec_t( &rFile );
    }

    /**
     * @brief reserve memory for the vector
     *
     * @param uiSize keep memory for at least uiSize many elements.
     * @param rVec the vector.
     */
    void reserve( size_t uiSize, vec_t& rVec )
    {
        if( uiSize > rVec.capacity( ) )
            rVec.reserve( std::max( uiSize, rVec.capacity( ) + rVec.capacity( ) / 5 ) );
    }

    /**
     * @brief shrink the capacity of the given vector to fit it's content
     * disabled at the moment.
     */
    void try_shrink_to_fit( vec_t& /*rVec*/ )
    {
        // intentionally ignore shrink to fit for avoiding endless reallocations
        // rVec.shrink_to_fit( );
    }

    /**
     * @brief generate a file for a given path.
     *
     * @param sPath path of the file.
     * @param bOpenInWriteMode whether the file shall be opened in write mode.
     * @return file_t the generated file.
     */
    file_t file( std::string sPath, bool bOpenInWriteMode )
    {
        return std::make_pair( sPath, bOpenInWriteMode );
    }
};

template <size_t D, size_t orthope, size_t orthotope_power>
using PowerDiskTypeDef = TypeDefs<default_coordinate_t, //
                                  default_val_t, //
                                  D, //
                                  default_class_key_t, //
                                  DiskVecGenerator, //
                                  DiskVecGenerator, //
                                  DiskVecGenerator, //
                                  DiskVecGenerator, //
                                  DiskVecGenerator, //
                                  DiskVecGenerator, //
                                  orthope, //
                                  orthotope_power, //
                                  EXPLAIN, //
                                  StdOutProgressStream>;

template <size_t D, size_t orthope, typename val_t>
using ValDiskTypeDef = TypeDefs<default_coordinate_t, //
                                val_t, //
                                D, //
                                default_class_key_t, //
                                DiskVecGenerator, //
                                DiskVecGenerator, //
                                DiskVecGenerator, //
                                DiskVecGenerator, //
                                DiskVecGenerator, //
                                DiskVecGenerator, //
                                orthope, //
                                DEFAULT_ORTHOTOPE_POWER, //
                                EXPLAIN, //
                                StdOutProgressStream>;
/**
 * @brief Type definitions for a Disk Index
 *
 * Index that loads all data from a file on startup and store it back to the file on shutdown.
 * Expect it to consume as much RAM as the filesize.
 *
 * @tparam D number of dimensions
 * @tparam orthope number of orthope dimensions
 */
template <size_t D, size_t orthope> using DiskTypeDef = ValDiskTypeDef<D, orthope, default_val_t>;


#ifdef WITH_STXXL

template <size_t D, size_t orthope, size_t orthotope_power>
using PowerCachedTypeDef = TypeDefs<default_coordinate_t, //
                                    default_val_t, //
                                    D, //
                                    default_class_key_t, //
                                    DiskVecGenerator, //
                                    DiskVecGenerator, //
                                    DiskVecGenerator, //
                                    DiskVecGenerator, //
                                    CachedVecGenerator, //
                                    CachedVecGenerator, //
                                    orthope, //
                                    orthotope_power, //
                                    EXPLAIN, //
                                    StdOutProgressStream>;

template <size_t D, size_t orthope, typename val_t>
using ValCachedTypeDef = TypeDefs<default_coordinate_t, //
                                  val_t, //
                                  D, //
                                  default_class_key_t, //
                                  DiskVecGenerator, //
                                  DiskVecGenerator, //
                                  DiskVecGenerator, //
                                  DiskVecGenerator, //
                                  CachedVecGenerator, //
                                  CachedVecGenerator, //
                                  orthope, //
                                  DEFAULT_ORTHOTOPE_POWER, //
                                  EXPLAIN, //
                                  StdOutProgressStream>;
/**
 * @brief Type definitions for a Cached Index
 *
 * Index that uses a cache to load data from and store data to a file dynamically as needed during runtime.
 * Expect this storage type to be slightly slower than the other two options.
 * For large datasets this storage is necessary, as it allows the RAM usage to be independent of the amount of data
 * stored.
 *
 * @tparam D number of dimensions
 * @tparam orthope number of orthope dimensions
 */
template <size_t D, size_t orthope> using CachedTypeDef = ValCachedTypeDef<D, orthope, default_val_t>;

#endif