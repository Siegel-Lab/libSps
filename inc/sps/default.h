#pragma once

#include "sps/type_defs.h"
#include <chrono>
#include <iostream>
#include <stxxl/io>
#include <stxxl/sort>
#include <stxxl/vector>


using namespace sps;

struct StdOutProgressStream
{
    std::chrono::time_point<std::chrono::high_resolution_clock> xLastPrint{ };
    Verbosity xVerb;
    Verbosity xCurr;

    bool printAgain( )
    {
        auto xCurr = std::chrono::high_resolution_clock::now( );
        std::chrono::duration<double, std::milli> xTime = xCurr - xLastPrint;
        if( xTime.count( ) <= 1000 )
            return false;
        xLastPrint = xCurr;
        return true;
    }

    bool active( ) const
    {
        return xVerb > xCurr;
    }

    template <typename T> StdOutProgressStream& operator<<( const T& sStr )
    {
        if( active( ) )
            std::cout << sStr << std::flush;
        return *this;
    }

    StdOutProgressStream& operator<<( const Verbosity& xCurr )
    {
        this->xCurr = xCurr;
        return *this;
    }

    StdOutProgressStream( size_t uiVerb )
        : xLastPrint( std::chrono::high_resolution_clock::now( ) ), xVerb( uiVerb ), xCurr( 0 )
    {}
};


template <typename val_t> struct RamVecGenerator
{
    using file_t = size_t;
    using vec_t = std::vector<val_t>;
    static const bool THREADSAVE = true;

    vec_t vec( file_t& )
    {
        return vec_t( );
    }

    file_t file( std::string, bool )
    {
        return 0;
    }
};


template <typename it_t, typename cmp_t> struct RamVectorSorter
{
    void operator( )( it_t xBegin, it_t xEnd, cmp_t xComp ) const
    {
        std::sort( xBegin, xEnd, xComp );
    }
};

using default_coordinate_t = uint32_t;
using default_val_t = uint32_t;
using default_class_key_t = uint16_t;
static const bool EXPLAIN = false;

template <size_t D, bool dependant_dim, size_t orthope>
using InMemTypeDef = TypeDefs<default_coordinate_t, //
                              default_val_t, //
                              D, //
                              default_class_key_t, //
                              RamVecGenerator, //
                              RamVectorSorter, //
                              dependant_dim, //
                              orthope, //
                              EXPLAIN, //
                              StdOutProgressStream>;


template <typename val_t> struct CachedVecGenerator
{
    using file_t = stxxl::syscall_file;

    using vec_t = typename stxxl::VECTOR_GENERATOR<val_t, 1, 4096 * 16, 4096>::result;

    static const bool THREADSAVE = false;

    vec_t vec( file_t& rFile )
    {
        return vec_t( &rFile );
    }

    file_t file( std::string sPath, bool bOpenInWriteMode )
    {
        return file_t( sPath, ( bOpenInWriteMode ? stxxl::file::open_mode::RDWR | stxxl::file::open_mode::CREAT
                                                 : stxxl::file::open_mode::RDONLY | stxxl::file::open_mode::NO_LOCK ) |
                                  stxxl::file::open_mode::DIRECT );
    }
};

template <typename it_t, typename cmp_t> struct CachedVectorSorter
{
    void operator( )( it_t xBegin, it_t xEnd, cmp_t xComp ) const
    {
        stxxl::sort( xBegin, xEnd, xComp, 1024 * 1024 * 1024 );
    }
};

template <size_t D, bool dependant_dim, size_t orthope>
using CachedTypeDef = TypeDefs<default_coordinate_t, //
                               default_val_t, //
                               D, //
                               default_class_key_t, //
                               CachedVecGenerator, //
                               CachedVectorSorter, //
                               dependant_dim, //
                               orthope, //
                               EXPLAIN, //
                               StdOutProgressStream>;


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

template <typename val_t> struct DiskVecGenerator
{

    using file_t = std::pair<std::string, bool>;
    using vec_t = DiskVec<val_t>;
    static const bool THREADSAVE = true;


    vec_t vec( file_t& rFile )
    {
        return vec_t( &rFile );
    }

    file_t file( std::string sPath, bool bOpenInWriteMode )
    {
        return std::make_pair( sPath, bOpenInWriteMode );
    }
};

template <size_t D, bool dependant_dim, size_t orthope>
using DiskTypeDef = TypeDefs<default_coordinate_t, //
                             default_val_t, //
                             D, //
                             default_class_key_t, //
                             DiskVecGenerator, //
                             RamVectorSorter, //
                             dependant_dim, //
                             orthope, //
                             EXPLAIN, //
                             StdOutProgressStream>;