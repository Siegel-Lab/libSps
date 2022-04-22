#pragma once

#include <stxxl/io>
#include <stxxl/sort>
#include <stxxl/vector>
#include <chrono>
#include <iostream>
#include "sps/type_defs.h"


using namespace sps;

struct StdOutProgressStream
{
    std::chrono::time_point<std::chrono::high_resolution_clock> xLastPrint{};
    Verbosity xVerb;
    Verbosity xCurr;

    bool printAgain()
    {
        auto xCurr = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> xTime = xCurr - xLastPrint;
        if(xTime.count() <= 1000)
            return false;
        xLastPrint = xCurr;
        return true;
    }

    bool active() const
    {
        return xVerb > xCurr;
    }

    template<typename T>
    StdOutProgressStream& operator<<(const T& sStr)
    {
        if(active())
            std::cout << sStr << std::flush;
        return *this;
    }
    
    StdOutProgressStream& operator<<(const Verbosity& xCurr)
    {
        this->xCurr = xCurr;
        return *this;
    }

    StdOutProgressStream(size_t uiVerb) :
        xLastPrint(std::chrono::high_resolution_clock::now()), xVerb(uiVerb), xCurr(0)
    {}
};


template <typename val_t, size_t ele_per_block> struct RamVecGenerator
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

template <size_t D, bool dependant_dim>
using InMemTypeDef = TypeDefs<default_coordinate_t, //
                              default_val_t, //
                              D, //
                              default_class_key_t, //
                              RamVecGenerator, //
                              RamVectorSorter, //
                              dependant_dim, //
                              EXPLAIN, //
                              StdOutProgressStream
                              >;
                            


template <typename val_t, size_t ele_per_block> struct CachedVecGenerator
{
    using file_t = stxxl::syscall_file;

    using vec_t = typename stxxl::VECTOR_GENERATOR<val_t, 1, 1024 * 16, 4096>::result;

    static const bool THREADSAVE = false;

    vec_t vec( file_t& rFile )
    {
        return vec_t( &rFile );
    }

    file_t file( std::string sPath, bool bOpenInWriteMode )
    {
        return file_t( sPath, (bOpenInWriteMode ?
                    stxxl::file::open_mode::RDWR | stxxl::file::open_mode::CREAT : 
                    stxxl::file::open_mode::RDONLY | stxxl::file::open_mode::NO_LOCK) 
                | stxxl::file::open_mode::DIRECT );
    }
};

template <typename it_t, typename cmp_t> struct CachedVectorSorter
{
    void operator( )( it_t xBegin, it_t xEnd, cmp_t xComp ) const
    {
        stxxl::sort( xBegin, xEnd, xComp, 1024 * 1024 * 1024 );
    }
};

template <size_t D, bool dependant_dim>
using OnDiskTypeDef = TypeDefs<default_coordinate_t, //
                              default_val_t, //
                              D, //
                              default_class_key_t, //
                              CachedVecGenerator, //
                              CachedVectorSorter, //
                              dependant_dim, //
                              EXPLAIN, //
                              StdOutProgressStream
                              >;
