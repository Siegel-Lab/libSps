#pragma once

#include "kdpstree/psarray.h"
#include "kdpstree/tree.h"

#include <stxxl/io>
#include <stxxl/sort>
#include <stxxl/vector>
#include <chrono>
#include <iostream>


using namespace kdpstree;

struct StdOutProgressStream
{
    std::chrono::time_point<std::chrono::high_resolution_clock> xLastPrint{};

    bool printAgain()
    {
        auto xCurr = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> xTime = xCurr - xLastPrint;
        if(xTime.count() <= 1000)
            return false;
        xLastPrint = xCurr;
        return true;
    }

    template<typename T>
    StdOutProgressStream& operator<<(const T& sStr)
    {
        std::cout << sStr << std::flush;
        return *this;
    }

    StdOutProgressStream() :
        xLastPrint(std::chrono::high_resolution_clock::now())
    {}
};

template <typename val_t> struct TmpVecGenerator
{
    using vec_t = std::vector<val_t>;

    vec_t vec( )
    {
        return vec_t( );
    }
};

template <typename val_t, size_t ele_per_block> struct RamVecGenerator
{
    using file_t = size_t;
    using vec_t = std::vector<val_t>;

    vec_t vec( file_t& )
    {
        return vec_t( );
    }

    file_t file( std::string )
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

static const bool EXPLAIN = false;
static const size_t B = 170;

using default_coordinate_t = uint32_t;
using default_val_t = uint32_t;
using default_layer_t = uint8_t;
using default_class_key_t = uint16_t; // @todo should be uint8_t

template <size_t layers>
using InMemTypeDef = TypeDefs<default_coordinate_t, //
                              default_val_t, //
                              default_layer_t, //
                              layers, //
                              default_class_key_t, //
                              TmpVecGenerator, //
                              RamVecGenerator, //
                              RamVectorSorter, //
                              B, //
                              size_t, //
                              EXPLAIN, // explain
                              StdOutProgressStream
                              >;
                            
static const size_t B_TEST = 3;
template <size_t layers>
using TestTypeDef = TypeDefs<default_coordinate_t, //
                              default_val_t, //
                              default_layer_t, //
                              layers, //
                              default_class_key_t, //
                              TmpVecGenerator, //
                              RamVecGenerator, //
                              RamVectorSorter, //
                              B_TEST, //
                              size_t, //
                              EXPLAIN, // explain
                              StdOutProgressStream
                              >;


template <typename val_t, size_t ele_per_block> struct DiskVecGenerator
{
    using file_t = stxxl::syscall_file;

    using vec_t =
        typename stxxl::VECTOR_GENERATOR<val_t, 1, ( 1024 * 1024 * 1024 ) / ( sizeof( val_t ) * ele_per_block ),
                                         sizeof( val_t ) * ele_per_block>::result;

    vec_t vec( file_t& rFile )
    {
        return vec_t( &rFile );
    }

    file_t file( std::string sPath )
    {
        return file_t( sPath, stxxl::file::RDWR | stxxl::file::CREAT | stxxl::file::DIRECT );
    }
};

template <typename it_t, typename cmp_t> struct DiskVectorSorter
{
    void operator( )( it_t xBegin, it_t xEnd, cmp_t xComp ) const
    {
        stxxl::sort( xBegin, xEnd, xComp, 1024 * 1024 * 1024 );
    }
};

template <size_t layers>
using OnDiskTypeDef = TypeDefs<default_coordinate_t, //
                               default_val_t, //
                               default_layer_t, //
                               layers, //
                               default_class_key_t, //
                               TmpVecGenerator, //
                               DiskVecGenerator, //
                               DiskVectorSorter, //
                               B, //
                               size_t, //
                               EXPLAIN, // explain
                               StdOutProgressStream
                               >;
