#pragma once

#include "kdpstree/tree.h"
#include "kdpstree/psarray.h"

#include <stxxl/io>
#include <stxxl/sort>
#include <stxxl/vector>


using namespace kdpstree;

template <typename val_t> struct TmpVecGenerator
{
    using vec_t = std::vector<val_t>;

    vec_t vec(  )
    {
        return vec_t( );
    }
};

template <typename val_t> struct RamVecGenerator
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

using default_coordinate_t = uint32_t;
using default_val_t = uint32_t;
using default_layer_t = uint8_t;
using default_class_key_t = uint16_t;

template <size_t layers>
using InMemTypeDef = TypeDefs<default_coordinate_t, //
                              default_val_t, //
                              default_layer_t, //
                              layers, //
                              default_class_key_t, //
                              TmpVecGenerator, //
                              RamVecGenerator, //
                              RamVectorSorter, //
                              8, //
                              size_t, //
                              false // explain
                              >;


template <typename val_t> struct DiskVecGenerator
{
    using file_t = stxxl::syscall_file;
    
    using vec_t = typename stxxl::VECTOR_GENERATOR<val_t, 4, 2, sizeof(val_t) * 4 * 1024>::result;

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
                               8, //
                               size_t, //
                               false // explain
                               >;
