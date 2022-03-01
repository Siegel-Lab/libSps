#pragma once

#include "kdpstree/tree.h"
#include <map>

using namespace kdpstree;

template <typename val_t> struct RamVecGenerator
{
    std::vector<val_t> operator( )( )
    {
        return std::vector<val_t>( );
    }

    std::vector<val_t> operator( )( std::string )
    {
        return (*this)( );
    }
};

template <typename it_t, typename cmp_t> struct RamVectorSorter
{
    void operator( )( it_t xBegin, it_t xEnd, cmp_t xComp ) const
    {
        std::sort( xBegin, xEnd, xComp );
    }
};
template <typename key_t, typename val_t> struct RamMapGenerator
{
    std::map<key_t, val_t> operator( )( )
    {
        return std::map<key_t, val_t>( );
    }

    std::map<key_t, val_t> operator( )( std::string )
    {
        return (*this)( );
    }
};

using default_coordinate_t = uint32_t;
using default_val_t = uint32_t;
using default_class_key_t = uint16_t;

template <size_t d>
using InMemTypeDef = TypeDefs<default_coordinate_t, //
                              default_val_t, //
                              d, //
                              default_class_key_t, //
                              RamVecGenerator, //
                              RamMapGenerator, //
                              RamVectorSorter //
                              >;