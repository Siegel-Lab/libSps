#pragma once

#include "cstree/tree.h"

using default_coordinate_t = uint32_t;
using default_dimensions_t = uint32_t;
using default_cont_sum_val_t = uint32_t;
using default_points_vec_offset_t = uint16_t;
using default_sum_vec_offset_t = uint32_t;

using namespace cstree;

template <typename coordinate_t> struct SimpleBinCordsGenerator
{
    std::vector<coordinate_t> operator( )( coordinate_t uiBegin, coordinate_t uiEnd, coordinate_t uiDim,
                                           coordinate_t uiLayer )
    {
        std::vector<coordinate_t> vRet;
        for( coordinate_t uiX = uiBegin; uiX <= uiEnd; uiX++ )
            vRet.push_back( uiX );
        return vRet;
    }
};

template <typename vec_t> struct RamVecGenerator
{
    vec_t operator( )( std::string )
    {
        return vec_t( );
    }
};

template <typename it_t, typename cmp_t> struct RamVectorSorter
{
    void operator( )( it_t xBegin, it_t xEnd, cmp_t xComp )
    {
        std::sort( xBegin, xEnd, xComp );
    }
};

template <size_t d>
using InMemTypeDef = TypeDefs<default_coordinate_t, //
                              default_cont_sum_val_t, //
                              default_dimensions_t, //
                              default_points_vec_offset_t, //
                              default_sum_vec_offset_t, //
                              d, //
                              SimpleBinCordsGenerator, //
                              RamVecGenerator, //
                              std::vector, //
                              std::vector, //
                              std::vector, //
                              std::vector, //
                              RamVectorSorter //
                              >;