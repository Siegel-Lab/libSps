#include "cstree/tree.h"

using default_coordinate_t = uint32_t;
using default_cont_sum_val_t = uint32_t;
using default_points_vec_offset_t = uint16_t;
using default_sum_vec_offset_t = uint32_t;

template <typename coordinate_t>
std::vector<coordinate_t> simpleBinCordsGenerator( coordinate_t uiBegin, coordinate_t uiEnd, coordinate_t uiDim,
                                                   coordinate_t uiLayer )
{
    std::vector<coordinate_t> vRet;
    for( coordinate_t uiX = uiBegin; uiX <= uiEnd; uiX++ )
        vRet.push_back( uiX );
    return vRet;
}

template <typename vec_t> vec_t inMemVecGenerator( std::string )
{
    return vec_t( );
}

template <size_t d>
using InMemTypeDef = TypeDefs<default_coordinate_t, //
                              default_cont_sum_val_t, //
                              default_points_vec_offset_t, //
                              default_sum_vec_offset_t, //
                              d, //
                              simpleBinCordsGenerator, //
                              inMemVecGenerator, //
                              std::vector, //
                              std::vector, //
                              std::vector, //
                              std::vector, //
                              std::sort //
                              >;

PYBIND11_MODULE( libcstree, m )
{
    exportTree<InMemTypeDef>( m, "CsTree" );
}