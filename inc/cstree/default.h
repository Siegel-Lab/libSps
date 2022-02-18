#pragma once

#include "cstree/tree.h"

using namespace cstree;

template <typename coordinate_t> struct SimpleBinCordsGenerator
{
    std::vector<coordinate_t> binCoords( coordinate_t uiBinCoordsBegin, coordinate_t uiBinCoordsEnd,
                                         coordinate_t uiDim ) const
    {
        std::vector<coordinate_t> vRet;
        for( coordinate_t uiX = uiBinCoordsBegin; uiX <= uiBinCoordsEnd; uiX++ )
            vRet.push_back( uiX );
        return vRet;
    }

    coordinate_t indexForPos( coordinate_t uiBinCoordsBegin, coordinate_t uiBinCoordsEnd, coordinate_t uiDim,
                              coordinate_t uiPos ) const
    {
        std::vector<coordinate_t> vBins = binCoords( uiBinCoordsBegin, uiBinCoordsEnd, uiDim );
        auto iIt = std::lower_bound( vBins.begin( ), vBins.end( ), uiPos );
        if( iIt == vBins.end( ) )
            return std::numeric_limits<coordinate_t>::max( );
        return iIt - vBins.begin( );
    }

    coordinate_t posForIndex( coordinate_t uiBinCoordsBegin, coordinate_t uiBinCoordsEnd, coordinate_t uiDim,
                              coordinate_t uiIndex ) const
    {
        std::vector<coordinate_t> vBins = binCoords( uiBinCoordsBegin, uiBinCoordsEnd, uiDim );
        return vBins[ uiIndex ];
    }

    coordinate_t axisSize( coordinate_t uiBinCoordsBegin, coordinate_t uiBinCoordsEnd, coordinate_t uiDim ) const
    {
        return binCoords( uiBinCoordsBegin, uiBinCoordsEnd, uiDim ).size( );
    }

    coordinate_t minCoord( coordinate_t uiDim ) const
    {}

    coordinate_t maxCoord( coordinate_t uiDim ) const
    {}
};

#if WITH_PYTHON
template <typename type_defs> void exportBinCorrdsGen( pybind11::module& m, std::string sName )
{
    pybind11::class_<SimpleBinCordsGenerator<typename type_defs::coordinate_t>>( m, sName.c_str( ) )
        .def( pybind11::init<>( ) ) // constructor
        ;
}
#endif

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

using default_coordinate_t = uint32_t;
using default_dimensions_t = uint32_t;
using default_cont_sum_val_t = uint32_t;
using default_points_vec_offset_t = uint16_t;
using default_sum_vec_offset_t = uint32_t;

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