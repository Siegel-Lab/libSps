#pragma once

#include "kdpstree/desc.h"
#include "kdpstree/points.h"
#include "kdpstree/type_defs.h"
#include "kdpstree/sparse_coordinate.h"
#include "kdpstree/nd_grid.h"
#include "kdpstree/dataset.h"
#include <cassert>
#include <functional>
#include <string>

#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif

namespace kdpstree
{

template <typename type_defs> class Main
{
    EXTRACT_TYPE_DEFS; // macro call

    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;
    using desc_t = Desc<type_defs>;
    using overlay_t = Overlay<type_defs>;
    using sparse_coord_t = SparseCoord<type_defs>;
    using prefix_sum_grid_t = NDGrid<type_defs, val_t>;
    using overlay_grid_t = NDGrid<type_defs, overlay_t>;
    using dataset_t = Dataset<type_defs>;
    
    EXTRACT_VEC_GENERATOR( dataset, dataset_t ); // macro call

    
    points_t vPoints;
    desc_t vDesc;
    sparse_coord_t vSparseCoord;
    prefix_sum_grid_t vPrefixSumGrid;
    overlay_grid_t vOverlayGrid;
    
    dataset_file_t xFile;
    dataset_vec_t vDataSets;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter" // do not warn about vFrom and vTo

    template<size_t N>
    inline val_t countHelper( class_key_t uiDatasetIdx, pos_t& vCurr, pos_t vFrom, pos_t vTo, size_t uiDistTo ) const
    {
        if constexpr /* <- required to prevent infinite unrolling loop in compiler */(N == type_defs::D)
        {
            val_t uiRet = vDataSets[uiDatasetIdx].get(vOverlayGrid, vSparseCoord, vPrefixSumGrid, vCurr);
            return uiRet * (uiDistTo % 2 == 0 ? 1 : -1);
        }
        else
        {
            vCurr[N] = vFrom[N];
            val_t uiRet = countHelper<N+1>(uiDatasetIdx, vCurr, vFrom, vTo, uiDistTo + 1);
            vCurr[N] = vTo[N];
            return uiRet + countHelper<N+1>(uiDatasetIdx, vCurr, vFrom, vTo, uiDistTo);
        }
    }
#pragma GCC diagnostic pop

  public:
    Main( std::string sPrefix, bool bWrite = false ) :
        vPoints( sPrefix, bWrite ), 
        vDesc( sPrefix, bWrite ),
        vSparseCoord( sPrefix, bWrite ), 
        vPrefixSumGrid( sPrefix + ".values", bWrite ), 
        vOverlayGrid( sPrefix + ".overlays", bWrite ),
        xFile( dataset_vec_generator.file( sPrefix + ".datsets", bWrite ) ), 
        vDataSets( dataset_vec_generator.vec( xFile ) )
    {}

    void clear( )
    {
        vPoints.clear( );
        vDesc.clear( );
        vSparseCoord.clear( );
        vPrefixSumGrid.clear( );
        vOverlayGrid.clear( );
        vDataSets.clear( );
    }

    void addPoint( class_key_t, pos_t vPos, std::string sDesc )
    {
        vPoints.add( vPos, vDesc.add( sDesc ) );
    }

    coordinate_t numPoints( ) const
    {
        return vPoints.size( );
    }

    class_key_t generate( coordinate_t uiFrom, coordinate_t uiTo, std::optional<progress_stream_t> = { } )
    {
        typename points_t::Entry xPoints;
        xPoints.uiStartIndex = uiFrom;
        xPoints.uiEndIndex = uiTo;
        class_key_t uiRet = vDataSets.size();
        vDataSets.push_back(dataset_t(vOverlayGrid, vSparseCoord, vPrefixSumGrid, vPoints, uiRet, xPoints));
        return uiRet;
    }

    val_t count( class_key_t xDatasetId, pos_t vFrom, pos_t vTo ) const
    {
        pos_t vCurr;
        return countHelper<0>(xDatasetId, vCurr, vFrom, vTo, 0);
    }

};

}


namespace std
{


} // namespace std

#if WITH_PYTHON
template <typename type_defs> void exportStream( pybind11::module& m, std::string sName )
{
    pybind11::class_<typename type_defs::progress_stream_t>( m, sName.c_str( ) );
}
template <typename type_defs> void exportMain( pybind11::module& m, std::string sName )
{
    pybind11::class_<kdpstree::Main<type_defs>>( m, sName.c_str( ) )
        .def( pybind11::init<std::string, bool>( ), pybind11::arg( "path" ),
              pybind11::arg( "write_mode" ) = false ) // constructor
        .def( "add_point", &kdpstree::Main<type_defs>::addPoint )
        .def( "generate", &kdpstree::Main<type_defs>::generate, 
                pybind11::arg( "uiFrom" ), pybind11::arg( "uiTo" ),
                pybind11::arg( "xProg" ) = typename type_defs::progress_stream_t( ) )
        // pybind11::arg( "xProg" ) = std::optional<typename type_defs::progress_stream_t>( ) )
        .def( "count", &kdpstree::Main<type_defs>::count )
        //.def( "get", &kdpstree::Main<type_defs>::get, "" )
        //.def( "__str__", &kdpstree::Main<type_defs>::print )
        .def( "__len__", &kdpstree::Main<type_defs>::numPoints )
        .def( "clear", &kdpstree::Main<type_defs>::clear )

        ;
}
#endif