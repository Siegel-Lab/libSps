#pragma once

#include "sps/util.h"
#include "sps/desc.h"
#include "sps/points.h"
#include "sps/type_defs.h"
#include "sps/sparse_coordinate.h"
#include "sps/nd_grid.h"
#include "sps/dataset.h"
#include <stxxl/vector>
#include <cassert>
#include <functional>
#include <string>

#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif

namespace sps
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
    inline val_t countHelper( class_key_t uiDatasetIdx, pos_t& vCurr, pos_t vFrom, pos_t vTo, size_t uiDistTo,
                              std::optional<progress_stream_t> xProg = { } ) const
    {
        if constexpr /* <- required to prevent infinite unrolling loop in compiler */(N == type_defs::D)
        {
            xProg << "query: " << uiDatasetIdx << " " << vCurr << "\n";
            val_t uiRet = vDataSets[uiDatasetIdx].get(vOverlayGrid, vSparseCoord, vPrefixSumGrid, vCurr, xProg);

            xProg << "is " << (uiDistTo % 2 == 0 ? "+" : "-") << uiRet << "\n";
            return uiRet * (uiDistTo % 2 == 0 ? 1 : -1);
        }
        else
        {
            val_t uiRet = 0;
            if(vFrom[N] > 0)
            {
                vCurr[N] = vFrom[N] - 1;
                uiRet += countHelper<N+1>(uiDatasetIdx, vCurr, vFrom, vTo, uiDistTo + 1, xProg);
            }
            if(vTo[N] > 0)
            {
                vCurr[N] = vTo[N] - 1;
                uiRet += countHelper<N+1>(uiDatasetIdx, vCurr, vFrom, vTo, uiDistTo, xProg);
            }
            return uiRet;
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

    void addPoint( pos_t vPos, std::string sDesc )
    {
        vPoints.add( vPos, vDesc.add( sDesc ) );
    }

    coordinate_t numPoints( ) const
    {
        return vPoints.size( );
    }

    class_key_t generate( coordinate_t uiFrom, coordinate_t uiTo, std::optional<progress_stream_t> xProg = { } )
    {
        typename points_t::Entry xPoints;
        xPoints.uiStartIndex = uiFrom;
        xPoints.uiEndIndex = uiTo;
        class_key_t uiRet = vDataSets.size();
        vDataSets.push_back(dataset_t(vOverlayGrid, vSparseCoord, vPrefixSumGrid, vPoints, xPoints, xProg));
        return uiRet;
    }

    val_t count( class_key_t xDatasetId, pos_t vFrom, pos_t vTo, std::optional<progress_stream_t> xProg = { } ) const
    {
        pos_t vCurr;
        return countHelper<0>(xDatasetId, vCurr, vFrom, vTo, 0, xProg);
    }

    std::string str() const
    {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }

    friend std::ostream& operator<<( std::ostream& os, const Main& rMain )
    {
        os << "vDataSets: ";
        os << rMain.vDataSets << std::endl;
        os << "vOverlayGrid: ";
        os << rMain.vOverlayGrid << std::endl;
        os << "vSparseCoord: ";
        os << rMain.vSparseCoord << std::endl;
        os << "vPrefixSumGrid: ";
        os << rMain.vPrefixSumGrid << std::endl;
        os << "vPoints: ";
        os << rMain.vPoints << std::endl;
        os << "vDesc: ";
        os << rMain.vDesc << std::endl;

        os << "Pretty Print: ";
        for(size_t uiI = 0; uiI < rMain.vDataSets.size(); uiI++)
            rMain.vDataSets[uiI].stream(os, rMain.vOverlayGrid, rMain.vSparseCoord, rMain.vPrefixSumGrid,
                                        rMain.vPoints, rMain.vDesc) << std::endl;

        return os;
    }
};

}


#if WITH_PYTHON
template <typename type_defs> void exportStream( pybind11::module& m, std::string sName )
{
    pybind11::class_<typename type_defs::progress_stream_t>( m, sName.c_str( ) );
}
template <typename type_defs> void exportMain( pybind11::module& m, std::string sName )
{
    pybind11::class_<sps::Main<type_defs>>( m, sName.c_str( ) )
        .def( pybind11::init<std::string, bool>( ), 
              pybind11::arg( "path" ),
              pybind11::arg( "write_mode" ) = false ) // constructor
        .def( "add_point", &sps::Main<type_defs>::addPoint )
        .def( "generate", &sps::Main<type_defs>::generate, 
                pybind11::arg( "uiFrom" ), pybind11::arg( "uiTo" ),
                pybind11::arg( "xProg" ) = typename type_defs::progress_stream_t( ) )
        // pybind11::arg( "xProg" ) = std::optional<typename type_defs::progress_stream_t>( ) )
        .def( "count", &sps::Main<type_defs>::count, 
                pybind11::arg( "uiDatasetId" ), pybind11::arg( "uiFrom" ), pybind11::arg( "uiTo" ),
                pybind11::arg( "xProg" ) = typename type_defs::progress_stream_t( ) )
        //.def( "get", &sps::Main<type_defs>::get, "" )
        .def( "__str__", &sps::Main<type_defs>::str )
        .def( "__len__", &sps::Main<type_defs>::numPoints )
        .def( "clear", &sps::Main<type_defs>::clear )

        ;
}
#endif