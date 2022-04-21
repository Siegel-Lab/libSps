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
    
    //template<int s> struct CheckSizeOfOverlay;
    //CheckSizeOfOverlay<sizeof(Overlay<type_defs>)> xCheckSizeOfOverlay;

    using overlay_t = AlignedPower2<Overlay<type_defs>>;
    using sparse_coord_t = SparseCoord<type_defs>;
    using prefix_sum_grid_t = NDGrid<type_defs, val_t>;
    using overlay_grid_t = NDGrid<type_defs, overlay_t>;
    using dataset_t = AlignedPower2<Dataset<type_defs>>;

    EXTRACT_VEC_GENERATOR( dataset, dataset_t ); // macro call

    points_t vPoints;
    desc_t vDesc;
    sparse_coord_t vSparseCoord;
    prefix_sum_grid_t vPrefixSumGrid;
    overlay_grid_t vOverlayGrid;
    
    dataset_file_t xFile;
    dataset_vec_t vDataSets;


  public:
    Main( std::string sPrefix, bool bWrite = false ) :
        vPoints( sPrefix, bWrite ), 
        vDesc( sPrefix, bWrite ),
        vSparseCoord( sPrefix, bWrite ), 
        vPrefixSumGrid( sPrefix + ".prefix_sums", bWrite ), 
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

    class_key_t generate( coordinate_t uiFrom, coordinate_t uiTo, progress_stream_t& xProg )
    {
        typename points_t::Entry xPoints;
        xPoints.uiStartIndex = uiFrom;
        xPoints.uiEndIndex = uiTo;
        class_key_t uiRet = vDataSets.size();
        vDataSets.push_back(dataset_t(vOverlayGrid, vSparseCoord, vPrefixSumGrid, vPoints, xPoints, xProg));
        return uiRet;
    }

    val_t count( class_key_t xDatasetId, pos_t vFrom, pos_t vTo, progress_stream_t& xProg ) const
    {
        val_t uiRet = 0;
        forAllCombinations<pos_t>([&](pos_t vPos, size_t uiDistToTo){
            for(size_t uiI = 0; uiI < D; uiI++)
                --vPos[uiI];

            xProg << "query: " << xDatasetId << " " << vPos << "\n";
            val_t uiCurr = vDataSets[xDatasetId].get(vOverlayGrid, vSparseCoord, vPrefixSumGrid, vPos, xProg);

            xProg << "is " << (uiDistToTo % 2 == 0 ? "+" : "-") << uiCurr << "\n";
            uiRet += uiCurr * (uiDistToTo % 2 == 0 ? 1 : -1);

        }, vFrom, vTo, [](coordinate_t uiPos){ return uiPos > 0; } );
        return uiRet;
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

    std::vector<typename dataset_t::OverlayInfo> getOverlayInfo(class_key_t xDatasetId) const
    {
        return vDataSets[xDatasetId].getOverlayInfo(vOverlayGrid, vSparseCoord, vPoints); 
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
    using OI = typename sps::Dataset<type_defs>::OverlayInfo;

    pybind11::class_<OI>( m, ("__"+sName+"_OverlayInfo").c_str( ) )
        .def_readonly("bottom_left", &OI::vBottomLeft)
        .def_readonly("top_right", &OI::vTopRight)
        .def_readonly("grid_pos", &OI::vGridPos)
        .def_readonly("index", &OI::uiIdx)
        .def_readonly("pred_indices", &OI::vPredIds)
        .def_readonly("points", &OI::vvPoints)

        ;

    pybind11::class_<sps::Main<type_defs>>( m, sName.c_str( ) )
        .def( pybind11::init<std::string, bool>( ), 
              pybind11::arg( "path" ),
              pybind11::arg( "write_mode" ) = false ) // constructor
        .def( "add_point", &sps::Main<type_defs>::addPoint )
        .def( "generate", &sps::Main<type_defs>::generate, 
                pybind11::arg( "uiFrom" ), pybind11::arg( "uiTo" ),
                pybind11::arg( "xProg" ) = typename type_defs::progress_stream_t( 1 ) )
                //pybind11::arg( "xProg" ) = typename type_defs::progress_stream_t( 5 ) )
        // pybind11::arg( "xProg" ) = std::optional<typename type_defs::progress_stream_t>( ) )
        .def( "count", &sps::Main<type_defs>::count, 
                pybind11::arg( "uiDatasetId" ), pybind11::arg( "uiFrom" ), pybind11::arg( "uiTo" ),
                pybind11::arg( "xProg" ) = typename type_defs::progress_stream_t( 0 ) )
                //pybind11::arg( "xProg" ) = typename type_defs::progress_stream_t( 5 ) )
        //.def( "get", &sps::Main<type_defs>::get, "" )
        .def( "__str__", &sps::Main<type_defs>::str )
        .def( "__len__", &sps::Main<type_defs>::numPoints )
        .def( "clear", &sps::Main<type_defs>::clear )
        .def( "get_overlay_info", &sps::Main<type_defs>::getOverlayInfo )

        ;
}
#endif