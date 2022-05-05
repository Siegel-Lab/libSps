#pragma once

#include "sps/dataset.h"
#include "sps/desc.h"
#include "sps/nd_grid.h"
#include "sps/points.h"
#include "sps/sparse_coordinate.h"
#include "sps/type_defs.h"
#include "sps/util.h"
#include <cassert>
#include <functional>
#include <string>
#include <stxxl/vector>

#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif

namespace sps
{

class AbstractMain
{
    public:
    // make AbstractMain downcastable to Main<something something...> in python by making it abstract
    virtual void dummy() {}
};

template <typename type_defs> class Main: public AbstractMain
{
    EXTRACT_TYPE_DEFS; // macro call

    using point_t = Point<type_defs>;

    using points_t = Points<type_defs>;
    using desc_t = Desc<type_defs>;

    // template<int s> struct CheckSizeOfOverlay;
    // CheckSizeOfOverlay<sizeof(Overlay<type_defs>)> xCheckSizeOfOverlay;

    using overlay_t = AlignedPower2<Overlay<type_defs>>;
    using sparse_coord_t = SparseCoord<type_defs>;
    using prefix_sum_grid_t = NDGrid<type_defs, sps_t>;
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
    Main( std::string sPrefix, bool bWrite = false )
        : vPoints( sPrefix, bWrite ),
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

    template<bool trigger = !IS_ORTHOTOPE>
    typename std::enable_if_t<trigger> addPoint( ret_pos_t vPos, std::string sDesc = "" )
    {
        vPoints.add( vPos, vDesc.add( sDesc ) );
    }

    std::array<pos_t, 2> addDims(ret_pos_t vStart, ret_pos_t vEnd, bool bStartZero) const
    {
        std::array<pos_t, 2> vRet;
        for( size_t uiI = 0; uiI < D - ORTHOTOPE_DIMS; uiI++ )
        {
            vRet[0][uiI] = vStart[uiI];
            vRet[1][uiI] = vEnd[uiI];
        }
        for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
        {
            if(bStartZero)
                vRet[0][uiI + D - ORTHOTOPE_DIMS] = 0;
            else
                vRet[0][uiI + D - ORTHOTOPE_DIMS] = vEnd[uiI] - vStart[uiI];
            vRet[1][uiI + D - ORTHOTOPE_DIMS] = vEnd[uiI] - vStart[uiI];
        }
        return vRet;
    }

    template<bool trigger = IS_ORTHOTOPE>
    typename std::enable_if_t<trigger> addPoint( ret_pos_t vStart, ret_pos_t vEnd, std::string sDesc = "" )
    {
        for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
            if(vStart[uiI] > vEnd[uiI])
                throw std::invalid_argument("end must be larger-equal than start for orthotope dimensions.");
        for( size_t uiI = ORTHOTOPE_DIMS; uiI < D - ORTHOTOPE_DIMS; uiI++ )
            if(vStart[uiI] != vEnd[uiI])
                throw std::invalid_argument("end must equal start for non-orthotope dimensions.");
        auto vP = addDims(vStart, vEnd, false);
        vPoints.add( vP[0], vP[1], vDesc.add( sDesc ) );
    }

    coordinate_t numPoints( ) const
    {
        return vPoints.size( );
    }

    class_key_t generate( coordinate_t uiFrom, coordinate_t uiTo, size_t uiVerbosity=1 )
    {
        progress_stream_t xProg(uiVerbosity);
        typename points_t::Entry xPoints;
        xPoints.uiStartIndex = uiFrom;
        xPoints.uiEndIndex = uiTo;
        class_key_t uiRet = vDataSets.size( );
        vDataSets.push_back( dataset_t( vOverlayGrid, vSparseCoord, vPrefixSumGrid, vPoints, xPoints, xProg ) );
        return uiRet;
    }

    val_t count( class_key_t xDatasetId, ret_pos_t vFromR, ret_pos_t vToR, size_t uiVerbosity=0 ) const
    {
        for( size_t uiI = 0; uiI < D - ORTHOTOPE_DIMS; uiI++ )
            if(vFromR[uiI] > vToR[uiI])
                throw std::invalid_argument("end must be larger-equal than start.");
        progress_stream_t xProg(uiVerbosity);
        auto vP = addDims(vFromR, vToR, true);
        pos_t& vFrom = vP[0];
        pos_t& vTo = vP[1];
        val_t uiRet = 0;
#pragma GCC diagnostic push
// uiD unused if IS_ORTHOTOPE = false
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
        forAllCombinations<pos_t>(
            [ & ]( size_t uiD, pos_t vPos, size_t uiDistToTo ) {
                for( size_t uiI = 0; uiI < D; uiI++ )
                    --vPos[ uiI ];

                xProg << "query: " << xDatasetId << " " << vPos << "\n";
                sps_t vCurrArr = vDataSets[ xDatasetId ].get( vOverlayGrid, vSparseCoord, vPrefixSumGrid, vPos, xProg );

                val_t uiCurr;
                if constexpr(IS_ORTHOTOPE)
                    uiCurr = vCurrArr[uiD / (1 << (D - ORTHOTOPE_DIMS))]; // uiD / 2 ^ (D - ORTHOTOPE_DIMS)
                else
                    uiCurr = vCurrArr;

                xProg << "is " << ( uiDistToTo % 2 == 0 ? "+" : "-" ) << uiCurr << " [" 
                      << uiD << "/" << (1 << (D - ORTHOTOPE_DIMS)) << "]" << "\n";
                uiRet += uiCurr * ( uiDistToTo % 2 == 0 ? 1 : -1 );
            },
            vFrom, vTo, []( coordinate_t uiPos ) { return uiPos > 0; } );
#pragma GCC diagnostic pop
        return uiRet;
    }

    std::string str( ) const
    {
        std::stringstream ss;
        ss << *this;
        return ss.str( );
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
        for( size_t uiI = 0; uiI < rMain.vDataSets.size( ); uiI++ )
            rMain.vDataSets[ uiI ].stream( os, rMain.vOverlayGrid, rMain.vSparseCoord, rMain.vPrefixSumGrid,
                                           rMain.vPoints, rMain.vDesc )
                << std::endl;

        return os;
    }

    std::vector<typename dataset_t::OverlayInfo> getOverlayInfo( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getOverlayInfo( vOverlayGrid, vSparseCoord, vPoints );
    }
};

} // namespace sps


#if WITH_PYTHON
template <typename type_defs> std::string exportMain( pybind11::module& m, std::string sName, std::string sDesc )
{
    using OI = typename sps::Dataset<type_defs>::OverlayInfo;

    pybind11::class_<OI>( m, ( "__" + sName + "_OverlayInfo" ).c_str( ) )
        .def_readonly( "bottom_left", &OI::vBottomLeft )
        .def_readonly( "top_right", &OI::vTopRight )
        .def_readonly( "grid_pos", &OI::vGridPos )
        .def_readonly( "index", &OI::uiIdx )
        .def_readonly( "pred_indices", &OI::vPredIds )
        .def_readonly( "points", &OI::vvPoints )

        ;


    pybind11::class_<sps::Main<type_defs>, sps::AbstractMain> xMain( 
            m, sName.c_str( ), 
            (R"pbdoc(
    Prefix sum index )pbdoc" + sDesc + R"pbdoc(.
    
    .. automethod:: add_point
    .. automethod:: generate
    .. automethod:: count
    .. automethod:: __str__
    .. automethod:: __len__
    .. automethod:: clear

)pbdoc").c_str()
        );
    

    std::string sPointName;
    if(type_defs::ORTHOTOPE_DIMS == 0)
        sPointName = "point";
    else if(type_defs::ORTHOTOPE_DIMS == 1)
        sPointName = "interval";
    else if(type_defs::ORTHOTOPE_DIMS == 2)
        sPointName = "rectangle";
    else if(type_defs::ORTHOTOPE_DIMS == 3)
        sPointName = "cube";
    else
        sPointName = "" + std::to_string(type_defs::ORTHOTOPE_DIMS) + "-orthotope";


    if constexpr(!type_defs::IS_ORTHOTOPE)
        xMain.def( "add_point", 
            [](sps::Main<type_defs>& rM, typename type_defs::pos_t vPos, std::string sDesc){
                rM.addPoint(vPos, sDesc);
            },
            pybind11::arg( "pos" ),
            pybind11::arg( "desc" ) = "",
(R"pbdoc(
    Append a point to the data structure.
    
    :param pos: The position of the point.
    :type pos: list[int[)pbdoc" + std::to_string(type_defs::D - type_defs::ORTHOTOPE_DIMS) + R"pbdoc(]]
    
    :param desc: A description for the Point, defaults to "".
    :type desc: str

    The point will not be queryable until generate is called.
)pbdoc").c_str()
        );
    else
        xMain.def( "add_point", 
            [](sps::Main<type_defs>& rM, typename type_defs::ret_pos_t vStart, 
                typename type_defs::ret_pos_t vEnd, std::string sDesc){
                rM.addPoint(vStart, vEnd, sDesc);
            },
            pybind11::arg( "start" ),
            pybind11::arg( "end" ),
            pybind11::arg( "desc" ) = "",
(R"pbdoc(
    Append a )pbdoc" + sPointName + R"pbdoc( to the data structure.
    
    :param start: The bottom left position of the )pbdoc" + sPointName + R"pbdoc(.
    :type start: list[int[)pbdoc" + std::to_string(type_defs::D - type_defs::ORTHOTOPE_DIMS) + R"pbdoc(]]
    
    :param end: The top right position of the )pbdoc" + sPointName + R"pbdoc(.
    :type end: list[int[)pbdoc" + std::to_string(type_defs::D - type_defs::ORTHOTOPE_DIMS) + R"pbdoc(]]
    
    :param desc: A description for the Point, defaults to "".
    :type desc: str

    The )pbdoc" + sPointName + R"pbdoc( will not be queryable until generate is called.

    Dimensions 1 - )pbdoc" + std::to_string(type_defs::ORTHOTOPE_DIMS) + R"pbdoc( of start and end may have different values, where the value of end must be larger equal to the value of start.
    Dimensions )pbdoc" + std::to_string(type_defs::ORTHOTOPE_DIMS + 1) + " - " + 
    std::to_string(type_defs::D - type_defs::ORTHOTOPE_DIMS) + R"pbdoc( of start and end must have equal values.

    Note that this function will add one point for each outside corner of the given )pbdoc" + sPointName + ".").c_str() );


    xMain.def( pybind11::init<std::string, bool>( ),
              pybind11::arg( "path" ) = "",
              pybind11::arg( "write_mode" ) = false,
R"pbdoc(
    Create a new Index.
    
    :param path: Prefix path of the index on the filesystem (multiple files with different endings will be created), defaults to "".
    :type path: str
    
    :param write_mode: Open the index in write mode (if this is set to False no changes can be made to the index), defaults to True.
    :type write_mode: str
)pbdoc" ) // constructor
        .def( "generate", &sps::Main<type_defs>::generate, pybind11::arg( "from_points" ), pybind11::arg( "to_points" ),
              pybind11::arg( "verbosity" ) = 1,
R"pbdoc(
    Generate a new dataset.
    
    :param from_points: Index of the first point that shall be part of this dataset.
    :type from_points: int
    
    :param to_points: Index of the last point that shall be part of this dataset.
    :type to_points: int
    
    :param verbosity: Degree of verbosity while creating the dataset, defaults to 1.
    :type verbosity: int

    :return: The id of the generated dataset.
    :rtype: int

    This may take a long time to compute.

    Use len(index) to determine the index of the first and last point, as add_point may add multiple points per call.

    This function is multithreaded.
)pbdoc" )
        .def( "count", &sps::Main<type_defs>::count, pybind11::arg( "dataset_id" ), pybind11::arg( "from_pos" ),
                                                     pybind11::arg( "to_pos" ), //
               pybind11::arg( "verbosity" ) = 0 ,
               (R"pbdoc(
    Count the number of points between from and to and in the given dataset.
    
    :param dataset_id: The id of the dataset to query
    :type dataset_id: int
    
    :param from_pos: The bottom left position of the query region.
    :type from_pos: list[int[)pbdoc" + std::to_string(type_defs::D - type_defs::ORTHOTOPE_DIMS) + R"pbdoc(]]
    
    :param to_pos: The top right position of the query region.
    :type to_pos: list[int[)pbdoc" + std::to_string(type_defs::D - type_defs::ORTHOTOPE_DIMS) + R"pbdoc(]]
    
    :param verbosity: Degree of verbosity while counting, defaults to 0.
    :type verbosity: int

    :return: The number of points in dataset_id between from_pos and to_pos.
    :rtype: int

    to_pos must be larger equal than from_pos in each dimension.
)pbdoc").c_str() )
        //.def( "get", &sps::Main<type_defs>::get, "" )
        .def( "__str__", &sps::Main<type_defs>::str, "Return a string describing the index." )
        .def( "__len__", &sps::Main<type_defs>::numPoints, "Total number of points (among all datasets)." )
        .def( "clear", &sps::Main<type_defs>::clear, "Clear the complete index." )
        .def( "__get_overlay_info", &sps::Main<type_defs>::getOverlayInfo )

        ;
    return "    " + sName + "\n";
}
#endif