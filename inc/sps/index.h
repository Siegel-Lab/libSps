/**
 * @file index.h
 * @brief Implements the main Datastructure.
 * @author Markus Schmidt
 */


#pragma once

#include "sps/abstract_index.h"
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

#ifdef DO_PROFILE
#include <gperftools/profiler.h>
#endif

namespace sps
{

/**
 * @brief An Enum for Querying the index.
 *
 * Which orthotopes to count, depending on how they intersect the queried area.
 * Only relevant for the Index.count() function.
 */
enum IntersectionType
{
    /// @brief count orthotopes that are fully enclosed by the queried area
    enclosed,
    /// @brief count orthotopes that fully enclose by the queried area
    encloses,
    /// @brief count orthotopes that overlap the queried area
    overlaps,
    /// @brief count orthotopes that have their bottom-left-front-.. corner in the queried area
    first,
    /// @brief count orthotopes that have their top-right-back-.. corner in the queried area
    last,
    /// @brief count orthotopes that are point-like and in the queried area
    points_only,
    /// @brief place the orthotope at a position accroding to the size in its orthotope dimensions (used for insert)
    slice,
};

/**
 * @brief The main sparse prefix sum index class.
 *
 * @tparam type_defs An instance of TypeDefs, that defines all compiletime parameters
 */
template <typename type_defs> class Index : public AbstractIndex
{
  public:
    EXTRACT_TYPE_DEFS; // macro call
  private:
    using point_t = Point<type_defs>;

    using points_t = Points<type_defs>;
    using desc_t = Desc<type_defs>;

    // template<int s> struct CheckSizeOfOverlay;
    // CheckSizeOfOverlay<sizeof(Overlay<type_defs>)> xCheckSizeOfOverlay;

    using overlay_t = AlignedPower2<Overlay<type_defs>>;
    using sparse_coord_t = SparseCoord<type_defs>;
    using overlay_grid_t = typename Overlay<type_defs>::overlay_grid_t;
    using prefix_sum_grid_t = typename Overlay<type_defs>::prefix_sum_grid_t;
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
    /**
     * @brief Construct a new Index object
     *
     * @param sPrefix Prefix path of the index on the filesystem (multiple files with different endings will be
     * created), defaults to "".
     * @param bWrite Open the index in write mode (if this is set to False no changes can be made to the index),
     * defaults to false.
     */
    Index( std::string sPrefix = "", bool bWrite = false )
        : vPoints( sPrefix, bWrite ),
          vDesc( sPrefix, bWrite ),
          vSparseCoord( sPrefix, bWrite ),
          vPrefixSumGrid( sPrefix + ".prefix_sums", bWrite ),
          vOverlayGrid( sPrefix + ".overlays", bWrite ),
          xFile( dataset_vec_generator.file( sPrefix + ".datsets", bWrite ) ),
          vDataSets( dataset_vec_generator.vec( xFile ) )
    {}

    /**
     * @brief Clear the complete index.
     *
     * Clears all datasets.
     */
    void clear( )
    {
        vPoints.clear( );
        vDesc.clear( );
        vSparseCoord.clear( );
        vPrefixSumGrid.clear( );
        vOverlayGrid.clear( );
        vDataSets.clear( );
    }

    /**
     * @brief Append a point to the data structure.
     *
     * The point will not be queryable until generate is called.
     *
     * @tparam trigger This function is only active if IS_ORTHOTOPE = false
     * @param vPos The position of the point.
     * @param sDesc A description for the Point, defaults to "".
     */
    template <bool trigger = !IS_ORTHOTOPE>
    typename std::enable_if_t<trigger> addPoint( ret_pos_t vPos, std::string sDesc = "" )
    {
        vPoints.add( vPos, vDesc.add( sDesc ) );
    }

    std::array<pos_t, 2> addDims( ret_pos_t vStart, ret_pos_t vEnd,
                                  IntersectionType xIntersectionType = IntersectionType::slice ) const
    {
        std::array<pos_t, 2> vRet;
        for( size_t uiI = 0; uiI < D - ORTHOTOPE_DIMS; uiI++ )
        {
            vRet[ 0 ][ uiI ] = vStart[ uiI ];
            vRet[ 1 ][ uiI ] = vEnd[ uiI ];
        }
        for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
        {
            if( xIntersectionType == IntersectionType::slice || xIntersectionType == IntersectionType::encloses )
                vRet[ 0 ][ uiI + D - ORTHOTOPE_DIMS ] = vEnd[ uiI ] - vStart[ uiI ];
            else
                vRet[ 0 ][ uiI + D - ORTHOTOPE_DIMS ] = 0;

            if( xIntersectionType == IntersectionType::slice || xIntersectionType == IntersectionType::enclosed )
                vRet[ 1 ][ uiI + D - ORTHOTOPE_DIMS ] = vEnd[ uiI ] - vStart[ uiI ];
            else if( xIntersectionType == IntersectionType::points_only )
                vRet[ 1 ][ uiI + D - ORTHOTOPE_DIMS ] = 1;
            else
                vRet[ 1 ][ uiI + D - ORTHOTOPE_DIMS ] = std::numeric_limits<coordinate_t>::max( );
        }
        return vRet;
    }

    /**
     * @brief Append an orthotope to the data structure.
     *
     * The orthotope will not be queryable until generate is called.
     *
     * @tparam trigger This function is only active if IS_ORTHOTOPE = true
     * @param vStart The bottom left position of the orthotope.
     * @param vEnd The top right position of the orthotope.
     * @param sDesc A description for the orthotope, defaults to "".
     */
    template <bool trigger = IS_ORTHOTOPE>
    typename std::enable_if_t<trigger> addPoint( ret_pos_t vStart, ret_pos_t vEnd, std::string sDesc = "" )
    {
        for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
            if( vStart[ uiI ] > vEnd[ uiI ] )
                throw std::invalid_argument( "end must be larger-equal than start for orthotope dimensions." );
        for( size_t uiI = ORTHOTOPE_DIMS; uiI < D - ORTHOTOPE_DIMS; uiI++ )
            if( vStart[ uiI ] != vEnd[ uiI ] )
                throw std::invalid_argument( "end must equal start for non-orthotope dimensions." );
        auto vP = addDims( vStart, vEnd );
        vPoints.add( vP[ 0 ], vP[ 1 ], vDesc.add( sDesc ) );
    }

    /**
     * @brief Total number of points (among all datasets).
     */
    coordinate_t numPoints( ) const
    {
        return vPoints.size( );
    }

    /**
     * @brief Generate a new dataset.
     *
     * This may take a long time to compute.
     *
     * Use len(index) to determine the index of the first and last point, as add_point may add multiple points per call.
     *
     * This function is multithreaded.
     *
     * @param uiFrom Size of the index before adding the first point of this dataset, defaults to zero.
     * @param uiTo Size of the index after adding all points of this dataset, defaults to the current size of the index.
     * @param uiVerbosity Degree of verbosity while creating the dataset, defaults to 1.
     * @return class_key_t The id of the generated dataset.
     */
    class_key_t generate( coordinate_t uiFrom = 0,
                          coordinate_t uiTo = std::numeric_limits<coordinate_t>::max( ),
                          double fFac = 0,
                          size_t uiVerbosity = 1 )
    {
#ifdef DO_PROFILE
        ProfilerStart( "gperftools.generate.prof" );
#endif
        if( uiTo == std::numeric_limits<coordinate_t>::max( ) )
            uiTo = numPoints( );

        progress_stream_t xProg( uiVerbosity );
        typename points_t::Entry xPoints;
        xPoints.uiStartIndex = uiFrom;
        xPoints.uiEndIndex = uiTo;
        // generate the dataset in ram then push it into the index to make sure that the cache of the vector
        // does not unload the memory half way through the initialization. (not relevant for std::vector
        // implementations)
        dataset_t xNew( vOverlayGrid, vSparseCoord, vPrefixSumGrid, vPoints, xPoints, fFac, xProg );
        class_key_t uiRet = vDataSets.size( );
        vDataSets.push_back( xNew );

#ifdef DO_PROFILE
        ProfilerFlush( );
        ProfilerStop( );
#endif
#ifndef NDEBUG
        xProg << Verbosity( 1 ) << "\n\nMaximal prefix sum value: " << maxPrefixSumValue( ) << ".\n";
#endif
        return uiRet;
    }

  private:
#pragma GCC diagnostic push
// vPosTopRightActual not used with DEPENDANT_DIMENSION == false
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
    static inline __attribute__( ( always_inline ) ) void
    countSizeLimitedInvariant( size_t uiD, pos_t vPos, size_t uiDistToTo, IntersectionType xInterType,
                               const dataset_vec_t& vDataSets, const sparse_coord_t& vSparseCoord,
                               const prefix_sum_grid_t& vPrefixSumGrid, const overlay_grid_t& vOverlayGrid,
                               const class_key_t xDatasetId, val_t& uiRet )
    {
        for( size_t uiI = 0; uiI < D; uiI++ )
            --vPos[ uiI ];

#if GET_PROG_PRINTS
        xProg << "query: " << xDatasetId << " " << vPos << "\n";
#endif

        size_t uiCornerIndex = 0;
        if constexpr( IS_ORTHOTOPE )
        {
            size_t uiDLookup;
            switch( xInterType )
            {
                default:
                case IntersectionType::points_only:
                case IntersectionType::first:
                    uiDLookup = 0;
                    break;
                case IntersectionType::last:
                    uiDLookup = ( ( 1 << D ) - 1 );
                    break;
                case IntersectionType::enclosed:
                case IntersectionType::encloses:
                    uiDLookup = uiD;
                    break;
                case IntersectionType::overlaps:
                    // invert the last D bits of uiD
                    // set all other bits to 0
                    uiDLookup = ( ~uiD ) & ( ( 1 << D ) - 1 );
                    break;
            }
            // uiDLookup / 2 ^ (D - ORTHOTOPE_DIMS)
            uiCornerIndex = uiDLookup / ( 1 << ( D - ORTHOTOPE_DIMS ) );
        }
        val_t uiCurr = vDataSets[ xDatasetId ].get( vOverlayGrid, vSparseCoord, vPrefixSumGrid, vPos, uiCornerIndex
#if GET_PROG_PRINTS
                                                    ,
                                                    xProg
#endif
        );

        val_t uiFac = ( uiDistToTo % 2 == 0 ? 1 : -1 );
#if GET_PROG_PRINTS
        xProg << "is " << ( uiFac == 1 ? "+" : "-" ) << uiCurr << " [" << uiD << "/" << ( 1 << ( D - ORTHOTOPE_DIMS ) )
              << "]"
              << "\n";
#endif
        uiRet += uiCurr * uiFac;
    }

#pragma GCC diagnostic pop

    static inline __attribute__( ( always_inline ) ) bool
    countSizeLimitedInvariantCond( coordinate_t uiPos, size_t /*uiD*/, bool /*bIsFrom*/ )
    {
        return uiPos > 0;
    }

  public:
    /**
     * @brief Count the number of points between from and to in the given dataset.
     *
     * As opposed to count, this function allows specifying the start and end positions for all dimensions in the
     * datastructure.
     * This is only relevant for indices with orthotope dimensions.
     *
     * @param xDatasetId The id of the dataset to query
     * @param vFrom The bottom left position of the query region.
     * @param vTo The top right position of the query region.
     * @param xInterType The used intersection type, defaults to enclosed.
     * @param uiVerbosity Degree of verbosity while counting, defaults to 0.
     * @return val_t The number of points in dataset_id between from_pos and to_pos.
     */
    val_t countSizeLimited( class_key_t xDatasetId, pos_t vFrom, pos_t vTo,
                            IntersectionType xInterType = IntersectionType::enclosed,
                            size_t
#if GET_PROG_PRINTS
                                uiVerbosity
#endif
                            = 0 ) const
    {
#if GET_PROG_PRINTS
        progress_stream_t xProg( uiVerbosity );
#endif
        val_t uiRet = 0;
#pragma GCC diagnostic push
// uiD unused if IS_ORTHOTOPE = false
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
        forAllCombinationsTmpl<pos_t>( vFrom, vTo, countSizeLimitedInvariant, countSizeLimitedInvariantCond, xInterType,
                                       vDataSets, vSparseCoord, vPrefixSumGrid, vOverlayGrid, xDatasetId, uiRet );
#pragma GCC diagnostic pop
        return uiRet;
    }

    /**
     * @brief Count the number of points between from and to and in the given dataset.
     *
     * to_pos must be larger equal than from_pos in each dimension.
     *
     * @param xDatasetId The id of the dataset to query
     * @param vFromR The bottom left position of the query region.
     * @param vToR The top right position of the query region.
     * @param xInterType The used intersection type, defaults to enclosed. Ignored if there are no orthotope dimensions.
     * @param uiVerbosity Degree of verbosity while counting, defaults to 0.
     * @return val_t The number of points in dataset_id between from_pos and to_pos.
     */
    val_t count( class_key_t xDatasetId, ret_pos_t vFromR, ret_pos_t vToR,
                 IntersectionType xInterType = IntersectionType::enclosed, size_t uiVerbosity = 0 ) const
    {
        for( size_t uiI = 0; uiI < D - ORTHOTOPE_DIMS; uiI++ )
            if( vFromR[ uiI ] > vToR[ uiI ] )
                throw std::invalid_argument( "end must be larger-equal than start." );
        auto vP = addDims( vFromR, vToR, xInterType );
        pos_t vFrom = vP[ 0 ];
        pos_t vTo = vP[ 1 ];
        return countSizeLimited( xDatasetId, vFrom, vTo, xInterType, uiVerbosity );
    }

    /**
     * @brief Count the number of points between from and to and in the given dataset.
     *
     * Counts for multiple regions.
     * to_pos must be larger equal than from_pos in each dimension.
     *
     * @param xDatasetId The id of the dataset to query
     * @param vRegions The bottom left and top right positions of the queried regions.
     * @param xInterType The used intersection type, defaults to enclosed. Ignored if there are no orthotope dimensions.
     * @param uiVerbosity Degree of verbosity while counting, defaults to 0.
     * @return std::vector<val_t> The number of points in dataset_id between from_pos and to_pos for each given region.
     */
    std::vector<val_t> countMultiple( class_key_t xDatasetId, std::vector<std::pair<ret_pos_t, ret_pos_t>> vRegions,
                                      IntersectionType xInterType = IntersectionType::enclosed,
                                      size_t uiVerbosity = 0 ) const
    {
#ifdef DO_PROFILE
        ProfilerStart( "gperftools.countMultiple.prof" );
#endif

        std::vector<val_t> vRet( vRegions.size( ) );

        for( size_t uiI = 0; uiI < vRegions.size( ); uiI++ )
            vRet[ uiI ] = count( xDatasetId, vRegions[ uiI ].first, vRegions[ uiI ].second, xInterType, uiVerbosity );
#ifdef DO_PROFILE
        ProfilerFlush( );
        ProfilerStop( );
#endif
        return vRet;
    }

    val_t maxPrefixSumValue( ) const
    {
        val_t uiMax = 0;
        for( const sps_t& rSps : vPrefixSumGrid.vData )
        {
            if constexpr( IS_ORTHOTOPE )
                for( size_t uiD = 0; uiD < ORTHOTOPE_DIMS; uiD++ )
                    uiMax = std::max( uiMax, rSps[ uiD ] );
            else
                uiMax = std::max( uiMax, rSps );
        }
        return uiMax;
    }

    coordinate_t getNumInternalPrefixSums( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getNumInternalPrefixSums( vOverlayGrid, vSparseCoord );
    }

    coordinate_t getNumOverlayPrefixSums( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getNumOverlayPrefixSums( vOverlayGrid, vSparseCoord );
    }

    coordinate_t getNumInternalSparseCoords( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getNumInternalSparseCoords( vOverlayGrid );
    }

    coordinate_t getNumOverlaySparseCoords( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getNumOverlaySparseCoords( vOverlayGrid );
    }

    static std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>
    estimateDataStructureElements( pos_t uiCoordinateSizes, size_t uiNumPoints, double fFac )
    {
        return dataset_t::estimateDataStructureElements( uiCoordinateSizes, uiNumPoints, fFac );
    }

    static uint64_t estimateDataStructureSize( pos_t uiCoordinateSizes, size_t uiNumPoints, double fFac )
    {
        return dataset_t::estimateDataStructureSize( uiCoordinateSizes, uiNumPoints, fFac );
    }

    static uint64_t pickNumOverlays( pos_t uiCoordinateSizes, size_t uiNumPoints )
    {
        return dataset_t::pickNumOverlays( uiCoordinateSizes, uiNumPoints );
    }
    
    static double toFactor( pos_t uiCoordinateSizes, uint64_t uiNumOverlays)
    {
        return dataset_t::toFactor( uiCoordinateSizes, uiNumOverlays );
    }

    /**
     * @brief Return a string describing the index.
     *
     * Very slow for large datasets.
     */
    std::string str( ) const
    {
        std::stringstream ss;
        ss << *this;
        return ss.str( );
    }

    friend std::ostream& operator<<( std::ostream& os, const Index& rMain )
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

    /**
     * @brief Returns the bottom-left-front-... and top-right-back-... position of all overlays.
     */
    std::vector<std::array<pos_t, 3>> getOverlayGrid( class_key_t xDatasetId ) const
    {
        if( vDataSets.size( ) <= xDatasetId )
            return std::vector<std::array<pos_t, 3>>{ };
        return vDataSets[ xDatasetId ].getOverlayGrid( vOverlayGrid, vSparseCoord );
    }
};

} // namespace sps


#if WITH_PYTHON
void exportEnum( pybind11::module& m )
{
    pybind11::enum_<IntersectionType>( m, "IntersectionType", R"pbdoc(
    An Enum for Querying the index.
    
    Which orthotopes to count, depending on how they intersect the queried area.
    Only relevant for the Index.count() function.
)pbdoc" )
        .value( "enclosed", IntersectionType::enclosed, "count orthotopes that are fully enclosed by the queried area" )
        .value( "encloses", IntersectionType::encloses, "count orthotopes that fully enclose by the queried area" )
        .value( "overlaps", IntersectionType::overlaps, "count orthotopes that overlap the queried area" )
        .value( "first", IntersectionType::first,
                "count orthotopes that have their bottom-left-front-.. corner in the queried area" )
        .value( "last", IntersectionType::last,
                "count orthotopes that have their top-right-back-.. corner in the queried area" )
        .value( "points_only", IntersectionType::points_only,
                "count orthotopes that are point-like and in the queried area" )
        .export_values( );
}

template <typename type_defs> std::string exportIndex( pybind11::module& m, std::string sName, std::string sDesc )
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


    pybind11::class_<sps::Index<type_defs>, sps::AbstractIndex> xMain( m, sName.c_str( ),
                                                                       ( R"pbdoc(
    Prefix sum index )pbdoc" + sDesc + R"pbdoc(.
    
    .. automethod:: add_point
    .. automethod:: generate
    .. automethod:: count
    .. automethod:: count_multiple
    .. automethod:: __str__
    .. automethod:: __len__
    .. automethod:: clear
)pbdoc" )
                                                                           .c_str( ) );


    std::string sPointName;
    if( type_defs::ORTHOTOPE_DIMS == 0 )
        sPointName = "point";
    else if( type_defs::ORTHOTOPE_DIMS == 1 )
        sPointName = "interval";
    else if( type_defs::ORTHOTOPE_DIMS == 2 )
        sPointName = "rectangle";
    else if( type_defs::ORTHOTOPE_DIMS == 3 )
        sPointName = "cube";
    else
        sPointName = "" + std::to_string( type_defs::ORTHOTOPE_DIMS ) + "-orthotope";


    if constexpr( !type_defs::IS_ORTHOTOPE )
        xMain.def(
            "add_point",
            []( sps::Index<type_defs>& rM, typename type_defs::pos_t vPos, std::string sDesc ) {
                rM.addPoint( vPos, sDesc );
            },
            pybind11::arg( "pos" ),
            pybind11::arg( "desc" ) = "",
            ( R"pbdoc(
    Append a point to the data structure.
    
    :param pos: The position of the point.
    :type pos: list[int[)pbdoc" +
              std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]
    
    :param desc: A description for the Point, defaults to "".
    :type desc: str

    The point will not be queryable until generate is called.
)pbdoc" )
                .c_str( ) );
    else
        xMain.def(
            "add_point",
            []( sps::Index<type_defs>& rM, typename type_defs::ret_pos_t vStart, typename type_defs::ret_pos_t vEnd,
                std::string sDesc ) { rM.addPoint( vStart, vEnd, sDesc ); },
            pybind11::arg( "start" ),
            pybind11::arg( "end" ),
            pybind11::arg( "desc" ) = "",
            ( R"pbdoc(
    Append a )pbdoc" +
              sPointName + R"pbdoc( to the data structure.
    
    :param start: The bottom left position of the )pbdoc" +
              sPointName + R"pbdoc(.
    :type start: list[int[)pbdoc" +
              std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]
    
    :param end: The top right position of the )pbdoc" +
              sPointName + R"pbdoc(.
    :type end: list[int[)pbdoc" +
              std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]
    
    :param desc: A description for the Point, defaults to "".
    :type desc: str

    The )pbdoc" +
              sPointName + R"pbdoc( will not be queryable until generate is called.

    Dimensions 1 - )pbdoc" +
              std::to_string( type_defs::ORTHOTOPE_DIMS ) +
              R"pbdoc( of start and end may have different values, where the value of end must be larger equal to the value of start.
    Dimensions )pbdoc" +
              std::to_string( type_defs::ORTHOTOPE_DIMS + 1 ) + " - " +
              std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) +
              R"pbdoc( of start and end must have equal values.

    Note that this function will add one point for each outside corner of the given )pbdoc" +
              sPointName + "." )
                .c_str( ) );


    xMain
        .def( pybind11::init<std::string, bool>( ),
              pybind11::arg( "path" ) = "",
              pybind11::arg( "write_mode" ) = false,
              R"pbdoc(
    Create a new Index.
    
    :param path: Prefix path of the index on the filesystem (multiple files with different endings will be created), defaults to "".
    :type path: str
    
    :param write_mode: Open the index in write mode (if this is set to False no changes can be made to the index), defaults to False.
    :type write_mode: str
)pbdoc" ) // constructor
        .def( "generate",
              &sps::Index<type_defs>::generate,
              pybind11::arg( "from_points" ) = 0,
              pybind11::arg( "to_points" ) = std::numeric_limits<typename type_defs::coordinate_t>::max( ),
              pybind11::arg( "factor" ) = 0,
              pybind11::arg( "verbosity" ) = 1,
              R"pbdoc(
    Generate a new dataset.
    
    :param from_points: Size of the index before adding the first point of this dataset, defaults to zero.
    :type from_points: int
    
    :param to_points: Size of the index after adding all points of this dataset, defaults to the current size of the index.
    :type to_points: int
    
    :param verbosity: Degree of verbosity while creating the dataset, defaults to 1.
    :type verbosity: int

    :return: The id of the generated dataset.
    :rtype: int

    This may take a long time to compute.

    Use len(index) to determine the index of the first and last point, as add_point may add multiple points per call.

    This function is multithreaded.
)pbdoc" )
        .def( "count", &sps::Index<type_defs>::count, pybind11::arg( "dataset_id" ), pybind11::arg( "from_pos" ),
              pybind11::arg( "to_pos" ), //
              pybind11::arg( "intersection_type" ) = IntersectionType::enclosed, //
              pybind11::arg( "verbosity" ) = 0,
              ( R"pbdoc(
    Count the number of points between from and to in the given dataset.
    
    :param dataset_id: The id of the dataset to query
    :type dataset_id: int
    
    :param from_pos: The bottom left position of the query region.
    :type from_pos: list[int[)pbdoc" +
                std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]
    
    :param to_pos: The top right position of the query region.
    :type to_pos: list[int[)pbdoc" +
                std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]
    
    :param verbosity: Degree of verbosity while counting, defaults to 0.
    :type verbosity: int

    :return: The number of points in dataset_id between from_pos and to_pos.
    :rtype: int

    to_pos must be larger equal than from_pos in each dimension.
)pbdoc" )
                  .c_str( ) )
        .def( "count_multiple", &sps::Index<type_defs>::countMultiple, pybind11::arg( "dataset_id" ),
              pybind11::arg( "regions" ), //
              pybind11::arg( "intersection_type" ) = IntersectionType::enclosed, //
              pybind11::arg( "verbosity" ) = 0,
              ( R"pbdoc(
    Count the number of points between from and to in the given dataset.

    Counts for multiple regions.
    
    :param dataset_id: The id of the dataset to query
    :type dataset_id: int
    
    :param regions: The bottom left and top right positions of the queried regions. Given as a list (individual regions) of tuples (bottom-left, top-right) of lists (individual coordinates). The top right of each region must be larger equal the bottom left.
    :type regions: list[tuple[list[int[)pbdoc" +
                std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]], list[int[)pbdoc" +
                std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]]]
    
    :param verbosity: Degree of verbosity while counting, defaults to 0.
    :type verbosity: int

    :return: The number of points in dataset_id between from_pos and to_pos.
    :rtype: int

    to_pos must be larger equal than from_pos in each dimension.
)pbdoc" )
                  .c_str( ) )
        .def( "count_size_limited", &sps::Index<type_defs>::countSizeLimited, pybind11::arg( "dataset_id" ),
              pybind11::arg( "from_pos" ), pybind11::arg( "to_pos" ), //
              pybind11::arg( "intersection_type" ) = IntersectionType::enclosed, //
              pybind11::arg( "verbosity" ) = 0,
              ( R"pbdoc(
    Count the number of points between from and to and in the given dataset.

    As opposed to count, this function allows specifying the start and end positions for all dimensions in the Datastructure. 
    This is only relevant for indices with orthotope dimensions.
    
    :param dataset_id: The id of the dataset to query
    :type dataset_id: int
    
    :param from_pos: The bottom left position of the query region.
    :type from_pos: list[int[)pbdoc" +
                std::to_string( type_defs::D ) + R"pbdoc(]]
    
    :param to_pos: The top right position of the query region.
    :type to_pos: list[int[)pbdoc" +
                std::to_string( type_defs::D ) + R"pbdoc(]]
    
    :param verbosity: Degree of verbosity while counting, defaults to 0.
    :type verbosity: int

    :return: The number of points in dataset_id between from_pos and to_pos.
    :rtype: int

    to_pos must be larger equal than from_pos in each dimension.
)pbdoc" )
                  .c_str( ) )
        //.def( "get", &sps::Index<type_defs>::get, "" )
        .def( "__str__", &sps::Index<type_defs>::str,
              "Return a string describing the index. Very slow for large datasets." )
        .def( "max_prefix_value", &sps::Index<type_defs>::maxPrefixSumValue,
              "Return the maximal stored prefix sum. Intended for storage space optimization purposes." )
        .def( "__len__", &sps::Index<type_defs>::numPoints, "Total number of points (among all datasets)." )
        .def( "clear", &sps::Index<type_defs>::clear, "Clear the complete index." )
        .def( "__get_overlay_info", &sps::Index<type_defs>::getOverlayInfo )
        .def( "get_overlay_grid", &sps::Index<type_defs>::getOverlayGrid, pybind11::arg( "dataset_id" ),
              "Returns the bottom-left-front-... and top-right-back-... position of all overlays." )
        .def_static( "estimate_num_elements", &sps::Index<type_defs>::estimateDataStructureElements )
        .def_static( "estimate_size", &sps::Index<type_defs>::estimateDataStructureSize )
        .def_static( "pick_num_overlays", &sps::Index<type_defs>::pickNumOverlays )
        .def_static( "to_factor", &sps::Index<type_defs>::toFactor )
        .def( "get_num_internal_prefix_sums", &sps::Index<type_defs>::getNumInternalPrefixSums )
        .def( "get_num_overlay_prefix_sums", &sps::Index<type_defs>::getNumOverlayPrefixSums )
        .def( "get_num_internal_sparse_coords", &sps::Index<type_defs>::getNumInternalSparseCoords )
        .def( "get_num_overlay_sparse_coords", &sps::Index<type_defs>::getNumOverlaySparseCoords )

        ;
    return "    " + sName + "\n";
}
#endif