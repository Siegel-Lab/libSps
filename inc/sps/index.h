/**
 * @file index.h
 * @brief Implements the main Datastructure.
 * @author Markus Schmidt
 */


#pragma once

#include "sps/abstract_index.h"
#include "sps/corners.h"
#include "sps/dataset.h"
#include "sps/desc.h"
#include "sps/nd_grid.h"
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
    using corner_t = Corner<type_defs>;

    using corners_t = Corners<type_defs>;

    using overlay_t = AlignedPower2<Overlay<type_defs>>;
    using sparse_coord_t = SparseCoord<type_defs>;
    using overlay_grid_t = typename Overlay<type_defs>::overlay_grid_t;
    using prefix_sum_grid_t = typename Overlay<type_defs>::prefix_sum_grid_t;
    using dataset_t = AlignedPower2<Dataset<type_defs>>;

    EXTRACT_VEC_GENERATOR( dataset, dataset_t ); // macro call

    corners_t vCorners;
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
        : vSparseCoord( sPrefix, bWrite ),
          vPrefixSumGrid( sPrefix + ".prefix_sums", bWrite ),
          vOverlayGrid( sPrefix + ".overlays", bWrite ),
          xFile( dataset_vec_generator.file( sPrefix + ".datsets", bWrite ) ),
          vDataSets( dataset_vec_generator.vec( xFile ) )
    {
        // corners should never be persistent...
        clearCorners( );
    }

    /**
     * @brief Clear the complete index.
     *
     * Clears all datasets and all points.
     */
    void clear( )
    {
        clearCorners( );
        clearKeepPoints( );
    }

    void reserve( size_t uiCorners, size_t uiCoords, size_t uiPrefixSums, size_t uiOverlays, size_t uiDatasets )
    {
        vCorners.vData.reserve( uiCorners );
        vSparseCoord.vData.reserve( uiCoords );
        vPrefixSumGrid.vData.reserve( uiPrefixSums );
        vOverlayGrid.vData.reserve( uiOverlays );
        vDataSets.reserve( uiDatasets );
    }

    /**
     * @brief Clear the complete index.
     *
     * @details
     * Clears all datasets and all points.
     */
    void clearCorners( )
    {
        vCorners.clear( );
    }

    /**
     * @brief Clear everything but the points.
     */
    void clearKeepPoints( )
    {
        vSparseCoord.clear( );
        vPrefixSumGrid.clear( );
        vOverlayGrid.clear( );
        vDataSets.clear( );
    }

    void popDataset( )
    {
        size_t uiRmSp =
            getNumInternalSparseCoords( vDataSets.size( ) - 1 ) + getNumOverlaySparseCoords( vDataSets.size( ) - 1 );
        vSparseCoord.vData.resize( vSparseCoord.vData.size( ) >= uiRmSp ? vSparseCoord.vData.size( ) - uiRmSp : 0 );

        size_t uiRmPs =
            getNumInternalPrefixSums( vDataSets.size( ) - 1 ) + getNumOverlayPrefixSums( vDataSets.size( ) - 1 );
        vPrefixSumGrid.vData.resize( vPrefixSumGrid.vData.size( ) >= uiRmPs ? vPrefixSumGrid.vData.size( ) - uiRmPs
                                                                            : 0 );

        size_t uiRmOG = getNumOverlays( vDataSets.size( ) - 1 );
        vOverlayGrid.vData.resize( vOverlayGrid.vData.size( ) >= uiRmOG ? vOverlayGrid.vData.size( ) - uiRmOG : 0 );
        vDataSets.pop_back( );
    }

    /**
     * @brief Append a point to the data structure.
     *
     * @details
     * The point will not be queryable until generate is called.
     *
     * @tparam trigger This function is only active if IS_ORTHOTOPE = false
     * @param vPos The position of the point.
     * @param uiVal The value of the point, defaults to 1.
     */
    template <bool trigger = !IS_ORTHOTOPE>
    typename std::enable_if_t<trigger> addPoint( ret_pos_t vPos, val_t uiVal = 1 )
    {
        vCorners.add( vPos, uiVal );
    }

  private:
    coordinate_t addDimFrom( coordinate_t uiSize, IntersectionType xIntersectionType ) const
    {
        if( xIntersectionType == IntersectionType::slice )
            return uiSize;
        else if( xIntersectionType == IntersectionType::encloses )
            return 1 + uiSize;
        else
            return 0;
    }

    coordinate_t addDimTo( coordinate_t uiSize, IntersectionType xIntersectionType ) const
    {

        if( xIntersectionType == IntersectionType::slice || xIntersectionType == IntersectionType::enclosed )
            return uiSize;
        else if( xIntersectionType == IntersectionType::points_only )
            return 1;
        else
            return std::numeric_limits<coordinate_t>::max( );
    }

    std::array<pos_t, 2> addDims( ret_pos_t vStart, ret_pos_t vEnd, isect_arr_t vIntersectionTypes ) const
    {
        std::array<pos_t, 2> vRet;
        for( size_t uiI = 0; uiI < D - ORTHOTOPE_DIMS; uiI++ )
        {
            vRet[ 0 ][ uiI ] = vStart[ uiI ];
            vRet[ 1 ][ uiI ] = vEnd[ uiI ];
        }
        for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
        {
            vRet[ 0 ][ uiI + D - ORTHOTOPE_DIMS ] =
                addDimFrom( vEnd[ uiI ] - vStart[ uiI ], vIntersectionTypes[ uiI ] );
            vRet[ 1 ][ uiI + D - ORTHOTOPE_DIMS ] = addDimTo( vEnd[ uiI ] - vStart[ uiI ], vIntersectionTypes[ uiI ] );
        }
        return vRet;
    }


    std::array<pos_t, 2> addDims( ret_pos_t vStart, ret_pos_t vEnd, 
                                  IntersectionType xInterType = IntersectionType::slice ) const
    {
        isect_arr_t vInterTypes;
        for( size_t uiI = 0; uiI < D - ORTHOTOPE_DIMS; uiI++ )
            vInterTypes[ uiI ] = xInterType;

        return addDims(vStart, vEnd, vInterTypes);
    }

    std::array<std::array<coordinate_t, ORTHOTOPE_DIMS>, 2>
    addDims( ret_pos_t vSize, IntersectionType xIntersectionType = IntersectionType::slice ) const
    {
        std::array<std::array<coordinate_t, ORTHOTOPE_DIMS>, 2> vRet;
        for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
        {
            vRet[ 0 ][ uiI ] = addDimFrom( vSize[ uiI ], xIntersectionType );
            vRet[ 1 ][ uiI ] = addDimTo( vSize[ uiI ], xIntersectionType );
        }
        return vRet;
    }

    typename corners_t::Entry makeEntry( )
    {
        typename corners_t::Entry xCorners;
        xCorners.uiStartIndex = 0;
        xCorners.uiEndIndex = vCorners.size( );
        return xCorners;
    }

  public:
    /**
     * @brief Append an orthotope to the data structure.
     *
     * The orthotope will not be queryable until generate is called.
     *
     * @tparam trigger This function is only active if IS_ORTHOTOPE = true
     * @param vStart The bottom left position of the orthotope.
     * @param vEnd The top right position of the orthotope.
     * @param uiVal The value of the orthotope, defaults to 1.
     */
    template <bool trigger = IS_ORTHOTOPE>
    typename std::enable_if_t<trigger> addPoint( ret_pos_t vStart, ret_pos_t vEnd, val_t uiVal = 1 )
    {
        for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
            if( vStart[ uiI ] > vEnd[ uiI ] )
                throw std::invalid_argument( "end must be larger-equal than start for orthotope dimensions." );
        for( size_t uiI = ORTHOTOPE_DIMS; uiI < D - ORTHOTOPE_DIMS; uiI++ )
            if( vStart[ uiI ] != vEnd[ uiI ] )
                throw std::invalid_argument( "end must equal start for non-orthotope dimensions." );
        auto vP = addDims( vStart, vEnd );
        vCorners.add( vP[ 0 ], vP[ 1 ], uiVal );
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
     * @param fFac Overlay size factor, defaults to -1.
     *             -1: factor is picked so that the estimated size of the datastructure is minimal
     * @param uiVerbosity Degree of verbosity while creating the dataset, defaults to 1.
     * @param uiNumOverlaySamples number of overlays to sample, default to 10000.
     * @param uiNumPointSamples number of points to sample per overlay, default to 10000.
     * @return class_key_t The id of the generated dataset.
     */
    class_key_t generate( double fFac = -1, size_t uiVerbosity = 1 )
    {

        progress_stream_t xProg( uiVerbosity );
        // generate the dataset in ram then push it into the index to make sure that the cache of the vector
        // does not unload the memory half way through the initialization. (not relevant for std::vector
        // implementations)
        dataset_t xNew( vOverlayGrid, vSparseCoord, vPrefixSumGrid, vCorners, makeEntry( ), fFac, xProg );
        class_key_t uiRet = vDataSets.size( );
        vDataSets.push_back( xNew );

        vCorners.clear( );

#ifndef NDEBUG
        xProg << Verbosity( 1 ) << "\n\nMaximal prefix sum value: " << maxPrefixSumValue( ) << ".\n";
#endif
        return uiRet;
    }

  private:
    template <size_t uiD, size_t uiDistToTo, size_t> struct SizeLimitedInvariant
    {
        static void count( pos_t vPos, const isect_arr_t& vInterTypes, const dataset_vec_t& vDataSets,
                           const sparse_coord_t& vSparseCoord, const prefix_sum_grid_t& vPrefixSumGrid,
                           const overlay_grid_t& vOverlayGrid, const class_key_t xDatasetId, 
                           const size_t uiIntersectTypeFac, val_t& uiRet
#if GET_PROG_PRINTS
                           ,
                           progress_stream_t& xProg
#endif
        )
        {
#if GET_PROG_PRINTS
            xProg << "query: " << xDatasetId << " vPos: " << vPos << " uiD: " << uiD << "\n";
#endif

            val_t uiCurr =
                vDataSets[ xDatasetId ].get( vOverlayGrid, vSparseCoord, vPrefixSumGrid, vPos,
                                             Overlay<type_defs>::template intersectionTypeToCornerIndex<uiD>( vInterTypes )
#if GET_PROG_PRINTS
                                                 ,
                                             xProg
#endif
                );

#ifndef NDEBUG
#if DU_UNREALISTIC_VALUE_CHECK
            if( uiCurr >= std::numeric_limits<val_t>::max( ) / 2 )
                throw std::runtime_error( "unrealistic value for uiCurr" );
#endif
#endif

            val_t uiFac = ( uiDistToTo % 2 == 0 ? 1 : -1 ) * uiIntersectTypeFac;
#if GET_PROG_PRINTS
            xProg << "is " << ( uiFac == 1 ? "+" : "-" ) << uiCurr << " [" << uiD << "/"
                  << ( 1 << ( D - ORTHOTOPE_DIMS ) ) << "]"
                  << "\n";
#endif
            uiRet += uiCurr * uiFac;
        }
    };


    static inline bool countSizeLimitedInvariantCond( coordinate_t uiPos, size_t /*uiD*/, bool /*bIsFrom*/ )
    {
        return uiPos != std::numeric_limits<coordinate_t>::max( );
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
     * @return val_t The number of points in dataset_id between from_pos and to_pos.
     */
    val_t countSizeLimited( class_key_t xDatasetId, pos_t vFrom, pos_t vTo, isect_arr_t vInterTypes,
                            size_t
#if GET_PROG_PRINTS
                                uiVerbosity
#endif
                            = 0 ) const
    {
        if( vDataSets[ xDatasetId ].getNumOverlays( ) == 0 )
            return val_t{ };
#if GET_PROG_PRINTS
        progress_stream_t xProg( uiVerbosity );
        xProg << "countSizeLimited " << vFrom << " to " << vTo << "\n";
#endif
        val_t uiRet = 0;
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            --vFrom[ uiI ];
            --vTo[ uiI ];
        }
        forAllCombinationsTmpl<pos_t, SizeLimitedInvariant>( vFrom, vTo, countSizeLimitedInvariantCond, vInterTypes,
                                                             vDataSets, vSparseCoord, vPrefixSumGrid, vOverlayGrid,
                                                             xDatasetId, 
                                                             Overlay<type_defs>::intersectionTypeToFactor(vInterTypes), 
                                                             uiRet
#if GET_PROG_PRINTS
                                                             ,
                                                             xProg
#endif
        );

#if GET_PROG_PRINTS
        xProg << "countSizeLimited uiRet=" << uiRet << "\n";
#if DU_UNREALISTIC_VALUE_CHECK
        if( uiRet >= std::numeric_limits<val_t>::max( ) / 2 )
            throw std::runtime_error( "unrealistic value for uiRet" );
#endif
#endif

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
    val_t count( class_key_t xDatasetId, ret_pos_t vFromR, ret_pos_t vToR, isect_arr_t vInterTypes,
                 size_t uiVerbosity = 0 ) const
    {
        for( size_t uiI = 0; uiI < D - ORTHOTOPE_DIMS; uiI++ )
            if( vFromR[ uiI ] > vToR[ uiI ] )
                throw std::invalid_argument( "end must be larger-equal than start." );

        auto vP = addDims( vFromR, vToR, vInterTypes );
        pos_t vFrom = vP[ 0 ];
        pos_t vTo = vP[ 1 ];
        return countSizeLimited( xDatasetId, vFrom, vTo, vInterTypes, uiVerbosity );
    }

    val_t count( class_key_t xDatasetId, ret_pos_t vFromR, ret_pos_t vToR,
                 IntersectionType xInterType = IntersectionType::enclosed, size_t uiVerbosity = 0 ) const
    {
        isect_arr_t vInterTypes;
        for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
            vInterTypes[ uiI ] = xInterType;

        return count( xDatasetId, vFromR, vToR, vInterTypes, uiVerbosity );
    }

    /**
     * @brief Count the number of points between from and to and in the given dataset.
     *
     * to_pos must be larger equal than from_pos in each dimension.
     *
     * @todo this interface needs to change:
     * instead of vPos, vSize, and vNum, there should be for each dimension a vector of the coordinates that need to be
     * queried.
     * The the resulting vector returns the values of the cells in between
     *
     * @param xDatasetId The id of the dataset to query
     * @param vGrid The coordinates of the grid-lines for each dimensions.
     * @param xInterType The used intersection type, defaults to enclosed. Ignored if there are no orthotope dimensions.
     * @param uiVerbosity Degree of verbosity while counting, defaults to 0.
     * @return val_t the values of the grid cells.
     */
    std::vector<val_t> gridCount( class_key_t xDatasetId, [[maybe_unused]] std::array<std::vector<coordinate_t>, D> vGrid,
                                  [[maybe_unused]] IntersectionType xInterType = IntersectionType::enclosed,
                                  [[maybe_unused]] size_t uiVerbosity = 0 ) const
    {
        if( vDataSets[ xDatasetId ].getNumOverlays( ) == 0 )
            return std::vector<val_t>{ };
        return std::vector<val_t>{ };
#if 0
        auto vP = addDims( vSize, xInterType );
        std::array<coordinate_t, ORTHOTOPE_DIMS> vOrthoFrom = vP[ 0 ];
        std::array<coordinate_t, ORTHOTOPE_DIMS> vOrthoTo = vP[ 1 ];
#if GET_PROG_PRINTS
        progress_stream_t xProg( uiVerbosity );
        xProg << "gridCount pos " << vPos << " size " << vSize << " amount " << vNum << " vOrthoFrom " << vOrthoFrom << " vOrthoTo " << vOrthoTo << "\n";
#endif
        for( size_t uiI = 0; uiI < ORTHOTOPE_DIMS; uiI++ )
        {
            --vOrthoFrom[ uiI ];
            --vOrthoTo[ uiI ];
        }
        for( size_t uiI = 0; uiI < D - ORTHOTOPE_DIMS; uiI++ )
            --vPos[ uiI ];

        return vDataSets[ xDatasetId ].grid( vOverlayGrid, vSparseCoord, vPrefixSumGrid, vPos, vSize, vNum, vOrthoFrom,
                                             vOrthoTo, xInterType
#if GET_PROG_PRINTS
                                             ,
                                             xProg
#endif
        );
#endif
    }


  private:
#define COUNT_MULTIPLE( args... )                                                                                      \
    std::vector<val_t> vRet;                                                                                           \
    vRet.reserve( vRegions.size( ) );                                                                                  \
                                                                                                                       \
    for( size_t uiI = 0; uiI < vRegions.size( ); uiI++ )                                                               \
    {                                                                                                                  \
        if( PyErr_CheckSignals( ) != 0 )                                                                               \
            throw pybind11::error_already_set( );                                                                      \
        vRet.push_back( count( args ) );                                                                               \
    }                                                                                                                  \
    return vRet;

  public:
    /**
     * @brief Count the number of points between from and to and in the given dataset.
     *
     * Counts for multiple regions.
     * to_pos must be larger equal than from_pos in each dimension.
     *
     * @param vRegions Id of the dataset to query, The bottom left and top right positions of the queried regions.
     * @param xInterType The used intersection type, defaults to enclosed. Ignored if there are no orthotope dimensions.
     * @param uiVerbosity Degree of verbosity while counting, defaults to 0.
     * @return std::vector<val_t> The number of points in dataset_id between from_pos and to_pos for each given region.
     */
    std::vector<val_t> countMultiple( std::vector<std::tuple<class_key_t, ret_pos_t, ret_pos_t>> vRegions,
                                      IntersectionType xInterType = IntersectionType::enclosed,
                                      size_t uiVerbosity = 0 ) const
    {
        COUNT_MULTIPLE( std::get<0>( vRegions[ uiI ] ), std::get<1>( vRegions[ uiI ] ), std::get<2>( vRegions[ uiI ] ),
                        xInterType, uiVerbosity );
    }

    std::vector<val_t> countMultiple( class_key_t uiDatasetId, std::vector<std::tuple<ret_pos_t, ret_pos_t>> vRegions,
                                      IntersectionType xInterType = IntersectionType::enclosed,
                                      size_t uiVerbosity = 0 ) const
    {
        COUNT_MULTIPLE( uiDatasetId, std::get<0>( vRegions[ uiI ] ), std::get<1>( vRegions[ uiI ] ), xInterType,
                        uiVerbosity );
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

    /**
     * @brief Count the number of internal prefix sums stored in the dataset with id dataset_id.
     *
     * @param xDatasetId The id of the dataset to query
     * @return coordinate_t number of internal prefix sums
     */
    coordinate_t getNumInternalPrefixSums( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getNumInternalPrefixSums( vOverlayGrid, vSparseCoord );
    }

    /**
     * @brief Count the number of overlay prefix sums stored in the dataset with id dataset_id.
     *
     * @param xDatasetId The id of the dataset to query
     * @return coordinate_t number of overlay prefix sums
     */
    coordinate_t getNumOverlayPrefixSums( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getNumOverlayPrefixSums( vOverlayGrid, vSparseCoord );
    }

    /**
     * @brief Count the number of internal sparse coordinates stored in the dataset with id dataset_id.
     *
     * @param xDatasetId The id of the dataset to query
     * @return coordinate_t number of internal sparse coordinates
     */
    coordinate_t getNumInternalSparseCoords( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getNumInternalSparseCoords( vOverlayGrid );
    }

    /**
     * @brief Count the number of overlay sparse coordinates stored in the dataset with id dataset_id.
     *
     * @param xDatasetId The id of the dataset to query
     * @return coordinate_t number of overlay sparse coordinates
     */
    coordinate_t getNumOverlaySparseCoords( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getNumOverlaySparseCoords( vOverlayGrid );
    }

    /**
     * @brief Count the number of overlay blocks.
     *
     * @param xDatasetId The id of the dataset to query
     * @return coordinate_t number of overlay blocks
     */
    coordinate_t getNumOverlays( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getNumOverlays( );
    }

    /**
     * @brief Returns the size of the dataset in bytes.
     *
     * @param xDatasetId id of the queried dataset
     * @return coordinate_t size in bytes
     */
    coordinate_t getSize( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getSize( vOverlayGrid, vSparseCoord );
    }

    /**
     * @brief Predict the number of data structure elements stored in a dataset.
     *
     * @details
     *
     * Here f is proportional to the number of boxes used in the data structure.
     * For any dataset, there is an optimal value for f that leads to the minimal data structure size.
     * Since the data structure size is proportional to the time required to build the datastructure this also minimizes
     * construction time.
     *
     * The purpose of this function is to find the optimal value for f.
     *
     * Uses a statistical approach to predict the number of elements.
     * For details see the corresponding github page or our manuscript.
     *
     * There are five values predicted:
     * - The number of internal prefix sums
     * - The number of overlay prefix sums
     * - The number of internal sparse coordinates
     * - The number of overlay sparse coordinates
     * - The number of global sparse coordinates
     * - The total size in bytes
     *
     *
     * @param vFac list of factors that are proportional to the number of boxes within the data structure
     * @param uiNumOverlaySamples number of overlays to sample, default to 10000.
     * @param uiNumPointSamples number of points to sample per overlay, default to 10000.
     * @return std:vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>>
     *         The predicted number of dataset structure elements for each factor
     */
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>>
    estimateDataStructureElements( std::vector<double> vFac )
    {
        return dataset_t::estimateDataStructureElements( vCorners, makeEntry( ), vFac );
    }


    /**
     * @brief Get the axis sizes of the inserted points.
     *
     * @return pos_t an array containing the axis size for each dimension.
     */
    pos_t gridSize( )
    {
        return dataset_t::generateCoordSizes( vCorners, makeEntry( ) )[ 0 ];
    }

    /**
     * @brief Predict the best f.
     *
     * @details
     *
     * Here f is a factor proportional to the number of boxes used in the data structure.
     * See estimate_num_elements for a detailed description.
     *
     * @param uiNumOverlaySamples number of overlays to sample, default to 10000.
     * @param uiNumPointSamples number of points to sample per overlay, default to 10000.
     * @param uiVerbosity Degree of verbosity while creating the dataset, defaults to 1.
     * @return uint64_t The predicted best value for f
     */
    uint64_t pickNumOverlays( size_t uiVerbosity = 0 )
    {
        progress_stream_t xProg( uiVerbosity );
        return dataset_t::pickNumOverlays( vCorners, makeEntry( ), xProg );
    }

    /**
     * @brief Convert a given number of overlay blocks to the factor f.
     *
     * @details
     *
     * Here f is a factor proportional to the number of boxes used in the data structure.
     * See estimate_num_elements for a detailed description.
     *
     * @param uiNumOverlays number of overlay blocks
     * @return double f
     */
    double toFactor( uint64_t uiNumOverlays )
    {
        return dataset_t::toFactor( vCorners, makeEntry( ), uiNumOverlays );
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
        os << "vCorners: ";
        os << rMain.vCorners << std::endl;

        os << "Pretty Print: ";
        for( size_t uiI = 0; uiI < rMain.vDataSets.size( ); uiI++ )
            rMain.vDataSets[ uiI ].stream( os, rMain.vOverlayGrid, rMain.vSparseCoord, rMain.vPrefixSumGrid,
                                           rMain.vCorners )
                << std::endl;

        return os;
    }

    /**
     * @brief Return information about the overlay boxes.
     *
     * @param xDatasetId The id of the dataset to query
     * @return std::vector<typename dataset_t::OverlayInfo> information about all overlays
     */
    std::vector<typename dataset_t::OverlayInfo> getOverlayInfo( class_key_t xDatasetId ) const
    {
        return vDataSets[ xDatasetId ].getOverlayInfo( vOverlayGrid, vSparseCoord, vCorners );
    }

    /**
     * @brief Returns the bottom-left-front-... and top-right-back-... position of all overlays.
     */
    std::vector<std::array<pos_t, 3>> getOverlayGrid( class_key_t xDatasetId ) const
    {
        if( vDataSets.size( ) <= xDatasetId )
            return std::vector<std::array<pos_t, 3>>{ };
        return vDataSets[ xDatasetId ].getOverlayGrid( vOverlayGrid );
    }

    val_t getMaxPrefixSum( ) const
    {
        val_t uiMax = 0;
        for( const sps_t& rSps : vPrefixSumGrid.vData )
        {
            if constexpr( IS_ORTHOTOPE )
                for( size_t uiX = 0; uiX < 1 << ORTHOTOPE_DIMS; uiX++ )
                    uiMax = std::max( uiMax, rSps[ uiX ] );
            else
                uiMax = std::max( uiMax, rSps );
        }
        return uiMax;
    }
};

} // namespace sps


#if WITH_PYTHON
void exportEnum( pybind11::module& m )
{
    pybind11::enum_<sps::IntersectionType>( m, "IntersectionType", R"pbdoc(
    An Enum for Querying the index.
    
    Which orthotopes to count, depending on how they intersect the queried area.
    Only relevant for the Index.count() function.
)pbdoc" )
        .value( "enclosed", sps::IntersectionType::enclosed,
                "count orthotopes that are fully enclosed by the queried area" )
        .value( "encloses", sps::IntersectionType::encloses, "count orthotopes that fully enclose by the queried area" )
        .value( "overlaps", sps::IntersectionType::overlaps, "count orthotopes that overlap the queried area" )
        .value( "first", sps::IntersectionType::first,
                "count orthotopes that have their bottom-left-front-.. corner in the queried area" )
        .value( "last", sps::IntersectionType::last,
                "count orthotopes that have their top-right-back-.. corner in the queried area" )
        .value( "points_only", sps::IntersectionType::points_only,
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
        .def_readonly( "points", &OI::vvCorners )

        ;


    pybind11::class_<sps::Index<type_defs>, sps::AbstractIndex> xMain( m, sName.c_str( ),
                                                                       ( R"pbdoc(
    Prefix sum index )pbdoc" + sDesc + R"pbdoc(.
    
    .. automethod:: add_point
    .. automethod:: generate
    .. automethod:: count
    .. automethod:: count_size_limited
    .. automethod:: count_multiple
    .. automethod:: __str__
    .. automethod:: clear
    .. automethod:: clear_keep_points
    .. automethod:: get_overlay_grid
    .. automethod:: estimate_num_elements
    .. automethod:: pick_num_overlays
    .. automethod:: to_factor
    .. automethod:: get_num_internal_prefix_sums
    .. automethod:: get_num_overlay_prefix_sums
    .. automethod:: get_num_internal_sparse_coords
    .. automethod:: get_num_overlay_sparse_coords
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
            []( sps::Index<type_defs>& rM, typename type_defs::pos_t vPos, typename type_defs::val_t uiVal ) {
                rM.addPoint( vPos, uiVal );
            },
            pybind11::arg( "pos" ),
            pybind11::arg( "val" ) = 1,
            ( R"pbdoc(
    Append a point to the data structure.
    
    :param pos: The position of the point.
    :type pos: list[int[)pbdoc" +
              std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]
    
    :param val: The value of the point.
    :type val: int
    
    :param desc: A description for the Point, defaults to "".
    :type desc: str

    The point will not be queryable until generate is called.
)pbdoc" )
                .c_str( ) );
    else
        xMain.def(
            "add_point",
            []( sps::Index<type_defs>& rM, typename type_defs::ret_pos_t vStart, typename type_defs::ret_pos_t vEnd,
                typename type_defs::val_t uiVal ) { rM.addPoint( vStart, vEnd, uiVal ); },
            pybind11::arg( "start" ),
            pybind11::arg( "end" ),
            pybind11::arg( "val" ) = 1,
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

    :param val: The value of the )pbdoc" +
              sPointName + R"pbdoc(.
    :type val: int
    
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
              pybind11::arg( "factor" ) = -1,
              pybind11::arg( "verbosity" ) = 1,
              R"pbdoc(
    Generate a new dataset from the previously added points.

    :param verbosity: Degree of verbosity while creating the dataset, defaults to 1.
    :type verbosity: int

    :return: The id of the generated dataset.
    :rtype: int

    This may take a long time to compute.

    Use len(index) to determine the index of the first and last point, as add_point may add multiple points per call.

    This function is multithreaded.
)pbdoc" )
        .def(
            "count",
            []( const sps::Index<type_defs>& rM, typename type_defs::class_key_t xDatasetId,
                typename type_defs::ret_pos_t vFromR, typename type_defs::ret_pos_t vToR,
                sps::IntersectionType xInterType,
                size_t uiVerbosity ) { return rM.count( xDatasetId, vFromR, vToR, xInterType, uiVerbosity ); },
            pybind11::arg( "dataset_id" ), pybind11::arg( "from_pos" ),
            pybind11::arg( "to_pos" ), //
            pybind11::arg( "intersection_type" ) = sps::IntersectionType::enclosed, //
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
    
    :param intersection_type: Which data elements to count (enclosed, overlapping, enclosing, etc...).
    :type intersection_type: IntersectionType
    
    :param verbosity: Degree of verbosity while counting, defaults to 0.
    :type verbosity: int

    :return: The number of points in dataset_id between from_pos and to_pos.
    :rtype: int

    to_pos must be larger equal than from_pos in each dimension.
)pbdoc" )
                .c_str( ) )
        .def(
            "count",
            []( const sps::Index<type_defs>& rM, typename type_defs::class_key_t xDatasetId,
                typename type_defs::ret_pos_t vFromR, typename type_defs::ret_pos_t vToR,
                typename type_defs::isect_arr_t vInterTypes,
                size_t uiVerbosity ) { return rM.count( xDatasetId, vFromR, vToR, vInterTypes, uiVerbosity ); },
            pybind11::arg( "dataset_id" ), pybind11::arg( "from_pos" ),
            pybind11::arg( "to_pos" ), //
            pybind11::arg( "intersection_types" ), //
            pybind11::arg( "verbosity" ) = 0,
            ( R"pbdoc(
    Count the number of points between from and to in the given dataset.
    This overload allows you to specify, for each dimension, how to deal with data-hyperrectangles.
    E.g. you can specify to count data-hyperrectangles that overlap the query-hyperrectangle in dimension 1
    but are fully enclosed by the query-hyperrectangle in dimension 2,
    by providing intersection_types = [overlaps, encloses].
    
    :param dataset_id: The id of the dataset to query
    :type dataset_id: int
    
    :param from_pos: The bottom left position of the query region.
    :type from_pos: list[int[)pbdoc" +
              std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]
    
    :param to_pos: The top right position of the query region.
    :type to_pos: list[int[)pbdoc" +
              std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]
    
    :param intersection_types: Which data elements to count (enclosed, overlapping, enclosing, etc...), separated for each dimension.
    :type intersection_types: IntersectionType

    :param verbosity: Degree of verbosity while counting, defaults to 0.
    :type verbosity: int

    :return: The number of points in dataset_id between from_pos and to_pos.
    :rtype: int

    to_pos must be larger equal than from_pos in each dimension.
)pbdoc" )
                .c_str( ) )
        .def( "grid_count", &sps::Index<type_defs>::gridCount, pybind11::arg( "dataset_id" ), pybind11::arg( "grid" ),
              pybind11::arg( "intersection_type" ) = sps::IntersectionType::enclosed, //
              pybind11::arg( "verbosity" ) = 0,
              ( R"pbdoc(
    Count the number of points between from and to in the given dataset.
    
    :param dataset_id: The id of the dataset to query
    :type dataset_id: int
    
    :param pos: The bottom left position of the grid query region.
    :type pos: list[int[)pbdoc" +
                std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]

    :param size: The size of one cell in the grid.
    :type size: list[int[)pbdoc" +
                std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]
    
    :param num: The number of grid cells per dimension.
    :type num: list[int[)pbdoc" +
                std::to_string( type_defs::D - type_defs::ORTHOTOPE_DIMS ) + R"pbdoc(]]
    
    :param verbosity: Degree of verbosity while counting, defaults to 0.
    :type verbosity: int

    :return: The number of points in dataset_id between from_pos and to_pos.
    :rtype: int

    to_pos must be larger equal than from_pos in each dimension.
)pbdoc" )
                  .c_str( ) )
        .def(
            "count_multiple",
            []( const sps::Index<type_defs>& rM,
                std::vector<std::tuple<typename type_defs::class_key_t, typename type_defs::ret_pos_t,
                                       typename type_defs::ret_pos_t>>
                    vRegions,
                sps::IntersectionType xInterType,
                size_t uiVerbosity ) { return rM.countMultiple( vRegions, xInterType, uiVerbosity ); },
            pybind11::arg( "regions" ), //
            pybind11::arg( "intersection_type" ) = sps::IntersectionType::enclosed, //
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
              pybind11::arg( "intersection_type" ) = sps::IntersectionType::enclosed, //
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
        .def( "clear", &sps::Index<type_defs>::clear, "Clear the complete index." )
        .def( "reserve", &sps::Index<type_defs>::reserve,
              "Reserve memory for the index corners, sparse coords and prefix sums. Does nothing a value less than the "
              "currently held memory is given." )
        .def( "clear_keep_points", &sps::Index<type_defs>::clearKeepPoints,
              "Clear the index datasets, but keep the points." )
        .def( "pop_dataset", &sps::Index<type_defs>::popDataset, "Remove the last dataset." )
        .def( "__get_overlay_info", &sps::Index<type_defs>::getOverlayInfo,
              "Return information about the overlay boxes." )
        .def( "get_overlay_grid", &sps::Index<type_defs>::getOverlayGrid, pybind11::arg( "dataset_id" ),
              "Returns the bottom-left-front-... and top-right-back-... position of all overlays." )
        .def( "get_size", &sps::Index<type_defs>::getSize, pybind11::arg( "dataset_id" ),
              "Returns the size of the dataset in bytes." )
        .def( "get_max_prefix_sum", &sps::Index<type_defs>::getMaxPrefixSum, "Get the largest stored prefix sum." )
        .def( "estimate_num_elements", &sps::Index<type_defs>::estimateDataStructureElements, pybind11::arg( "f" ),
              R"pbdoc(
    Predict the number of data structure elements stored in a dataset generated from the currently added points.

    Here f is proportional to the number of boxes used in the data structure.
    For any dataset, there is an optimal value for f that leads to the minimal data structure size.
    Since the data structure size is proportional to the time required to build the datastructure this also minimizes construction time.

    The purpose of this function is to find the optimal value for f.

    Uses a statistical approach to predict the number of elements.
    For details see the corresponding github page or our manuscript.

    There are five values predicted:
    - The number of internal prefix sums
    - The number of overlay prefix sums
    - The number of internal sparse coordinates
    - The number of overlay sparse coordinates
    - The number of global sparse coordinates
    - Total size of the datastructure

    :param f: list of factors that are proportional to the number of boxes within the data structure
    :type f: list[float]

    :return: The predicted number of dataset structure elements for each factor
    :rtype: list[tuple[int]]
)pbdoc" )
        .def( "pick_num_overlays", &sps::Index<type_defs>::pickNumOverlays, pybind11::arg( "verbosity" ) = 0,
              R"pbdoc(
    Predict the best factor f for the currently added points.

    Here f is a factor proportional to the number of boxes used in the data structure.
    See estimate_num_elements for a detailed description. 

    :return: The predicted best value for f.
    :rtype: float
)pbdoc" )
        .def( "to_factor", &sps::Index<type_defs>::toFactor, pybind11::arg( "num_overlays" ),
              R"pbdoc(
    Convert a given number of overlay blocks to the factor f for the currently added points.

    Here f is a factor proportional to the number of boxes used in the data structure.
    See estimate_num_elements for a detailed description. 

    :param num_overlays: number of overlay blocks
    :type f: int

    :return: f.
    :rtype: float
)pbdoc" )
        .def( "grid_size", &sps::Index<type_defs>::gridSize,
              R"pbdoc(
    Get the axis sizes of the area spanned by the currently added points.

    :return: a.
    :rtype: list[int]
)pbdoc" )
        .def( "get_num_internal_prefix_sums", &sps::Index<type_defs>::getNumInternalPrefixSums,
              pybind11::arg( "dataset_id" ),
              R"pbdoc(
    Count the number of internal prefix sums stored in the dataset with id dataset_id.

    :param dataset_id: The id of the dataset to query
    :type dataset_id: int

    :return: number of internal prefix sums.
    :rtype: int
)pbdoc" )
        .def( "get_num_overlay_prefix_sums", &sps::Index<type_defs>::getNumOverlayPrefixSums,
              pybind11::arg( "dataset_id" ),
              R"pbdoc(
    Count the number of overlay prefix sums stored in the dataset with id dataset_id.

    :param dataset_id: The id of the dataset to query
    :type dataset_id: int

    :return: number of overlay prefix sums.
    :rtype: int
)pbdoc" )
        .def( "get_num_internal_sparse_coords", &sps::Index<type_defs>::getNumInternalSparseCoords,
              pybind11::arg( "dataset_id" ),
              R"pbdoc(
    Count the number of internal sparse coordinates stored in the dataset with id dataset_id.

    :param dataset_id: The id of the dataset to query
    :type dataset_id: int

    :return: number of internal sparse coordinates.
    :rtype: int
)pbdoc" )
        .def( "get_num_overlay_sparse_coords", &sps::Index<type_defs>::getNumOverlaySparseCoords,
              pybind11::arg( "dataset_id" ),
              R"pbdoc(
    Count the number of overlay sparse coordinates stored in the dataset with id dataset_id.

    :param dataset_id: The id of the dataset to query
    :type dataset_id: int

    :return: number of overlay sparse coordinates.
    :rtype: int
)pbdoc" )


        ;
    return "    " + sName + "\n";
}
#endif