#pragma once

#include "sps/abstract_index.h"
#include "sps/desc.h"
#include "sps/point.h"
#include "sps/type_defs.h"
#include "sps/util.h"

#include <cassert>
#include <functional>
#include <string>

#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif


namespace sps
{

/**
 * @brief A naive, iteration-based index implementation for benchmarking.
 *
 * Ignores the number of orthotope dimensions ( always sets type_defs::ORTHOTOPE_DIMS = type_defs:D ).
 *
 * @tparam type_defs An instance of TypeDefs, that defines all compiletime parameters
 */
template <typename type_defs> class SimpleVector : public AbstractIndex
{
    EXTRACT_TYPE_DEFS; // macro call

    using base_corner_t = BaseCorner<type_defs>;
    using ortho_t = AlignedPower2<std::array<base_corner_t, 2>>;


  public:
    EXTRACT_VEC_GENERATOR( points, ortho_t ); // macro call
    points_file_t xFile;
    points_vec_t vData;

    /**
     * @brief Construct a new SimpleVector object
     *
     * @param sPrefix Prefix path of the vector on the filesystem (a file with the ending .orthos will be created),
     * defaults to "".
     * @param bWrite Open the vector in write mode (if this is set to False no changes can be made to the vector),
     * defaults to false.
     */
    SimpleVector( std::string sPrefix = "", bool bWrite = false )
        : xFile( points_vec_generator.file( sPrefix + ".orthos", bWrite ) ), vData( points_vec_generator.vec( xFile ) )
    {}

    /**
     * @brief Append an orthotope to the data structure.
     *
     * @param vStart The bottom-left of the orthotope.
     * @param vEnd The top-right of the orthotope.
     * @param uiVal Value of the points, default to 1.
     */
    void addPoint( pos_t vStart, pos_t vEnd, val_t uiVal = 1 )
    {
        ortho_t xNew;
        xNew[ 0 ] = base_corner_t( vStart, uiVal );
        xNew[ 1 ] = base_corner_t( vEnd, uiVal );
        vData.push_back( xNew );
    }

    /**
     * @brief Count the number of orthotopes within the given area.
     *
     * @param vStart The bottom-left of the queried area.
     * @param vEnd The top-right of the queried area.
     */
    val_t count( pos_t vStart, pos_t vEnd ) const
    {
        val_t uiRet = 0;
        for( const ortho_t& rP : vData )
        {
            bool bContinue = true;
            for( size_t uiI = 0; uiI < D; uiI++ )
                bContinue = bContinue && vStart[ uiI ] <= rP[ 0 ].vPos[ uiI ] && vEnd[ uiI ] >= rP[ 1 ].vPos[ uiI ];
            if( bContinue )
                uiRet += 1;
        }
        return uiRet;
    }
};

/**
 * @brief A simple vector implementation for benchmarking.
 *
 * @tparam type_defs An instance of TypeDefs, that defines all compiletime parameters
 */
template <typename type_defs> class SimpleValVector : public AbstractIndex
{
    EXTRACT_TYPE_DEFS; // macro call


  public:
    EXTRACT_VEC_GENERATOR( prefix_sums, val_t ); // macro call
    prefix_sums_file_t xFile;
    prefix_sums_vec_t vData;

    /**
     * @brief Construct a new SimpleValVector object
     *
     * @param sPrefix Prefix path of the vector on the filesystem (a file with the ending .vals will be created),
     * defaults to "".
     * @param bWrite Open the vector in write mode (if this is set to False no changes can be made to the vector),
     * defaults to false.
     */
    SimpleValVector( std::string sPrefix = "", bool bWrite = false )
        : xFile( prefix_sums_vec_generator.file( sPrefix + ".vals", bWrite ) ),
          vData( prefix_sums_vec_generator.vec( xFile ) )
    {}

    /**
     * @brief Append a value.
     *
     * @param uiV The value to append.
     */
    void add( val_t uiV )
    {
        vData.push_back( uiV );
    }

    /**
     * @brief query a value.
     *
     * @param uiX position to query.
     */
    val_t get( size_t uiX ) const
    {
        return vData[ uiX ];
    }

    std::vector<val_t> getMultiple( std::vector<size_t> vX ) const
    {
        std::vector<val_t> vRet( vX.size( ) );
        for( size_t uiI = 0; uiI < vX.size( ); uiI++ )
            vRet[ uiI ] = get( vX[ uiI ] );

        return vRet;
    }
};


} // namespace sps

#if WITH_PYTHON
template <typename type_defs> std::string exportSimpleVector( pybind11::module& m, std::string sName )
{
    pybind11::class_<sps::SimpleValVector<type_defs>>( m, sName.c_str( ),
                                                       R"pbdoc(
    Simple Vector for val_t entries.
    
    .. automethod:: add
    .. automethod:: get
    .. automethod:: get_multiple
)pbdoc" )
        .def( pybind11::init<std::string, bool>( ),
              pybind11::arg( "path" ) = "",
              pybind11::arg( "write_mode" ) = false,
              R"pbdoc(
    Create a new vector.

    :param path: Prefix path of the vector on the filesystem (a file with the ending .orthos will be created), defaults to "".
    :type path: str
    
    :param write_mode: Open the vector in write mode (if this is set to False no changes can be made to the vector), defaults to False.
    :type write_mode: str
)pbdoc" )
        .def( "add", &sps::SimpleValVector<type_defs>::add, pybind11::arg( "v" ),
              R"pbdoc(
    Append a value.

    :param v: The value.
    :type v: int
)pbdoc" )
        .def( "get", &sps::SimpleValVector<type_defs>::get, pybind11::arg( "x" ),
              R"pbdoc(
    Query a value.

    :param x: index of the value to query
    :type x: int
)pbdoc" )
        .def( "get_multiple", &sps::SimpleValVector<type_defs>::getMultiple, pybind11::arg( "xs" ), ( R"pbdoc(
    Query multiple values.
    Multithreaded if the underlying storage is threadsave.

    :param xs: list of indices that shall be queried.
    :type xs: list(int)
)pbdoc" ) )

        ;

    return "    " + sName + "\n ";
}
#endif
