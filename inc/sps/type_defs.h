#pragma once

#include "sps/operator_array.h"
#include <array>
#include <vector>

namespace sps
{

// compile time switch for the progress prints in the counting code (to properly optimize speed)
#ifdef NDEBUG
#define GET_PROG_PRINTS 0
#else
#define GET_PROG_PRINTS 1
#endif

#define DU_UNREALISTIC_VALUE_CHECK 0

#define ALWAYS_SIMULATE_GRID_QUERY true

class Verbosity
{
  public:
    size_t uiI;
    Verbosity( size_t uiI ) : uiI( uiI )
    {}

    bool operator>( const Verbosity& rOther ) const
    {
        return uiI > rOther.uiI;
    }
};

#define NAMED_VEC_GEN_AND_SORTER_TEMPLATE( name ) template <typename> typename _##name##_vec_generator

#define NAMED_VEC_GEN_AND_SORTER( name )                                                                               \
    template <typename val_type_t> using name##_vec_generator_t = _##name##_vec_generator<val_type_t>;


#define NAMED_VEC_GEN_AND_SORTER_EXTRACT( name )                                                                       \
    template <typename val_type_t>                                                                                     \
    using name##_tmpl_vec_generator_t = typename type_defs::template name##_vec_generator_t<val_type_t>;


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
 * @brief Definition of all compiletime parameters for Index.
 *
 * The Index class takes this as a template parameter.
 *
 * @todo a lot of the template parameters are not actually variable most of the time
 *      - having them drastically increases compiletime
 *      - also the TypeDefs class is kind of a bad idea since it needs to be instanciated every time
 *      - better approach:
 *          - create one big class/struct with the template parameters
 *          - then make everything a subclass of this one
 *          - split the commonly used and unused parameters into tow subclasses
 *
 * @tparam _coordinate_t type of coordinates, expected to be an unsigned int.
 * @tparam _val_t type of prefix sums, expected to be an unsigned int.
 * @tparam _D number of dimensions.
 * @tparam _class_key_t type of the dataset keys, expected to be an unsigned int.
 * @tparam _vec_generator generator object that implements the file and vec methods. See DiskVecGenerator.
 * @tparam _binary_search_based_sparse whether sparse space lookup tables shall be binary search based
 * @tparam _orthotope_dims number or orthotope dimensions
 * @tparam _explain debugging parameter. Be verbose while creating and querying the index.
 * @tparam _progress_stream_t object that catches all the print output. See StdOutProgressStream.
 */
template <typename _coordinate_t, //
          typename _val_t, //
          size_t _D, //
          typename _class_key_t, //
          NAMED_VEC_GEN_AND_SORTER_TEMPLATE( dataset ), //
          NAMED_VEC_GEN_AND_SORTER_TEMPLATE( coord ), //
          NAMED_VEC_GEN_AND_SORTER_TEMPLATE( overlays ), //
          NAMED_VEC_GEN_AND_SORTER_TEMPLATE( desc ), //
          NAMED_VEC_GEN_AND_SORTER_TEMPLATE( points ), //
          NAMED_VEC_GEN_AND_SORTER_TEMPLATE( prefix_sums ), //
          bool _binary_search_based_sparse, //
          size_t _orthotope_dims, //
          bool _explain, //
          typename _progress_stream_t>
class TypeDefs
{
  public:
    /// @brief individual coordinate position
    using coordinate_t = _coordinate_t;

    /// @brief prefix sum value
    using val_t = _val_t;
    static constexpr coordinate_t D = _D;

    /// @brief position of a point (with dependent dimensions; used internally)
    using pos_t = std::array<coordinate_t, D>;

    /// @brief key of a dataset
    using class_key_t = _class_key_t;

    NAMED_VEC_GEN_AND_SORTER( dataset )
    NAMED_VEC_GEN_AND_SORTER( coord )
    NAMED_VEC_GEN_AND_SORTER( overlays )
    NAMED_VEC_GEN_AND_SORTER( desc )
    NAMED_VEC_GEN_AND_SORTER( points )
    NAMED_VEC_GEN_AND_SORTER( prefix_sums )

    static constexpr bool EXPLAIN = _explain;


    static constexpr bool BINARY_SEARCH_BASED_SPARSE = _binary_search_based_sparse;

    static constexpr coordinate_t ORTHOTOPE_DIMS = _orthotope_dims;

    /// @brief position of a point
    using ret_pos_t = std::array<coordinate_t, D - ORTHOTOPE_DIMS>;

    static constexpr bool IS_ORTHOTOPE = ORTHOTOPE_DIMS > 0;

    using sps_t = typename std::conditional<IS_ORTHOTOPE, std::op_array<val_t, 1 << ORTHOTOPE_DIMS>, val_t>::type;

    using progress_stream_t = _progress_stream_t;

    using isect_arr_t = std::array<IntersectionType, ORTHOTOPE_DIMS>;
};

#define EXTRACT_TYPE_DEFS                                                                                              \
    /** @brief individual coordinate position */                                                                       \
    using coordinate_t = typename type_defs::coordinate_t;                                                             \
                                                                                                                       \
    /** @brief prefix sum value */                                                                                     \
    using val_t = typename type_defs::val_t;                                                                           \
                                                                                                                       \
    static constexpr coordinate_t D = type_defs::D;                                                                    \
                                                                                                                       \
    /** @brief position of a point */                                                                                  \
    using ret_pos_t = typename type_defs::ret_pos_t;                                                                   \
                                                                                                                       \
    /** @brief position of a point (with dependent dimensions; used internally) */                                     \
    using pos_t = typename type_defs::pos_t;                                                                           \
                                                                                                                       \
    /** @brief key of a dataset */                                                                                     \
    using class_key_t = typename type_defs::class_key_t;                                                               \
                                                                                                                       \
    NAMED_VEC_GEN_AND_SORTER_EXTRACT( dataset )                                                                        \
    NAMED_VEC_GEN_AND_SORTER_EXTRACT( coord )                                                                          \
    NAMED_VEC_GEN_AND_SORTER_EXTRACT( overlays )                                                                       \
    NAMED_VEC_GEN_AND_SORTER_EXTRACT( desc )                                                                           \
    NAMED_VEC_GEN_AND_SORTER_EXTRACT( points )                                                                         \
    NAMED_VEC_GEN_AND_SORTER_EXTRACT( prefix_sums )                                                                    \
                                                                                                                       \
    static constexpr bool EXPLAIN = type_defs::EXPLAIN;                                                                \
                                                                                                                       \
    static constexpr bool BINARY_SEARCH_BASED_SPARSE = type_defs::BINARY_SEARCH_BASED_SPARSE;                          \
                                                                                                                       \
    static constexpr coordinate_t ORTHOTOPE_DIMS = type_defs::ORTHOTOPE_DIMS;                                          \
                                                                                                                       \
    static constexpr bool IS_ORTHOTOPE = type_defs::IS_ORTHOTOPE;                                                      \
                                                                                                                       \
    static constexpr bool IS_ORTHOTOPE_CAT = type_defs::IS_ORTHOTOPE_CAT;                                              \
                                                                                                                       \
    using sps_t = typename type_defs::sps_t;                                                                           \
                                                                                                                       \
    using progress_stream_t = typename type_defs::progress_stream_t;                                                   \
                                                                                                                       \
    using isect_arr_t = typename type_defs::isect_arr_t;


#define EXTRACT_VEC_GENERATOR( name, content_t )                                                                       \
    static_assert( 4096 % sizeof( content_t ) == 0 );                                                                  \
                                                                                                                       \
    using name##_vec_generator_t = name##_tmpl_vec_generator_t<content_t>;                                             \
                                                                                                                       \
    name##_vec_generator_t name##_vec_generator = name##_vec_generator_t( );                                           \
                                                                                                                       \
    using name##_vec_t = typename name##_vec_generator_t::vec_t;                                                       \
                                                                                                                       \
    using name##_file_t = typename name##_vec_generator_t::file_t;                                                     \
                                                                                                                       \
    static const bool name##_THREADSAVE = name##_vec_generator_t::THREADSAVE;                                          \
                                                                                                                       \
    template <typename it_t, typename cmp_t>                                                                           \
    using name##_sort_func_t = typename name##_vec_generator_t::template sorter_t<it_t, cmp_t>;


} // namespace sps