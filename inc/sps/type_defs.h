#pragma once

#include <array>
#include <vector>

namespace sps
{

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

template <typename _coordinate_t, //
          typename _val_t, //
          size_t _D, //
          typename _class_key_t, //
          template <typename> typename _vec_generator, //
          template <typename, typename> typename _sort_func_t, //
          bool _dependant_dim, //
          bool _explain, //
          typename _progress_stream_t>
class TypeDefs
{
  public:
    using coordinate_t = _coordinate_t;
    using val_t = _val_t;
    static constexpr coordinate_t D = _D;
    using pos_t = std::array<coordinate_t, D>;
    using class_key_t = _class_key_t;

    template <typename val_type_t> using vec_generator_t = _vec_generator<val_type_t>;

    template <typename it_t, typename cmp_t> using sort_func_t = _sort_func_t<it_t, cmp_t>;

    static constexpr bool EXPLAIN = _explain;

    static constexpr bool DEPENDANT_DIMENSION = _dependant_dim;

    using progress_stream_t = _progress_stream_t;
};

#define EXTRACT_TYPE_DEFS                                                                                              \
    using coordinate_t = typename type_defs::coordinate_t;                                                             \
                                                                                                                       \
    using val_t = typename type_defs::val_t;                                                                           \
                                                                                                                       \
    static constexpr coordinate_t D = type_defs::D;                                                                    \
                                                                                                                       \
    using pos_t = typename type_defs::pos_t;                                                                           \
                                                                                                                       \
    using class_key_t = typename type_defs::class_key_t;                                                               \
                                                                                                                       \
    template <typename val_type_t> using vec_generator_t = typename type_defs::template vec_generator_t<val_type_t>;   \
                                                                                                                       \
    template <typename it_t, typename cmp_t>                                                                           \
    using sort_func_t = typename type_defs::template sort_func_t<it_t, cmp_t>;                                         \
                                                                                                                       \
    static constexpr bool EXPLAIN = type_defs::EXPLAIN;                                                                \
                                                                                                                       \
    static constexpr bool DEPENDANT_DIMENSION = type_defs::DEPENDANT_DIMENSION;                                        \
                                                                                                                       \
    using progress_stream_t = typename type_defs::progress_stream_t;


#define EXTRACT_VEC_GENERATOR( name, content_t )                                                                       \
    static_assert( 4096 % sizeof( content_t ) == 0 );                                                                  \
                                                                                                                       \
    using name##_vec_generator_t = vec_generator_t<content_t>;                                                         \
                                                                                                                       \
    name##_vec_generator_t name##_vec_generator = name##_vec_generator_t( );                                           \
                                                                                                                       \
    using name##_vec_t = typename name##_vec_generator_t::vec_t;                                                       \
                                                                                                                       \
    using name##_file_t = typename name##_vec_generator_t::file_t;                                                     \
                                                                                                                       \
    static const bool name##_THREADSAVE = name##_vec_generator_t::THREADSAVE;


} // namespace sps