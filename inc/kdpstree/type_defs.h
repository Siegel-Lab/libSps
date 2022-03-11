#pragma once

#include <array>
#include <vector>

namespace kdpstree
{

template <typename _coordinate_t, //
          typename _val_t, //
          typename _layers_t, //
          size_t _layers, //
          typename _class_key_t, //
          template <typename> typename _tmp_vec_generator, //
          template <typename> typename _vec_generator, //
          template <typename, typename> typename _sort_func_t, //
          size_t _b, //
          typename _offset_t, //
          bool _explain //
          >
class TypeDefs
{
  public:
    using coordinate_t = _coordinate_t;
    using val_t = _val_t;
    using layers_t = _layers_t;
    static const layers_t LAYERS = _layers;
    // using cnt_t = std::array<val_t, 2>;
    using data_t = std::array<val_t, LAYERS>;
    static const coordinate_t d = 2; // @todo remove
    using pos_t = std::array<coordinate_t, d>;
    using class_key_t = _class_key_t;

    using overlay_key_t = std::pair<class_key_t, pos_t>;

    template <typename val_type_t> using tmp_vec_generator = _tmp_vec_generator<val_type_t>;
    template <typename val_type_t> using vec_generator_t = _vec_generator<val_type_t>;

    template <typename it_t, typename cmp_t> using sort_func_t = _sort_func_t<it_t, cmp_t>;

    using overlay_entry_t = std::pair<coordinate_t, data_t>;

    static constexpr bool EXPLAIN_QUERY = _explain;

    static const size_t b = _b;

    using offset_t = _offset_t;
};

#define EXTRACT_TYPE_DEFS                                                                                              \
    using coordinate_t = typename type_defs::coordinate_t;                                                             \
                                                                                                                       \
    using val_t = typename type_defs::val_t;                                                                           \
                                                                                                                       \
    using layers_t = typename type_defs::layers_t;                                                                     \
                                                                                                                       \
    /*using cnt_t = typename type_defs::cnt_t;*/                                                                       \
                                                                                                                       \
    using data_t = typename type_defs::data_t;                                                                         \
                                                                                                                       \
    static const size_t LAYERS = type_defs::LAYERS;                                                                    \
                                                                                                                       \
    static const coordinate_t d = type_defs::d;                                                                        \
                                                                                                                       \
    using pos_t = typename type_defs::pos_t;                                                                           \
                                                                                                                       \
    using class_key_t = typename type_defs::class_key_t;                                                               \
                                                                                                                       \
    using overlay_key_t = typename type_defs::overlay_key_t;                                                           \
                                                                                                                       \
    template <typename val_type_t>                                                                                     \
    using tmp_vec_generator = typename type_defs::template tmp_vec_generator<val_type_t>;                              \
                                                                                                                       \
    template <typename val_type_t> using vec_generator_t = typename type_defs::template vec_generator_t<val_type_t>;   \
                                                                                                                       \
    template <typename it_t, typename cmp_t>                                                                           \
    using sort_func_t = typename type_defs::template sort_func_t<it_t, cmp_t>;                                         \
                                                                                                                       \
    using overlay_entry_t = typename type_defs::overlay_entry_t;                                                       \
                                                                                                                       \
    static constexpr bool EXPLAIN_QUERY = type_defs::EXPLAIN_QUERY;                                                    \
                                                                                                                       \
    static const size_t b = type_defs::b;                                                                              \
                                                                                                                       \
    using offset_t = typename type_defs::offset_t;


#define EXTRACT_VEC_GENERATOR( name, content_t )                                                                       \
                                                                                                                       \
    using name##_vec_generator_t = vec_generator_t<content_t>;                                                         \
                                                                                                                       \
    name##_vec_generator_t name##_vec_generator = name##_vec_generator_t( );                                           \
                                                                                                                       \
    using name##_vec_t = typename name##_vec_generator_t::vec_t;                                                       \
                                                                                                                       \
    using name##_file_t = typename name##_vec_generator_t::file_t

#define EXTRACT_TMP_VEC_GENERATOR( name, content_t )                                                                   \
                                                                                                                       \
    using name##_vec_generator_t = tmp_vec_generator<content_t>;                                                       \
                                                                                                                       \
    name##_vec_generator_t name##_vec_generator = name##_vec_generator_t( );                                           \
                                                                                                                       \
    using name##_vec_t = typename name##_vec_generator_t::vec_t;

} // namespace kdpstree