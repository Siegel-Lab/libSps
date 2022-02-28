#pragma once

#include <array>
#include <vector>

namespace bps
{

template <typename _coordinate_t, //
          typename _val_t, //
          size_t _d, //
          typename _class_key_t, //
          template <typename> typename _vec_generator, //
          template <typename> typename _map_generator, //
          template <typename, typename> typename _sort_func_t //
          >
class TypeDefs
{
  public:
    using coordinate_t = _coordinate_t;
    using val_t = _val_t;
    static const _coordinate_t d = _d;
    using pos_t = std::array<coordinate_t, d>;
    using class_key_t = _class_key_t;

    using overlay_key_t = std::pair<class_key_t, pos_t>;

    template <typename val_type_t> using vec_generator = _vec_generator<val_type_t>;
    template <typename key_type_t, typename val_type_t> using map_generator = _map_generator<key_type_t, val_type_t>;

    template <typename it_t, typename cmp_t> using sort_func_t = _sort_func_t<it_t, cmp_t>;

    using overlay_entry = std::pair<coordinate_t, val_t>;
};

} // namespace bps