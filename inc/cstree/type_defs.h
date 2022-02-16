#pragma once

#include <array>
#include <vector>

namespace cstree
{

template <typename coordinate_t>
using bin_cords_generator_t = std::vector<coordinate_t> ( * )( coordinate_t, coordinate_t, coordinate_t, coordinate_t );

template <typename vec_t> using vec_generator_t = vec_t ( * )( std::string );

template <typename it_t, typename cmp_t> using vec_sorter_t = void ( * )( it_t, it_t, cmp_t );

template <typename _coordinate_t, //
          typename _cont_sum_val_t, //
          typename _dimensions_t, //
          typename _points_vec_offset_t, //
          typename _sums_vec_offset_t, //
          size_t _d, //
          template <typename> class _bin_cords_generator, //
          template <typename> class _vec_generator, //
          template <typename> class _data_points_vec_t, //
          template <typename> class _cont_sums_vec_t, //
          template <typename> class _point_desc_vec_t, //
          template <typename> class _bin_coords_vec_t, //
          template <typename> class _sort_func_t //
          >
class TypeDefs
{
  public:
    using coordinate_t = _coordinate_t;
    using cont_sum_val_t = _cont_sum_val_t;
    using dimensions_t = _dimensions_t;
    using points_vec_offset_t = _points_vec_offset_t;
    using sums_vec_offset_t = _sums_vec_offset_t;
    static const _coordinate_t d = _d;
    template <typename data_point_t> using data_points_vec_t = _data_points_vec_t<data_point_t>;
    using point_t = std::array<coordinate_t, d>;
    using point_desc_vec_t = _point_desc_vec_t<char>;
    template <typename cont_sum_t> using cont_sums_vec_t = _cont_sums_vec_t<cont_sum_t>;
    using bin_coords_vec_t = _bin_coords_vec_t<coordinate_t>;
    using bin_cords_generator = _bin_cords_generator<coordinate_t>;
    using vec_generator = _vec_generator;
    using sort_func_t = _sort_func_t;
}

} // namespace cstree