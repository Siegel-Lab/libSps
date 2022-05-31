#include <array>

namespace std
{

template <class T, std::size_t N> class op_array : public array<T, N>
{
  public:
    using array<T, N>::array;

    op_array& operator+=( const op_array& rOther )
    {
        for( size_t uiI = 0; uiI < N; uiI++ )
            ( *this )[ uiI ] += rOther[ uiI ];
        return *this;
    }

    friend op_array& operator*( op_array& rMe, T rOther )
    {
        for( size_t uiI = 0; uiI < N; uiI++ )
            rMe[ uiI ] *= rOther;
        return rMe;
    }
};

} // namespace std