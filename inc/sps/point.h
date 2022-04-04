#pragma once

#include "sps/type_defs.h"

namespace sps
{
template <typename type_defs> class Point
{
    EXTRACT_TYPE_DEFS; // macro call
    
    using desc_t = Desc<type_defs>;

  public:
    pos_t vPos;
    size_t uiDescOffset;

    Point( pos_t vPos, size_t uiDescOffset )
        : vPos( vPos ), uiDescOffset( uiDescOffset )
    {}

    Point( ) : vPos{ }, uiDescOffset( 0 )
    {}

    friend std::ostream& operator<<( std::ostream& os, const Point& xPoint )
    {
        os << xPoint.vPos << " d" << xPoint.uiDescOffset;
        return os;
    }

    std::ostream& stream( std::ostream& os, const desc_t& vDesc ) const
    {
        os << vPos << ": " << vDesc.get(uiDescOffset);
        return os;
    }
};


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter" // do not warn about vFrom and vTo

template<typename pos_t, size_t N>
inline void forAllCombinationsHelper( std::function<void(pos_t, size_t)> fDo, 
                                      pos_t& vCurr, 
                                      pos_t vFrom, 
                                      pos_t vTo, size_t uiDistTo, 
                                      std::function<bool(typename pos_t::value_type)> fCond )
{
    if constexpr /* <- required to prevent infinite unrolling loop in compiler */(N == vCurr.size())
        fDo(vCurr, uiDistTo);
    else
    {
        vCurr[N] = vFrom[N];
        if( fCond(vCurr[N]) )
            forAllCombinationsHelper<pos_t, N+1>(fDo, vCurr, vFrom, vTo, uiDistTo + 1, fCond);
        vCurr[N] = vTo[N];
        if( fCond(vCurr[N]) )
            forAllCombinationsHelper<pos_t, N+1>(fDo, vCurr, vFrom, vTo, uiDistTo, fCond);
    }
}
    
#pragma GCC diagnostic pop

template<typename pos_t>
inline void forAllCombinations( std::function<void(pos_t, size_t)> fDo, 
                                pos_t vFrom, 
                                pos_t vTo,
                                std::function<bool(typename pos_t::value_type)> fCond = 
                                    [](typename pos_t::value_type){ return true; } )
{
    pos_t vCurr;
    forAllCombinationsHelper<pos_t, 0>(fDo, vCurr, vFrom, vTo, 0, fCond);
}

} // namespace sps
