#pragma once

#include "sps/type_defs.h"

namespace sps
{

template <typename type_defs> class BaseCorner
{
    EXTRACT_TYPE_DEFS; // macro call

    using desc_t = Desc<type_defs>;

  public:
    pos_t vPos;
    val_t uiVal;

    BaseCorner( pos_t vPos, val_t uiVal ) : vPos( vPos ), uiVal( uiVal )
    {}

    BaseCorner( ) : vPos{ }, uiVal( 0 )
    {}

    friend std::ostream& operator<<( std::ostream& os, const BaseCorner& xCorner )
    {
        os << xCorner.vPos;
        return os;
    }

    std::ostream& stream( std::ostream& os ) const
    {
        os << vPos;
        return os;
    }

    void addTo( val_t& uiTo ) const
    {
        uiTo += uiVal;
    }
};


// Orthotope == HyperRectangle
template <typename type_defs> class OrthotopeCorner : public BaseCorner<type_defs>
{
    EXTRACT_TYPE_DEFS; // macro call

    using desc_t = Desc<type_defs>;

    uint8_t uiIdx;

  public:
    OrthotopeCorner( pos_t vPos, val_t uiVal, uint8_t uiIdx ) : BaseCorner<type_defs>( vPos, uiVal ), uiIdx( uiIdx )
    {}

    OrthotopeCorner( pos_t vPos, val_t uiVal ) : BaseCorner<type_defs>( vPos, uiVal ), uiIdx( 0 )
    {}

    OrthotopeCorner( ) : BaseCorner<type_defs>( ), uiIdx( 0 )
    {}

    std::ostream& stream( std::ostream& os ) const
    {
        BaseCorner<type_defs>::stream( os ) << " i" << (size_t)uiIdx;
        return os;
    }

    void addTo( sps_t& uiTo ) const
    {
        BaseCorner<type_defs>::addTo( uiTo[ uiIdx ] );
    }
};

#define USED_POINT std::conditional<type_defs::IS_ORTHOTOPE, OrthotopeCorner<type_defs>, BaseCorner<type_defs>>::type

// conditional inheritance required to force minimal memory usage (memory of disabled types is still allocated)
template <typename type_defs> class Corner : public USED_POINT
{
    using X = typename USED_POINT;

  public:
    using X::X;

}; // class


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter" // do not warn about vFrom and vTo

template <typename pos_t, size_t N, size_t NE>
inline void forAllCombinationsHelper( std::function<void( size_t, pos_t, size_t )> fDo, pos_t& vCurr, pos_t vFrom,
                                      pos_t vTo, size_t uiDistTo, size_t uiNum,
                                      std::function<bool( typename pos_t::value_type )> fCond )
{
    if constexpr /* <- required to prevent infinite unrolling loop in compiler */ ( N == NE )
        fDo( uiNum, vCurr, uiDistTo );
    else
    {
        vCurr[ N ] = vFrom[ N ];
        if( fCond( vCurr[ N ] ) )
            forAllCombinationsHelper<pos_t, N + 1, NE>( fDo, vCurr, vFrom, vTo, uiDistTo + 1, uiNum, fCond );
        vCurr[ N ] = vTo[ N ];
        if( fCond( vCurr[ N ] ) )
            forAllCombinationsHelper<pos_t, N + 1, NE>( fDo, vCurr, vFrom, vTo, uiDistTo,
                                                        uiNum + ( 1 << ( NE - ( N + 1 ) ) ), fCond );
    }
}

#pragma GCC diagnostic pop

template <typename pos_t, size_t N>
inline void forAllCombinationsN(
    std::function<void( size_t, pos_t, size_t )> fDo,
    pos_t vFrom,
    pos_t vTo,
    std::function<bool( typename pos_t::value_type )> fCond = []( typename pos_t::value_type ) { return true; } )
{
    pos_t vCurr{ };
    forAllCombinationsHelper<pos_t, 0, N>( fDo, vCurr, vFrom, vTo, 0, 0, fCond );
}

template <typename> struct array_size;
template <typename T, size_t N> struct array_size<std::array<T, N>>
{
    static size_t const size = N;
};
template <typename pos_t>
inline void forAllCombinations(
    std::function<void( size_t, pos_t, size_t )> fDo,
    pos_t vFrom,
    pos_t vTo,
    std::function<bool( typename pos_t::value_type )> fCond = []( typename pos_t::value_type ) { return true; } )
{
    forAllCombinationsN<pos_t, array_size<pos_t>::size>( fDo, vFrom, vTo, fCond );
}


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter" // do not warn about vFrom and vTo

template <typename pos_t, size_t N, size_t NE, typename fDo_t, typename fCond_t, typename... extra_t>
inline __attribute__( ( always_inline ) ) void forAllCombinationsHelperTmpl( pos_t& vCurr, pos_t& vFrom, pos_t& vTo,
                                                                             size_t uiDistTo, size_t uiNum, fDo_t&& fDo,
                                                                             fCond_t&& fCond, extra_t&&... rExtra )
{
    if constexpr /* <- required to prevent infinite unrolling of loop in compiler */ ( N == NE )
        fDo( uiNum, vCurr, uiDistTo, rExtra... );
    else
    {
        vCurr[ N ] = vFrom[ N ];
        if( fCond( vCurr[ N ], N, true ) )
            forAllCombinationsHelperTmpl<pos_t, N + 1, NE>( vCurr, vFrom, vTo, uiDistTo + 1, uiNum, fDo, fCond,
                                                            rExtra... );
        vCurr[ N ] = vTo[ N ];
        if( fCond( vCurr[ N ], N, false ) )
            forAllCombinationsHelperTmpl<pos_t, N + 1, NE>(
                vCurr, vFrom, vTo, uiDistTo, uiNum + ( 1 << ( NE - ( N + 1 ) ) ), fDo, fCond, rExtra... );
    }
}
#pragma GCC diagnostic pop

template <typename pos_t, size_t N, typename fDo_t, typename fCond_t, typename... extra_t>
inline __attribute__( ( always_inline ) ) void forAllCombinationsNTmpl( pos_t& vFrom, pos_t& vTo, fDo_t&& fDo,
                                                                        fCond_t&& fCond, extra_t&&... rExtra )
{
    pos_t vCurr{ };
    forAllCombinationsHelperTmpl<pos_t, 0, N>( vCurr, vFrom, vTo, 0, 0, fDo, fCond, rExtra... );
}

template <typename pos_t, typename fDo_t, typename fCond_t, typename... extra_t>
inline __attribute__( ( always_inline ) ) void forAllCombinationsTmpl( pos_t& vFrom, pos_t& vTo, fDo_t&& fDo,
                                                                       fCond_t&& fCond, extra_t&&... rExtra )
{
    forAllCombinationsNTmpl<pos_t, array_size<pos_t>::size>( vFrom, vTo, fDo, fCond, rExtra... );
}

} // namespace sps
