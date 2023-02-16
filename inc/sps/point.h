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


template <typename pos_t, size_t N, size_t NE>
inline void forAllCombinationsHelper( std::function<void( size_t, pos_t, size_t )> fDo, pos_t& vCurr,
                                      [[maybe_unused]] pos_t vFrom, [[maybe_unused]] pos_t vTo, size_t uiDistTo,
                                      size_t uiNum, std::function<bool( typename pos_t::value_type )> fCond )
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


template <typename pos_t, size_t N, size_t NE, template <size_t, size_t, size_t> class fDo_t, size_t uiDistTo,
          size_t uiNum, size_t uiFirstZero, typename fCond_t, typename... extra_t>
inline void forAllCombinationsHelperTmpl( pos_t& vCurr, [[maybe_unused]] pos_t& vFrom, [[maybe_unused]] pos_t& vTo,
                                          fCond_t&& fCond, extra_t&&... rExtra )
{
    if constexpr /* <- required to prevent infinite unrolling of loop in compiler */ ( N == NE )
        fDo_t<uiNum, uiDistTo, uiFirstZero>::count( vCurr, rExtra... );
    else
    {
        vCurr[ N ] = vFrom[ N ];
        if( fCond( vCurr[ N ], N, true ) )
        {
            if constexpr( uiFirstZero == NE )
                forAllCombinationsHelperTmpl<pos_t, N + 1, NE, fDo_t, uiDistTo + 1, uiNum, N>( vCurr, vFrom, vTo, fCond,
                                                                                               rExtra... );
            else
                forAllCombinationsHelperTmpl<pos_t, N + 1, NE, fDo_t, uiDistTo + 1, uiNum, uiFirstZero>(
                    vCurr, vFrom, vTo, fCond, rExtra... );
        }
        vCurr[ N ] = vTo[ N ];
        if( fCond( vCurr[ N ], N, false ) )
            forAllCombinationsHelperTmpl<pos_t, N + 1, NE, fDo_t, uiDistTo, uiNum + ( 1 << ( NE - ( N + 1 ) ) ),
                                         uiFirstZero>( vCurr, vFrom, vTo, fCond, rExtra... );
    }
}

template <typename pos_t, size_t N, template <size_t, size_t, size_t> class fDo_t, typename fCond_t,
          typename... extra_t>
inline void forAllCombinationsNTmpl( pos_t& vFrom, pos_t& vTo, fCond_t&& fCond, extra_t&&... rExtra )
{
    pos_t vCurr{ };
    forAllCombinationsHelperTmpl<pos_t, 0, N, fDo_t, 0, 0, N>( vCurr, vFrom, vTo, fCond, rExtra... );
}

template <typename pos_t, template <size_t, size_t, size_t> class fDo_t, typename fCond_t, typename... extra_t>
inline void forAllCombinationsTmpl( pos_t& vFrom, pos_t& vTo, fCond_t&& fCond, extra_t&&... rExtra )
{
    forAllCombinationsNTmpl<pos_t, array_size<pos_t>::size, fDo_t>( vFrom, vTo, fCond, rExtra... );
}

} // namespace sps
