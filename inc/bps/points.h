#pragma once

#include "bps/desc.h"
#include "bps/point.h"
#include "bps/type_defs.h"
#include <string>


namespace bps
{

template <typename type_defs> class Points
{
    using vec_generator = typename type_defs::vec_generator;
    using point = Point<type_defs>;
    auto points_vec_generator = vec_generator<point>( );
    using points_vec = typeof( points_vec_generator( "" ) );
    using points_it_t = typename points_vec::const_iterator;
    using coordinate_t = typename type_defs::coordinate_t;

    struct PointsComperator
    {
        const size_t uiDim;

        PointsComperator( size_t uiDim ) : uiDim( uiDim )
        {}

        bool operator( )( point a, point b ) const
        {
            return a.vPos[ uiDim ] < b.vPos[ uiDim ];
        }

        coordinate_t min_value( ) const
        {
            return point::min( ).vPos[ uiDim ];
        };

        coordinate_t max_value( ) const
        {
            return point::max( ).vPos[ uiDim ];
        };
    };

    using sort_points_t = typename points_vec::sort_func_t<points_it_t::iterator, typeof( PointsComperator )>;

    using pos_t = typename type_defs::pos_t;
    using desc_vec = Desc<type_desfs>;

    points_vec vData;

  public:
    Points( std::string sPrefix ) : vData( points_vec_generator( sPrefix + ".points" ) )
    {}

    size_t add( pos_t vPos, std::string sDesc, desc_vec& vDesc )
    {
        vData.push_back( Point( vPos, vDesc.add( sDesc ) ) );
    }

    void forRange( std::function<void( const point& )> fDo, size_t uiFrom, size_t uiTo ) const
    {
        assert( uiTo >= uiFrom );

        points_it_t itEnd = vData.begin( ) + uiTo;
        for( points_it_t cIter = vData.begin( ) + uiFrom; cIter != itEnd; cIter++ )
            fDo( *cIter );
    }

    void sortByDim( size_t uiDim, size_t uiFrom, size_t uiTo )
    {
        sort_points_t( vData.begin( ) + uiFrom, vData.begin( ) + uiTo, PointsComperator( uiDim ) );
    }

    size_t size( )
    {
        return vData.size( );
    }
};

} // namespace bps