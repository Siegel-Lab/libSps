#pragma once

#include "kdpstree/desc.h"
#include "kdpstree/point.h"
#include "kdpstree/type_defs.h"
#include <cassert>
#include <functional>
#include <string>



namespace kdpstree
{

template <typename type_defs> class Points
{
    EXTRACT_TYPE_DEFS; // macro call

    using point_t = Point<type_defs>;

    vec_generator_t<point_t> points_vec_generator = vec_generator_t<point_t>( );
    using points_vec_t = typeof( points_vec_generator( "" ) );
    using points_it_t = typename points_vec_t::iterator;
    using const_points_it_t = typename points_vec_t::const_iterator;

    struct PointsComperator
    {
        const size_t uiDim;

        PointsComperator( size_t uiDim ) : uiDim( uiDim )
        {}

        bool operator( )( const point_t& a, const point_t& b ) const
        {
            return a.vPos[ uiDim ] < b.vPos[ uiDim ];
        }

        coordinate_t min_value( ) const
        {
            return 0;
        };

        coordinate_t max_value( ) const
        {
            return std::numeric_limits<coordinate_t>::max( );
        };
    };

    sort_func_t<points_it_t, PointsComperator> sort_points = sort_func_t<points_it_t, PointsComperator>( );

  public:
    points_vec_t vData;

    Points( std::string sPrefix ) : vData( points_vec_generator( sPrefix + ".points" ) )
    {}

    size_t add( pos_t vPos, size_t uiDescOffset )
    {
        vData.push_back( point_t( vPos, uiDescOffset ) );
    }

    void forRange( std::function<bool( const point_t& )> fDo, size_t uiFrom, size_t uiTo ) const
    {
        assert( uiTo >= uiFrom );

        const_points_it_t itEnd = vData.begin( ) + uiTo;
        bool bContinue = true;
        for( const_points_it_t cIter = vData.begin( ) + uiFrom; cIter != itEnd && bContinue; cIter++ )
            bContinue = fDo( *cIter );
    }

    void sortByDim( size_t uiDim, size_t uiFrom, size_t uiTo )
    {
        sort_points( vData.begin( ) + uiFrom, vData.begin( ) + uiTo, PointsComperator( uiDim ) );
    }

    size_t size( ) const
    {
        return vData.size( );
    }
    
    std::string print( ) const
    {
        std::string sRet = "";
        size_t uiX = 0;
        for(const auto& rP : vData)
        {
            sRet += std::to_string(uiX++) + ": ( ";
            for( size_t uiI = 0; uiI < d; uiI++ )
                sRet += std::to_string(rP.vPos[uiI]) + ", ";
            sRet += ") : d" + std::to_string(rP.uiDescOffset) + "\n";
        }
        return sRet;
    }

    void clear()
    {
        vData.clear();
    }
};

} // namespace kdpstree