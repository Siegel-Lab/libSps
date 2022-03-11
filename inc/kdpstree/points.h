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

    EXTRACT_VEC_GENERATOR(points, point_t); // macro call

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

        point_t min_value( ) const
        {
            return point_t();
        };

        point_t max_value( ) const
        {
            point_t xRet {};
            xRet.vPos[uiDim] = std::numeric_limits<coordinate_t>::max( );
            return xRet;
        };
    };

    struct PointsLayerComperator: public PointsComperator 
    {
        using PointsComperator::PointsComperator;

        bool operator( )( const point_t& a, const point_t& b ) const
        {
            if(a.uiLayer == b.uiLayer)
                return PointsComperator::operator()(a, b);
            return a.uiLayer < b.uiLayer;
        }

        point_t max_value( ) const
        {
            point_t xRet = PointsComperator::max_value();
            xRet.uiLayer = LAYERS;
            return xRet;
        };
    };

    sort_func_t<points_it_t, PointsComperator> sort_points = sort_func_t<points_it_t, PointsComperator>( );
    sort_func_t<points_it_t, PointsLayerComperator> sort_layer_points = sort_func_t<points_it_t, PointsLayerComperator>( );

  public:
    points_file_t xFile;
    points_vec_t vData;

    Points( std::string sPrefix ) : xFile(points_vec_generator.file( sPrefix + ".points" )), vData(points_vec_generator.vec( xFile ))
    {}

    size_t add( pos_t vPos, size_t uiDescOffset, layers_t uiLayer )
    {
        vData.push_back( point_t( vPos, uiDescOffset, uiLayer ) );
    }

    void forRange( std::function<bool( const point_t& )> fDo, offset_t uiFrom, offset_t uiTo ) const
    {
        assert( uiTo >= uiFrom );

        const_points_it_t itEnd = vData.begin( ) + uiTo;
        bool bContinue = true;
        for( const_points_it_t cIter = vData.begin( ) + uiFrom; cIter != itEnd && bContinue; cIter++ )
            bContinue = fDo( *cIter );
    }

    void sortByDim( size_t uiDim, offset_t uiFrom, offset_t uiTo )
    {
        sort_points( vData.begin( ) + uiFrom, vData.begin( ) + uiTo, PointsComperator( uiDim ) );
    }
    void sortByLayerAndDim( size_t uiDim, offset_t uiFrom, offset_t uiTo )
    {
        sort_layer_points( vData.begin( ) + uiFrom, vData.begin( ) + uiTo, PointsLayerComperator( uiDim ) );
    }


    size_t size( ) const
    {
        return vData.size( );
    }

    void clear()
    {
        vData.clear();
    }
};

template <typename type_defs>
std::ostream& operator<<(std::ostream& os, const Points<type_defs>& vPoints)
{
    size_t uiX = 0;
    for(const auto& rP : vPoints.vData)
        os << uiX++ << ": " << rP << std::endl;
    return os;
}

} // namespace kdpstree