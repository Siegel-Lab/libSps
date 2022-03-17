#pragma once

#include "kdpstree/util.h"

#include "kdpstree/type_defs.h"
#include "kdpstree/desc.h"


#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif


namespace kdpstree
{

template <typename type_defs> class PsArray
{
    EXTRACT_TYPE_DEFS; // macro call

  public:
    struct Entry
    {
        class_key_t uiClass;
        layers_t uiLayer;

        coordinate_t uiPos;
        size_t uiDescOffset;

        Entry( class_key_t uiClass, layers_t uiLayer, coordinate_t uiPos, size_t uiDescOffset )
            : uiClass( uiClass ), uiLayer( uiLayer ), uiPos( uiPos ), uiDescOffset( uiDescOffset )
        {}

        Entry( ) : uiClass( 0 ), uiLayer( 0 ), uiPos( 0 ), uiDescOffset( 0 )
        {}

        friend std::ostream& operator<<( std::ostream& os, const Entry& rEntry )
        {
            os << "c" << rEntry.uiClass << ", l" << (size_t)rEntry.uiLayer << ", p" << rEntry.uiPos << ", d"
               << rEntry.uiDescOffset;
            return os;
        }
    };

    EXTRACT_VEC_GENERATOR( entry, Entry ); // macro call

    Desc<type_defs> vDesc;
    entry_file_t xFile;
    entry_vec_t vEntries;
    std::vector<std::array<std::pair<offset_t, offset_t>, LAYERS>> vClassLayerCache;

    struct EntryComperator
    {
        bool operator( )( const Entry& a, const Entry& b ) const
        {
            if( a.uiClass == b.uiClass )
            {
                if( a.uiLayer == b.uiLayer )
                    return a.uiPos < b.uiPos;
                return a.uiLayer < b.uiLayer;
            }
            return a.uiClass < b.uiClass;
        }

        Entry min_value( ) const
        {
            return Entry{ };
        };

        Entry max_value( ) const
        {
            Entry xRet{ };
            xRet.uiClass = std::numeric_limits<class_key_t>::max( );
            return xRet;
        };
    };

    typename entry_vec_t::const_iterator iteratorFor( class_key_t uiClass, layers_t uiLayer, coordinate_t uiPos,
                                                      offset_t uiStart, offset_t uiEnd ) const
    {
        auto xRet = std::lower_bound( vEntries.begin( ) + uiStart,
                                      uiEnd >= vEntries.size( ) ? vEntries.end( ) : vEntries.begin( ) + uiEnd,
                                      Entry( uiClass, uiLayer, uiPos, 0 ), EntryComperator( ) );

        if constexpr( EXPLAIN_QUERY )
            std::cerr << "iteratorFor: iterator pos: " << xRet - vEntries.begin( ) << std::endl;
        return xRet;
    }

    val_t countToZero( class_key_t uiClass, layers_t uiLayer, coordinate_t uiPos, offset_t uiStart,
                       offset_t uiEnd ) const
    {
        auto xIt = iteratorFor( uiClass, uiLayer, uiPos, uiStart, uiEnd );
        if( xIt == vEntries.end( ) )
            return vEntries.size( ) - uiStart;
        return (val_t)( xIt - vEntries.begin( ) ) - (val_t)uiStart;
    }

    val_t countHelper( class_key_t uiClass, layers_t uiLayer, coordinate_t uiFrom, coordinate_t uiTo, offset_t uiStart,
                       offset_t uiEnd ) const
    {
        return countToZero( uiClass, uiLayer, uiTo, uiStart, uiEnd ) -
               countToZero( uiClass, uiLayer, uiFrom, uiStart, uiEnd );
    }

    void genClassLayerCache( )
    {
        vClassLayerCache.clear( );
        if( vEntries.size( ) > 0 )
        {
            class_key_t uiLastClass = vEntries.back( ).uiClass;
            offset_t uiSum = 0;
            for( class_key_t uiI = 0; uiI <= uiLastClass; uiI++ )
            {
                vClassLayerCache.emplace_back( );
                for( class_key_t uiJ = 0; uiJ < LAYERS; uiJ++ )
                {
                    vClassLayerCache.back( )[ uiJ ].first = uiSum;
                    uiSum += countToZero( uiI, uiJ, std::numeric_limits<coordinate_t>::max( ), uiSum,
                                          vClassLayerCache.size( ) );
                    vClassLayerCache.back( )[ uiJ ].second = uiSum;
                }
            }
        }
    }

    PsArray( std::string sPrefix )
        : vDesc( sPrefix ),
          xFile( entry_vec_generator.file( sPrefix + ".entries" ) ),
          vEntries( entry_vec_generator.vec( xFile ) ),
          vClassLayerCache{ }
    {
        genClassLayerCache( );
    }

    void addPoint( class_key_t uiClass, coordinate_t uiPos, layers_t uiLayer, std::string sDesc )
    {
        vEntries.push_back( Entry( uiClass, uiLayer, uiPos, vDesc.add( sDesc ) ) );
    }

    void addPointArray( class_key_t uiClass, std::array<coordinate_t, 1> vPos, layers_t uiLayer, std::string sDesc )
    {
        addPoint( uiClass, vPos[ 0 ], uiLayer, sDesc );
    }

    sort_func_t<typename entry_vec_t::iterator, EntryComperator> sort_points =
        sort_func_t<typename entry_vec_t::iterator, EntryComperator>( );

    void generate( )
    {
        sort_points( vEntries.begin( ), vEntries.end( ), EntryComperator( ) );
        genClassLayerCache( );
    }

    void generate2( class_key_t, size_t, size_t, size_t )
    {
        generate( );
    }

    val_t countOne( class_key_t uiClass, layers_t uiLayer, coordinate_t uiFrom, coordinate_t uiTo ) const
    {
        return countHelper( uiClass, uiLayer, uiFrom, uiTo, vClassLayerCache[ uiClass ][ uiLayer ].first,
                            vClassLayerCache[ uiClass ][ uiLayer ].second );
    }

    std::array<val_t, LAYERS> count( class_key_t uiClass, coordinate_t uiFrom, coordinate_t uiTo ) const
    {
        std::array<val_t, LAYERS> vRet{ };
        for( class_key_t uiJ = 0; uiJ < LAYERS; uiJ++ )
            vRet[ uiJ ] = countOne( uiClass, uiJ, uiFrom, uiTo );
        return vRet;
    }

    std::array<val_t, LAYERS> countArray( class_key_t uiClass, std::array<coordinate_t, 1> vFrom,
                                          std::array<coordinate_t, 1> vTo ) const
    {
        return count( uiClass, vFrom[ 0 ], vTo[ 0 ] );
    }

    void clear( )
    {
        vClassLayerCache.clear( );
        vEntries.clear( );
        vDesc.clear( );
    }

    std::array<std::vector<std::pair<coordinate_t, std::string>>, LAYERS> get( class_key_t uiClass, coordinate_t uiFrom,
                                                                               coordinate_t uiTo ) const
    {
        std::array<std::vector<std::pair<coordinate_t, std::string>>, LAYERS> vRet{ };
        for( class_key_t uiJ = 0; uiJ < LAYERS; uiJ++ )
        {
            auto xIt = iteratorFor( uiClass, uiJ, uiFrom, vClassLayerCache[ uiClass ][ uiJ ].first,
                                    vClassLayerCache[ uiClass ][ uiJ ].second );
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "get: iterator pos: " << xIt - vEntries.begin( ) << std::endl;
            while( xIt != vEntries.end( ) && xIt->uiClass == uiClass && xIt->uiLayer == uiJ && xIt->uiPos < uiTo )
            {
                vRet[ uiJ ].emplace_back( xIt->uiPos, vDesc.get( xIt->uiDescOffset ) );
                ++xIt;
            }
        }
        return vRet;
    }

    std::array<std::vector<std::pair<std::array<coordinate_t, 1>, std::string>>, LAYERS>
    getArray( class_key_t uiClass, std::array<coordinate_t, 1> vFrom, std::array<coordinate_t, 1> vTo ) const
    {
        std::array<std::vector<std::pair<std::array<coordinate_t, 1>, std::string>>, LAYERS> vRet{ };
        for( class_key_t uiJ = 0; uiJ < LAYERS; uiJ++ )
        {
            auto xIt = iteratorFor( uiClass, uiJ, vFrom[ 0 ], vClassLayerCache[ uiClass ][ uiJ ].first,
                                    vClassLayerCache[ uiClass ][ uiJ ].second );

            if constexpr( EXPLAIN_QUERY )
                std::cerr << "getArray: iterator pos: " << xIt - vEntries.begin( ) << " " << uiClass << " " << uiJ
                          << std::endl;
            while( xIt != vEntries.end( ) && xIt->uiClass == uiClass && xIt->uiLayer == uiJ && xIt->uiPos < vTo[ 0 ] )
            {
                vRet[ uiJ ].emplace_back( std::array<coordinate_t, 1>{ xIt->uiPos }, vDesc.get( xIt->uiDescOffset ) );
                ++xIt;
            }
        }
        return vRet;
    }

    friend std::ostream& operator<<( std::ostream& os, const PsArray& rArray )
    {
        os << "Cache:" << std::endl;
        for( size_t uiI = 0; uiI < rArray.vClassLayerCache.size( ); uiI++ )
            os << rArray.vClassLayerCache[ uiI ] << std::endl;

        os << "Array:" << std::endl;
        for( size_t uiI = 0; uiI < rArray.vEntries.size( ); uiI++ )
            os << rArray.vEntries[ uiI ] << std::endl;

        os << "Desc:" << std::endl << rArray.vDesc << std::endl;
        return os;
    }

    std::string print( )
    {
        std::stringstream ss;
        ss << *this << std::endl;
        return ss.str( );
    }
};


} // namespace kdpstree

#if WITH_PYTHON
template <typename type_defs> void exportArray( pybind11::module& m, std::string sName )
{
    pybind11::class_<kdpstree::PsArray<type_defs>>( m, sName.c_str( ) )
        .def( pybind11::init<std::string>( ) ) // constructor
        .def( "add_point", &kdpstree::PsArray<type_defs>::addPoint )
        .def( "add_point", &kdpstree::PsArray<type_defs>::addPointArray )
        .def( "generate", &kdpstree::PsArray<type_defs>::generate )
        .def( "generate", &kdpstree::PsArray<type_defs>::generate2 )
        .def( "get", &kdpstree::PsArray<type_defs>::get )
        .def( "get", &kdpstree::PsArray<type_defs>::getArray )
        .def( "count", &kdpstree::PsArray<type_defs>::count )
        .def( "count", &kdpstree::PsArray<type_defs>::countArray )
        .def( "clear", &kdpstree::PsArray<type_defs>::clear )
        .def( "__str__", &kdpstree::PsArray<type_defs>::print )

        ;
}
#endif