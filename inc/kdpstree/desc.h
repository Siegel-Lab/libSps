#pragma once

#include "kdpstree/type_defs.h"
#include <string>

namespace kdpstree
{

template <typename type_defs> class Desc
{
    EXTRACT_TYPE_DEFS; // macro call

    vec_generator_t<char> vec_generator = vec_generator_t<char>( );
    using desc_vec_t = typeof( vec_generator( "" ) );

    desc_vec_t vData;
    char cEof;

  public:
    Desc( std::string sPrefix, char cEof ) : vData( vec_generator( sPrefix + ".desc" ) ), cEof( cEof )
    {}
    Desc( std::string sPrefix ) : Desc( sPrefix, std::char_traits<char>::eof( ) )
    {}

    size_t add( std::string sDesc )
    {
        size_t uiRet = vData.size( );
        for( char c : sDesc )
            vData.push_back( c );
        vData.push_back( cEof );
        return uiRet;
    }

    std::string get( size_t uiPos ) const
    {
        std::string sRet = "";
        typename desc_vec_t::const_iterator cIter = vData.begin( ) + uiPos;
        while( *cIter != cEof )
            sRet += *( cIter++ );
        return sRet;
    }
};

} // namespace kdpstree