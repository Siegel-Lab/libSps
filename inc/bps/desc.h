#pragma once

#include "bps/type_defs.h"
#include <string>

namespace bps
{

template <typename type_defs> class Desc
{
    using vec_generator = typename type_defs::vec_generator;
    using desc_vec = typeof( vec_generator<char>( )( "" ) );

    desc_vec vData;
    char cEof;

  public:
    Desc( std::string sPrefix, char cEof ) : vData( vec_generator<char>( )( sPrefix + ".desc" ) ), cEof( cEof )
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
        typename desc_vec::const_iterator cIter = vData.begin( ) + uiPos;
        while( *cIter != cEof )
            sRet += *( cIter++ );
        return sRet;
    }
};

} // namespace bps