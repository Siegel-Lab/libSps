#pragma once

#include "sps/type_defs.h"
#include <iostream>
#include <string>

namespace sps
{

template <typename type_defs> class Desc;

}

namespace std
{

template <typename type_defs> ostream& operator<<( ostream& os, const typename sps::Desc<type_defs>& rDesc );

} // namespace std

namespace sps
{

template <typename type_defs> class Desc
{
    EXTRACT_TYPE_DEFS; // macro call

    EXTRACT_VEC_GENERATOR( desc, char ); // macro call

    desc_file_t xFile;
    desc_vec_t vData;
    char cEof;

    friend std::ostream& std::operator<< <>( std::ostream& os, const Desc& rTree );

  public:
    Desc( std::string sPrefix, bool bWrite, char cEof )
        : xFile( desc_vec_generator.file( sPrefix + ".desc", bWrite ) ),
          vData( desc_vec_generator.vec( xFile ) ),
          cEof( cEof )
    {}
    Desc( std::string sPrefix, bool bWrite ) : Desc( sPrefix, bWrite, std::char_traits<char>::eof( ) )
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

    void clear( )
    {
        vData.clear( );
    }

    friend std::ostream& operator<<( std::ostream& os, const Desc& rDesc )
    {
        os << "0: ";
        size_t uiI = 0;
        auto cIter = rDesc.vData.begin( );
        while( cIter != rDesc.vData.end( ) )
        {
            uiI++;
            if( *cIter == rDesc.cEof )
                os << std::endl << uiI << ": ";
            else
                os << *cIter;
            ++cIter;
        }
        os << "<EoF>" << std::endl;
        return os;
    }
};

} // namespace sps
