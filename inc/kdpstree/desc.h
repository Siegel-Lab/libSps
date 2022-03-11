#pragma once

#include "kdpstree/type_defs.h"
#include <string>

namespace kdpstree
{
    
template <typename type_defs> class Desc;

template <typename type_defs>
std::ostream& operator<<(std::ostream& os, const Desc<type_defs>& rDesc);

template <typename type_defs> class Desc
{
    EXTRACT_TYPE_DEFS; // macro call

    EXTRACT_VEC_GENERATOR(desc, char); // macro call

    desc_file_t xFile;
    desc_vec_t vData;
    char cEof;
    
    friend std::ostream& operator<< <>(std::ostream& os, const Desc& rTree);

  public:
    Desc( std::string sPrefix, char cEof ) : 
        xFile( desc_vec_generator.file(sPrefix + ".desc") ), 
        vData( desc_vec_generator.vec( xFile ) ), 
        cEof( cEof )
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

    void clear()
    {
        vData.clear();
    }
};

template <typename type_defs>
std::ostream& operator<<(std::ostream& os, const Desc<type_defs>& rDesc)
{
    os << "0: ";
    size_t uiI = 0;
    auto cIter = rDesc.vData.begin( );
    while( cIter != rDesc.vData.end() )
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


} // namespace kdpstree