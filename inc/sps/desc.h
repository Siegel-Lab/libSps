#pragma once

#include "sps/type_defs.h"
#include <iostream>
#include <string>

namespace sps
{

template <template <typename> typename vec_gen_t> class DescImpl;

}

namespace std
{

template <template <typename> typename vec_gen_t>
ostream& operator<<( ostream& os, const typename sps::DescImpl<vec_gen_t>& rDesc );

} // namespace std

namespace sps
{
template <template <typename> typename vec_gen_t> class DescImpl
{
    using char_gen_t = vec_gen_t<char>;
    char_gen_t desc_vec_generator = char_gen_t( );

    typename char_gen_t::file_t xFile;
    typename char_gen_t::vec_t vData;
    char cEof;

    friend std::ostream& std::operator<< <>( std::ostream& os, const DescImpl& rTree );

  public:
    DescImpl( std::string sPrefix, bool bWrite, char cEof )
        : xFile( desc_vec_generator.file( sPrefix + ".desc", bWrite ) ),
          vData( desc_vec_generator.vec( xFile ) ),
          cEof( cEof )
    {}
    DescImpl( std::string sPrefix, bool bWrite ) : DescImpl( sPrefix, bWrite, std::char_traits<char>::eof( ) )
    {}

    DescImpl( )
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
        typename char_gen_t::vec_t::const_iterator cIter = vData.begin( ) + uiPos;
        while( *cIter != cEof )
            sRet += *( cIter++ );
        return sRet;
    }

    void clear( )
    {
        vData.clear( );
    }

    friend std::ostream& operator<<( std::ostream& os, const DescImpl& rDesc )
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

template <typename type_defs> using Desc = DescImpl<type_defs::template desc_vec_generator_t>;


} // namespace sps
