#pragma once

#include "kdpstree/util.h"

#include "kdpstree/type_defs.h"
#include <cassert>
#include <functional>
#include <iostream>
#include <string>


namespace kdpstree
{

template <typename type_defs> class OverlayEntries
{
    EXTRACT_TYPE_DEFS; // macro call

    EXTRACT_VEC_GENERATOR( entry, overlay_entry_t ); // macro call


    /**
     * @brief see https://algorithmica.org/en/eytzinger
     *
     * @param uiPos
     * @param uiBegin
     * @param uiSize
     * @return size_t
     */
    template <bool SILENT> size_t getIndex( coordinate_t uiPos, size_t uiBegin, size_t uiSize ) const
    {
        // prevent write I/O
        const entry_vec_t& vData = this->vData;

        size_t uiK = 1;
        size_t uiLastRight = 1;
        if constexpr( EXPLAIN_QUERY && !SILENT )
            std::cerr << "\t\t\tbin search for " << uiPos << " between " << uiBegin << " and " << uiBegin + uiSize
                      << std::endl;
        while( uiK <= uiSize )
        {
            if( vData[ uiK - 1 + uiBegin ].first == uiPos )
            {
                if constexpr( EXPLAIN_QUERY && !SILENT )
                    std::cerr << "\t\t\tbin search found at " << uiK - 1 + uiBegin << std::endl;
                return uiK - 1 + uiBegin;
            }
            else if( vData[ uiK - 1 + uiBegin ].first < uiPos )
            {
                if constexpr( EXPLAIN_QUERY && !SILENT )
                    std::cerr << "\t\t\tbin search going right due to val " << vData[ uiK - 1 + uiBegin ].first
                              << " at " << uiK - 1 + uiBegin << std::endl;
                uiLastRight = uiK;
                uiK = 2 * uiK + 1;
            }
            else
            {
                if constexpr( EXPLAIN_QUERY && !SILENT )
                    std::cerr << "\t\t\tbin search going left due to val " << vData[ uiK - 1 + uiBegin ].first << " at "
                              << uiK - 1 + uiBegin << std::endl;
                uiK = 2 * uiK;
            }
        }
        return uiLastRight - 1 + uiBegin;
    }

    void forRange( std::function<void( coordinate_t, const data_t& )>& fDo, size_t uiK, size_t uiBegin, size_t uiSize,
                   coordinate_t uiFrom, coordinate_t uiTo, bool bHaveGoneLeft ) const
    {
        if( uiK <= uiSize )
        {
            if( vData[ uiK - 1 + uiBegin ].first >= uiFrom )
                forRange( fDo, 2 * uiK, uiBegin, uiSize, uiFrom, uiTo, true );
            bool bLeftOk = vData[ uiK - 1 + uiBegin ].first < uiTo;
            bool bRightOk = vData[ uiK - 1 + uiBegin ].first >= uiFrom || 2 * uiK > uiSize;
            if( bLeftOk && bRightOk )
                fDo( vData[ uiK - 1 + uiBegin ].first, vData[ uiK - 1 + uiBegin ].second );
            if( vData[ uiK - 1 + uiBegin ].first < uiTo )
                forRange( fDo, 2 * uiK + 1, uiBegin, uiSize, uiFrom, uiTo, bHaveGoneLeft );
        }
    }

    const data_t zero{ };

  public:
    entry_file_t xFile;
    entry_vec_t vData;

    OverlayEntries( std::string sPrefix )
        : xFile( entry_vec_generator.file( sPrefix + ".overlay_entries" ) ), vData( entry_vec_generator.vec( xFile ) )
    {}

    void incSize( size_t uiNum )
    {
        vData.resize( vData.size( ) + uiNum );
    }

    const data_t& get( coordinate_t uiPos, size_t uiBegin, size_t uiSize ) const
    {
        // prevent write I/O
        const entry_vec_t& vData = this->vData;
        size_t uiIdx = getIndex<true>( uiPos, uiBegin, uiSize );
        if( vData[ uiIdx ].first > uiPos )
        {
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "\t\t\treturning zero from overlay since " << uiPos << " is to the bottom of all entries."
                          << std::endl;
            return zero;
        }
        return vData[ uiIdx ].second;
    }

    const data_t& getLast( size_t uiBegin, size_t uiSize ) const
    {
        return get( std::numeric_limits<coordinate_t>::max( ), uiBegin, uiSize );
    }

    bool has( coordinate_t uiPos, size_t uiBegin, size_t uiSize ) const
    {
        // prevent write I/O
        const entry_vec_t& vData = this->vData;
        return vData[ getIndex<true>( uiPos, uiBegin, uiSize ) ].first == uiPos;
    }

    /**
     * @brief allows editing access to valies
     *
     * @param uiPos
     * @param uiBegin
     * @param uiSize
     * @return data_t&
     */
    data_t& variableGet( coordinate_t uiPos, size_t uiBegin, size_t uiSize )
    {
        return vData[ getIndex<true>( uiPos, uiBegin, uiSize ) ].second;
    }

    void forRange( std::function<void( coordinate_t, const data_t& )> fDo, size_t uiBegin, size_t uiSize,
                   coordinate_t uiFrom, coordinate_t uiTo ) const
    {
        forRange( fDo, 1, uiBegin, uiSize, uiFrom, uiTo, false );
    }

    void forRange( std::function<void( coordinate_t, const data_t& )> fDo, size_t uiBegin, size_t uiSize ) const
    {
        forRange( fDo, uiBegin, uiSize, 0, std::numeric_limits<coordinate_t>::max( ) );
    }

    size_t size( ) const
    {
        return vData.size( );
    }

    void clear( )
    {
        vData.clear( );
    }
};


} // namespace kdpstree

namespace std
{

template <typename type_defs> ostream& operator<<( ostream& os, const kdpstree::OverlayEntries<type_defs>& rEntries )
{
    size_t uiI = 0;
    for( const typename type_defs::overlay_entry_t& rX : rEntries.vData )
        os << uiI++ << ": " << rX << endl;
    return os;
}

} // namespace std