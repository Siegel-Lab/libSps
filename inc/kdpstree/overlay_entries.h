#pragma once

#include "kdpstree/type_defs.h"
#include <functional>
#include <string>
#include <cassert>
#include <iostream>


namespace kdpstree
{

template <typename type_defs> class OverlayEntries
{
    EXTRACT_TYPE_DEFS; // macro call

    vec_generator_t<overlay_entry_t> overlay_entry_vec_generator = vec_generator_t<overlay_entry_t>( );
    using overlay_entry_vec_t = typeof( overlay_entry_vec_generator( "" ) );


    /**
     * @brief see https://algorithmica.org/en/eytzinger
     *
     * @param uiPos
     * @param uiBegin
     * @param uiSize
     * @return size_t
     */
    size_t getIndex( coordinate_t uiPos, size_t uiBegin, size_t uiSize ) const
    {
        // prevent write I/O
        const overlay_entry_vec_t& vData = this->vData;

        size_t uiK = 1;
        size_t uiLastRight = 1;
        if constexpr( EXPLAIN_QUERY )
            std::cerr << "\t\t\tbin search for " << uiPos << " between " << uiBegin << " and "
                      << uiBegin + uiSize << std::endl;
        while( uiK <= uiSize )
        {
            if( vData[ uiK - 1 + uiBegin ].first == uiPos )
            {
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\t\t\tbin search found at " << uiK - 1 + uiBegin << std::endl;
                return uiK - 1 + uiBegin;
            }
            else if( vData[ uiK - 1 + uiBegin ].first < uiPos )
            {
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\t\t\tbin search going right due to val " << vData[ uiK - 1 + uiBegin ].first 
                              << " at " << uiK - 1 + uiBegin << std::endl;
                uiLastRight = uiK;
                uiK = 2 * uiK + 1;
            }
            else
            {
                if constexpr( EXPLAIN_QUERY )
                    std::cerr << "\t\t\tbin search going left due to val " << vData[ uiK - 1 + uiBegin ].first 
                              << " at " << uiK - 1 + uiBegin << std::endl;
                uiK = 2 * uiK;
            }
        }
        return uiLastRight - 1 + uiBegin;
    }

    void forRange( std::function<void( coordinate_t, const val_t& )>& fDo, size_t uiK, size_t uiBegin, size_t uiSize,
                   coordinate_t uiFrom, coordinate_t uiTo ) const
    {
        if( uiK <= uiSize )
        {
            if( vData[ uiK - 1 + uiBegin ].first >= uiFrom )
                forRange( fDo, 2 * uiK, uiBegin, uiSize );
            if( vData[ uiK - 1 + uiBegin ].first >= uiFrom && vData[ uiK - 1 + uiBegin ].first < uiTo )
                fDo( vData[ uiK - 1 + uiBegin ].first, vData[ uiK - 1 + uiBegin ].second );
            if( vData[ uiK - 1 + uiBegin ].first < uiTo )
                forRange( fDo, 2 * uiK + 1, uiBegin, uiSize );
        }
    }

    const val_t zero = 0;

  public:
    overlay_entry_vec_t vData;

    OverlayEntries( std::string sPrefix ) : vData( overlay_entry_vec_generator( sPrefix + ".overlay_entries" ) )
    {}

    void incSize( size_t uiNum )
    {
        vData.resize( vData.size( ) + uiNum );
    }

    const val_t& get( coordinate_t uiPos, size_t uiBegin, size_t uiSize ) const
    {
        // prevent write I/O
        const overlay_entry_vec_t& vData = this->vData;
        size_t uiIdx = getIndex( uiPos, uiBegin, uiSize );
        if( vData[uiIdx].first > uiPos )
        {
            if constexpr( EXPLAIN_QUERY )
                std::cerr << "\t\t\treturning zero from overlay since " << uiPos << " is to the bottom of all entries." << std::endl;
            return zero;
        }
        return vData[ uiIdx ].second;
    }

    bool has( coordinate_t uiPos, size_t uiBegin, size_t uiSize ) const
    {
        // prevent write I/O
        const overlay_entry_vec_t& vData = this->vData;
        return vData[ getIndex( uiPos, uiBegin, uiSize ) ].first == uiPos;
    }

    /**
     * @brief allows editing access to valies
     *
     * @param uiPos
     * @param uiBegin
     * @param uiSize
     * @return val_t&
     */
    val_t& variableGet( coordinate_t uiPos, size_t uiBegin, size_t uiSize )
    {
        return vData[ getIndex( uiPos, uiBegin, uiSize ) ].second;
    }

    void forRange( std::function<void( coordinate_t, const val_t& )>& fDo, size_t uiBegin, 
                   size_t uiSize, coordinate_t uiFrom, coordinate_t uiTo ) const
    {
        forRange( fDo, 1, uiBegin, uiSize, uiFrom, uiTo );
    }

    void forRange( std::function<void( coordinate_t, const val_t& )>& fDo, size_t uiBegin, size_t uiSize ) const
    {
        forRange( fDo, uiBegin, uiSize, 0, std::numeric_limits<coordinate_t>::max() );
    }

    size_t size( ) const
    {
        return vData.size( );
    }
    
    std::string print( ) const
    {
        std::string sRet = "";
        size_t uiI = 0;
        for(const auto& rX : vData)
            sRet += std::to_string(uiI++) + ": " + std::to_string(rX.first) + "->" + std::to_string(rX.second) + "\n";
        return sRet;
    }

    void clear()
    {
        vData.clear();
    }
};

} // namespace kdpstree