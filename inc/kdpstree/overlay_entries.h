#pragma once

#include "kdpstree/type_defs.h"
#include <functional>
#include <string>


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
        while( uiK <= uiSize )
        {
            if( vData[ uiK - 1 + uiBegin ].first == uiPos )
                return uiK - 1 + uiBegin;
            else if( vData[ uiK - 1 + uiBegin ].first < uiPos )
                uiK = 2 * uiK + 1;
            else
            {
                uiLastRight = uiK;
                uiK = 2 * uiK;
            }
        }
        return uiLastRight - 1 + uiBegin;
    }

    size_t forRange( std::function<void( size_t, val_t& )>& fDo, size_t uiI, size_t uiK, size_t uiBegin, size_t uiSize )
    {
        if( uiK <= uiSize )
        {
            uiI = forRange( fDo, uiI, 2 * uiK + 1, uiBegin, uiSize );
            fDo( uiI++, vData[ uiK - 1 + uiBegin ].second );
            uiI = forRange( fDo, uiI, 2 * uiK, uiBegin, uiSize );
        }
        return uiI;
    }

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
        return vData[ getIndex( uiPos, uiBegin, uiSize ) ].second;
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

    void forRange( std::function<void( size_t, val_t& )>& fDo, size_t uiBegin, size_t uiSize )
    {
        forRange( fDo, 0, 1, uiBegin, uiSize );
    }

    size_t size( ) const
    {
        return vData.size( );
    }
};

} // namespace kdpstree