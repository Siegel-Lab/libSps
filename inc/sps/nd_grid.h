#pragma once

#include "sps/type_defs.h"
#include <cassert>
#include <functional>
#include <string>


namespace sps
{

template <typename type_defs, typename data_t> class NDGrid
{
    EXTRACT_TYPE_DEFS; // macro call

    EXTRACT_VEC_GENERATOR( data, data_t ); // macro call

  public:
    data_file_t xFile;
    data_vec_t vData;
    static const data_t uiZero;

    template <size_t N>
    struct Entry{
        std::array<coordinate_t, N> vAxisSizes; 
        coordinate_t uiStartIndex;
    };

    NDGrid( std::string sFileName, bool bWrite )
        : xFile( data_vec_generator.file( sFileName, bWrite ) ), vData( data_vec_generator.vec( xFile ) )
    {}

    template <size_t N>
    coordinate_t indexOf(const std::array<coordinate_t, N>& vX, const Entry<N>& rInfo) const
    {
        coordinate_t uiRet = 0;
        for(size_t uiI = 0; uiI < N; uiI++)
        {
            if(vX[uiI] == std::numeric_limits<coordinate_t>::max())
                return std::numeric_limits<coordinate_t>::max();
            uiRet = uiRet * rInfo.vAxisSizes[uiI] + vX[uiI];
        }
        return uiRet + rInfo.uiStartIndex;
    }

    template <size_t N>
    std::array<coordinate_t, N> posOf(coordinate_t uiIndex, const Entry<N>& rInfo) const
    {
        uiIndex -= rInfo.uiStartIndex;
        std::array<coordinate_t, N> vRet;
        for(size_t _uiI = 0; _uiI < N; _uiI++)
        {
            size_t uiI = N - 1 - _uiI;
            vRet[uiI] += uiIndex % rInfo.vAxisSizes[uiI];
            uiIndex /= rInfo.vAxisSizes[uiI];
        }
        return vRet;
    }

    template <size_t N>
    coordinate_t sizeOf(const Entry<N>& rInfo) const
    {
        coordinate_t uiRet = 1;
        for(size_t uiI = 0; uiI < N; uiI++)
            uiRet *= rInfo.vAxisSizes[uiI];
        return uiRet;
    }

    template <size_t N>
    const data_t& get(const std::array<coordinate_t, N>& vX, const Entry<N>& rInfo) const
    {
        auto uiIdx = indexOf<N>(vX, rInfo);
        if(uiIdx == std::numeric_limits<coordinate_t>::max())
            return uiZero;
        return vData[uiIdx];
    }

    template <size_t N>
    data_t& get(const std::array<coordinate_t, N>& vX, const Entry<N>& rInfo)
    {
        auto uiIdx = indexOf<N>(vX, rInfo);
        assert(uiIdx != std::numeric_limits<coordinate_t>::max());
        return vData[uiIdx];
    }

    template <size_t N>
    Entry<N> add( const std::array<coordinate_t, N>& vAxisSizes )
    {
        Entry<N> xRet{};
        xRet.vAxisSizes = vAxisSizes;
        xRet.uiStartIndex = vData.size();
        // make space for the new grid
        vData.resize(indexOf<N>(vAxisSizes, xRet));
        return xRet;
    }

    template <size_t N>
    void remove( const Entry<N>& xEntry )
    {
        assert(xEntry.uiStartIndex + indexOf(xEntry.vAxisSizes, xEntry) == vData.size());
        vData.resize(xEntry.uiStartIndex);
    }

    void clear( )
    {
        vData.clear( );
    }
};

template <typename type_defs, typename data_t>
const data_t NDGrid<type_defs, data_t>::uiZero{};


} // namespace sps
