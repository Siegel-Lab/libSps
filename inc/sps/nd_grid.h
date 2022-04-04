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
        coordinate_t uiStartIndex = std::numeric_limits<coordinate_t>::max();

        friend std::ostream& operator<<( std::ostream& os, const Entry& rEntry )
        {
            os << "s";
            os << rEntry.vAxisSizes;

            os << " i";
            os << rEntry.uiStartIndex;

            return os;
        }

        template<typename... TS>
        std::ostream& stream( std::ostream& os, const NDGrid& rGrid, TS&... args ) const
        {
            os << "{ "; 
            for(size_t uiI = 0; uiI < rGrid.sizeOf(*this); uiI++)
            {
                if(uiI > 0)
                    os << ", ";
                os << rGrid.posOf(uiI + uiStartIndex, *this) << ": ";
                rGrid.vData[uiI + uiStartIndex].stream(os, args...);
            }
            os << " }";

            return os;
        }

        std::ostream& streamOp( std::ostream& os, const NDGrid& rGrid ) const
        {
            os << "{ "; 
            for(size_t uiI = 0; uiI < rGrid.sizeOf(*this); uiI++)
            {
                if(uiI > 0)
                    os << ", ";
                os << rGrid.vData[uiI + uiStartIndex];
            }
            os << " }";

            return os;
        }
    };

    NDGrid( std::string sFileName, bool bWrite )
        : xFile( data_vec_generator.file( sFileName, bWrite ) ), vData( data_vec_generator.vec( xFile ) )
    {}

    template <size_t N>
    coordinate_t indexOf(const std::array<coordinate_t, N>& vX, const Entry<N>& rInfo) const
    {
        if(rInfo.uiStartIndex == std::numeric_limits<coordinate_t>::max())
            return std::numeric_limits<coordinate_t>::max();
        coordinate_t uiRet = 0;
        for(size_t uiI = 0; uiI < N; uiI++)
        {
            if(vX[uiI] == std::numeric_limits<coordinate_t>::max())
                return std::numeric_limits<coordinate_t>::max();
            assert(vX[uiI] < rInfo.vAxisSizes[uiI]);
            uiRet = uiRet * rInfo.vAxisSizes[uiI] + vX[uiI];
        }
        assert(uiRet < sizeOf(rInfo));
        return uiRet + rInfo.uiStartIndex;
    }

    template <size_t N>
    std::array<coordinate_t, N> posOf(coordinate_t uiIndexIn, const Entry<N>& rInfo) const
    {
        coordinate_t uiIndex = uiIndexIn - rInfo.uiStartIndex;
        std::array<coordinate_t, N> vRet;
        for(size_t _uiI = 0; _uiI < N; _uiI++)
        {
            size_t uiI = N - 1 - _uiI;
            vRet[uiI] = uiIndex % rInfo.vAxisSizes[uiI];
            uiIndex /= rInfo.vAxisSizes[uiI];
        }
        assert(uiIndex == 0);
        assert(indexOf(vRet, rInfo) == uiIndexIn);
        return vRet;
    }

    template <size_t N>
    coordinate_t sizeOf(const Entry<N>& rInfo) const
    {
        if(rInfo.uiStartIndex == std::numeric_limits<coordinate_t>::max())
            return 0;
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
        for(size_t uiI = 0; uiI < N; uiI++)
            xRet.vAxisSizes[uiI] = vAxisSizes[uiI];
        xRet.uiStartIndex = vData.size();
        // make space for the new grid
        vData.resize(xRet.uiStartIndex + sizeOf(xRet));
        // stxxl does not zero initialize the resized elements :'(
        for(size_t uiI = xRet.uiStartIndex; uiI < vData.size(); uiI++)
            vData[uiI] = data_t{};
        return xRet;
    }

    template <size_t N>
    void remove( const Entry<N>& xEntry )
    {
        assert(xEntry.uiStartIndex + sizeOf(xEntry) == vData.size());
        vData.resize(xEntry.uiStartIndex);
    }

    void clear( )
    {
        vData.clear( );
    }

    friend std::ostream& operator<<( std::ostream& os, const NDGrid& rGrid )
    {
        os << rGrid.vData;

        return os;
    }
};

template <typename type_defs, typename data_t>
const data_t NDGrid<type_defs, data_t>::uiZero{};


} // namespace sps
