#pragma once

#include "kdpstree/type_defs.h"
#include <cassert>
#include <functional>
#include <string>


namespace kdpstree
{

template <typename type_defs> class SparseCoord
{
    EXTRACT_TYPE_DEFS; // macro call

    using entry_t = std::pair<coordinate_t, coordinate_t>;

    EXTRACT_VEC_GENERATOR( coord, entry_t ); // macro call

  public:
    coord_file_t xFile;
    coord_vec_t vData;

    struct Entry{
        coordinate_t uiStartIndex;
        coordinate_t uiStartCord;
        coordinate_t uiEndCord;
    };

    SparseCoord( std::string sPrefix, bool bWrite )
        : xFile( coord_vec_generator.file( sPrefix + ".coords", bWrite ) ), vData( coord_vec_generator.vec( xFile ) )
    {}

    coordinate_t replace(coordinate_t uiX, const Entry& rInfo) const
    {
        if(uiX < rInfo.uiStartCord)
            return std::numeric_limits<coordinate_t>::max();
        if(uiX >= rInfo.uiEndCord)
            return vData[rInfo.uiStartIndex + rInfo.uiEndCord - 1 - rInfo.uiStartCord].second;
        return vData[rInfo.uiStartIndex + uiX - rInfo.uiStartCord].second;
    }
    
    template<size_t N>
    std::array<coordinate_t, N> axisSizes(const std::array<Entry, N>& vAxes) const
    {
        std::array<coordinate_t, N> vAxisSizes;
        for(size_t uiI = 0; uiI < N; uiI++)
            vAxisSizes[uiI] = replace(vAxes[uiI].uiEndCord, vAxes[uiI]) + 1;
        return vAxisSizes;
    }

    template<size_t N>
    std::array<coordinate_t, N> sparse(
                                        const std::array<coordinate_t, N>& vCoords,
                                        const std::array<Entry, N>& vAxes) const
    {
        std::array<coordinate_t, N> vRet;
        for(size_t uiI = 0; uiI < N; uiI++)
            vRet[uiI] = replace(vCoords[uiI], vAxes[uiI]);
        return vRet;
    }

    template<typename Iterator_t>
    Entry add( Iterator_t xBegin, const Iterator_t& xEnd )
    {
        Entry xRet{};
        xRet.uiStartIndex = vData.size();
        xRet.uiStartCord = *xBegin;
        auto uiLast = *xBegin;
        size_t uiI = 0;
        ++xBegin;
        while(xBegin != xEnd)
        {
            for(coordinate_t uiX = uiLast; uiX < *xBegin; uiX++)
                vData.push_back(entry_t(uiX, uiI));
            uiLast = *xBegin;
            ++xBegin;
            uiI++;
        }
        xRet.uiEndCord = uiLast;
        vData.push_back(entry_t(uiLast, uiI));

        return xRet;
    }

    class EntryIterator
    {
        const SparseCoord& rCord;
        const Entry& rInfo;
        size_t uiI;
    public:
        EntryIterator(const SparseCoord& rCord, const Entry& rInfo) :
            rCord(rCord),
            rInfo(rInfo),
            uiI(0)
        {}

        void operator++()
        {
            coordinate_t uiLast = (**this).second;
            while(uiI < rInfo.uiEndCord - rInfo.uiStartCord && (**this).second == uiLast)
                uiI++;
        }

        const std::pair<coordinate_t, coordinate_t>& operator*() const
        {
            return rCord.vData[uiI + rInfo.uiStartIndex];
        }

        const std::pair<coordinate_t, coordinate_t>* operator->() const
        {
            return &rCord.vData[uiI + rInfo.uiStartIndex];
        }

        bool operator!=(const EntryIterator& rOther) const
        {
            return uiI != rOther.uiI;
        }

        friend class SparseCoord;
    };

    EntryIterator cbegin(const Entry& rInfo) const
    {
        return EntryIterator(*this, rInfo);
    }

    EntryIterator cend(const Entry& rInfo) const
    {
        EntryIterator xRet(*this, rInfo);
        xRet.uiI += rInfo.uiEndCord - rInfo.uiStartCord;
        return xRet;
    }

    void iterate(std::function<void(coordinate_t, coordinate_t)> fDo, const Entry& rInfo) const
    {
        auto xIt = this->cbegin(rInfo);
        while(xIt != this->cend(rInfo))
            fDo(xIt->first, xIt->second);
    }

    template<size_t I, size_t N>
    inline typename std::enable_if<I != N, void>::type
    iterateHelper(
                 std::function<void(const std::array<coordinate_t, N>&, const std::array<coordinate_t, N>&)> fDo,
                 const std::array<Entry, N>& rInfos, 
                 std::array<coordinate_t, N>& rFrom, std::array<coordinate_t, N>& rTo) const
    {
        iterate(
            [&](coordinate_t uiFrom, coordinate_t uiTo){
                rFrom[I] = uiFrom;
                rTo[I] = uiTo;
                iterateHelper<I + 1, N>(fDo, rInfos, rFrom, rTo);
            },
            rInfos[I]
        );
    }

    template<size_t I, size_t N>
    inline typename std::enable_if<I == N, void>::type
    iterateHelper(
                 std::function<void(const std::array<coordinate_t, N>&, const std::array<coordinate_t, N>&)> fDo,
                 const std::array<Entry, N>& , std::array<coordinate_t, N>& rFrom,
                 std::array<coordinate_t, N>& rTo) const
    {
        fDo(rFrom, rTo);
    }

    template<size_t N>
    void iterate(
                 std::function<void(const std::array<coordinate_t, N>&, const std::array<coordinate_t, N>&)> fDo,
                 const std::array<Entry, N>& rInfos
                 ) const
    {
        std::array<coordinate_t, N> rFrom;
        std::array<coordinate_t, N> rTo;
        iterateHelper<0, N>(
            fDo,
            rInfos,
            rFrom,
            rTo
        );
    }

    void clear( )
    {
        vData.clear( );
    }
};


} // namespace kdpstree
