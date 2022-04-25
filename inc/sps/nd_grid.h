#pragma once

#include "sps/type_defs.h"
#include "sps/util.h"
#include <cassert>
#include <functional>
#include <string>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>


namespace sps
{

template <typename type_defs, typename data_t> class NDGrid
{
    EXTRACT_TYPE_DEFS; // macro call

    EXTRACT_VEC_GENERATOR( data, data_t ); // macro call

  public:
    static constexpr bool THREADSAVE = data_THREADSAVE;
    Lockable xLockable;
    data_file_t xFile;
    data_vec_t vData;
    static const data_t uiZero;

    template <size_t N> struct Entry
    {
        std::array<coordinate_t, N> vAxisSizes;
        coordinate_t uiStartIndex = std::numeric_limits<coordinate_t>::max( );

        friend std::ostream& operator<<( std::ostream& os, const Entry& rEntry )
        {
            os << "s";
            os << rEntry.vAxisSizes;

            os << " i";
            os << rEntry.uiStartIndex;

            return os;
        }

        template <typename... TS> std::ostream& stream( std::ostream& os, const NDGrid& rGrid, TS&... args ) const
        {
            os << "{ ";
            for( size_t uiI = 0; uiI < rGrid.sizeOf( *this ); uiI++ )
            {
                if( uiI > 0 )
                    os << ", ";
                os << rGrid.posOf( uiI + uiStartIndex, *this ) << ": ";
                rGrid.vData[ uiI + uiStartIndex ].stream( os, rGrid.posOf( uiI + uiStartIndex, *this ), args... );
            }
            os << " }";

            return os;
        }

        std::ostream& streamOp( std::ostream& os, const NDGrid& rGrid ) const
        {
            os << "{ ";
            for( size_t uiI = 0; uiI < rGrid.sizeOf( *this ); uiI++ )
            {
                if( uiI > 0 )
                    os << ", ";
                os << rGrid.vData[ uiI + uiStartIndex ];
            }
            os << " }";

            return os;
        }
    };

    NDGrid( std::string sFileName, bool bWrite )
        : xLockable(THREADSAVE ? std::numeric_limits<size_t>::max() : 1), 
          xFile( data_vec_generator.file( sFileName, bWrite ) ), vData( data_vec_generator.vec( xFile ) )
    {}

    template <size_t N> static coordinate_t indexOf( const std::array<coordinate_t, N>& vX, const Entry<N>& rInfo )
    {
        if( rInfo.uiStartIndex == std::numeric_limits<coordinate_t>::max( ) )
            return std::numeric_limits<coordinate_t>::max( );
        coordinate_t uiRet = 0;
        for( size_t uiI = 0; uiI < N; uiI++ )
        {
            if( vX[ uiI ] == std::numeric_limits<coordinate_t>::max( ) )
                return std::numeric_limits<coordinate_t>::max( );
            assert( vX[ uiI ] < rInfo.vAxisSizes[ uiI ] );
            uiRet = uiRet * rInfo.vAxisSizes[ uiI ] + vX[ uiI ];
        }
        assert( uiRet < sizeOf( rInfo ) );
        return uiRet + rInfo.uiStartIndex;
    }

    template <size_t N> static std::array<coordinate_t, N> posOf( coordinate_t uiIndexIn, const Entry<N>& rInfo )
    {
        coordinate_t uiIndex = uiIndexIn - rInfo.uiStartIndex;
        std::array<coordinate_t, N> vRet;
        for( size_t _uiI = 0; _uiI < N; _uiI++ )
        {
            size_t uiI = N - 1 - _uiI;
            vRet[ uiI ] = uiIndex % rInfo.vAxisSizes[ uiI ];
            uiIndex /= rInfo.vAxisSizes[ uiI ];
        }
        assert( uiIndex == 0 );
        assert( indexOf( vRet, rInfo ) == uiIndexIn );
        return vRet;
    }

    template <size_t N> static coordinate_t sizeOf( const Entry<N>& rInfo )
    {
        if( rInfo.uiStartIndex == std::numeric_limits<coordinate_t>::max( ) )
            return 0;
        coordinate_t uiRet = 1;
        for( size_t uiI = 0; uiI < N; uiI++ )
            uiRet *= rInfo.vAxisSizes[ uiI ];
        return uiRet;
    }

    template <size_t N> const data_t& get( const std::array<coordinate_t, N>& vX, const Entry<N>& rInfo ) const
    {
        auto uiIdx = indexOf<N>( vX, rInfo );
        if( uiIdx == std::numeric_limits<coordinate_t>::max( ) )
            return uiZero;
        return vData[ uiIdx ];
    }

    template <size_t N> data_t& get( const std::array<coordinate_t, N>& vX, const Entry<N>& rInfo )
    {
        auto uiIdx = indexOf<N>( vX, rInfo );
        assert( uiIdx != std::numeric_limits<coordinate_t>::max( ) );
        return vData[ uiIdx ];
    }

    template <size_t N> Entry<N> add( const std::array<coordinate_t, N>& vAxisSizes )
    {
        Entry<N> xRet{ };
        for( size_t uiI = 0; uiI < N; uiI++ )
            xRet.vAxisSizes[ uiI ] = vAxisSizes[ uiI ];
        xRet.uiStartIndex = vData.size( );
        // make space for the new grid

        // vData.resize(xRet.uiStartIndex + sizeOf(xRet));
        // stxxl does not zero initialize the resized elements :'(
        for( size_t uiI = 0; uiI < sizeOf( xRet ); uiI++ )
            vData.push_back( data_t{ } );
        return xRet;
    }

    template <size_t N> void remove( const Entry<N>& xEntry )
    {
        assert( xEntry.uiStartIndex + sizeOf( xEntry ) == vData.size( ) );
        vData.resize( xEntry.uiStartIndex );
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

    template <size_t N>
    class ParallelIterator
    {
        const Entry<N> xEntry;
        std::vector<int8_t> vNumPredComputed;
        std::queue<size_t> vNext;
        std::mutex xLock;
        std::condition_variable xVar;
        std::function<std::vector<size_t>(size_t, size_t,  const Entry<N>&)> fSuccessors;

        public:
            ParallelIterator(const Entry<N> xEntry, 
                             std::function<std::vector<size_t>(size_t, size_t, const Entry<N>&)> fSuccessors = 
                        [](size_t uiIdx, size_t uiDim, const Entry<N>& xEntry){
                            auto vPos = NDGrid::posOf(uiIdx, xEntry);
                            ++vPos[uiDim];
                            std::vector<size_t> vRet;
                            if(vPos[uiDim] < xEntry.vAxisSizes[uiDim])
                                vRet.push_back(NDGrid::indexOf(vPos, xEntry));
                            return vRet;
                        }) : 
                    xEntry(xEntry), vNumPredComputed(sizeOf(xEntry)), 
                    vNext(), xLock(), xVar(), fSuccessors(fSuccessors)
            {
                vNext.push(xEntry.uiStartIndex);
            }

            bool done()
            {
                std::unique_lock xGuard(xLock);
                return vNext.front() == std::numeric_limits<size_t>::max();
            }

            size_t getNext()
            {
                std::unique_lock xGuard(xLock);
                while(vNext.size() == 0)
                    xVar.wait(xGuard);
                assert(vNext.size() > 0);
                size_t uiRet = vNext.front();
                if(uiRet != std::numeric_limits<size_t>::max())
                    vNext.pop();
                return uiRet;
            }

            void doneWith(size_t uiIdx)
            {
                bool bHas = false;
                for(size_t uiD = 0; uiD < N; uiD++)
                    for(size_t uiIdx: fSuccessors(uiIdx, uiD, xEntry))
                    {
                        bHas = true;
                        size_t uiNumReq = 0;
                        auto vPos = NDGrid::posOf(uiIdx, xEntry);
                        for(size_t uiI = 0; uiI < N; uiI++)
                            if(vPos[uiI] > 0)
                                uiNumReq++;

                        size_t uiX;
                        {
                            std::unique_lock xGuard(xLock);
                            vNumPredComputed[uiIdx]++;
                            uiX = vNumPredComputed[uiIdx];
                            assert(uiX <= uiNumReq);
                            if(uiX == uiNumReq)
                                vNext.push(uiIdx);
                        }
                        if(uiX == uiNumReq)
                            xVar.notify_one();
                    }
                
                if(!bHas)
                {
                    {
                        std::unique_lock xGuard(xLock);
                        vNext.push(std::numeric_limits<size_t>::max());
                    }
                    xVar.notify_all();
                }
            }

            void process(size_t uiNumThreads, std::function<void(size_t)> fDo)
            {
                std::vector<std::thread> vThreads;
                for(size_t uiI = 0; uiI < uiNumThreads; uiI++)
                    vThreads.emplace_back(
                        [&](ParallelIterator* xIt, std::function<void(size_t)> fDo){
                        while(true)
                        {
                            size_t uiIdx = xIt->getNext();
                            if(uiIdx == std::numeric_limits<size_t>::max())
                                break;

                            fDo(uiIdx);

                            doneWith(uiIdx);
                        }
                    }, this, fDo);

                for(auto& xThread : vThreads)
                    xThread.join();
            }
    }; // class

    template <size_t N> ParallelIterator<N> genIterator(const Entry<N> xEntry) const
    {
        return ParallelIterator<N>(xEntry);
    }

    template <size_t N> ParallelIterator<N> genIterator(const Entry<N> xEntry, 
                                std::function<std::vector<size_t>(size_t, size_t, const Entry<N>&)> fSuccessors) const
    {
        return ParallelIterator<N>(xEntry, fSuccessors);
    }
};

template <typename type_defs, typename data_t> const data_t NDGrid<type_defs, data_t>::uiZero{ };


} // namespace sps
