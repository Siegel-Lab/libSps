#pragma once

#include "sps/type_defs.h"
#include "sps/util.h"
#include <cassert>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <string>
#include <thread>


namespace sps
{

template <typename type_defs, typename data_t, template <typename> typename data_tmpl_vec_generator_t> class NDGrid
{
    EXTRACT_TYPE_DEFS; // macro call

    EXTRACT_VEC_GENERATOR( data, data_t ); // macro call

  public:
    static constexpr bool THREADSAVE = data_THREADSAVE;
    std::mutex xResizeLock;
    std::mutex xRWLock;
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
        : xFile( data_vec_generator.file( sFileName, bWrite ) ), vData( data_vec_generator.vec( xFile ) )
    {}

    template <size_t N, bool SANITY = true>
    static coordinate_t indexOf( const std::array<coordinate_t, N>& vX, const Entry<N>& rInfo )
    {
        if constexpr( SANITY )
            if( rInfo.uiStartIndex == std::numeric_limits<coordinate_t>::max( ) )
                return std::numeric_limits<coordinate_t>::max( );
        assert(rInfo.uiStartIndex != std::numeric_limits<coordinate_t>::max( ));
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

    template <size_t N, bool SANITY_UNINITIALIZED = true, bool SANITY_OUT_OF_BOUND = true>
    const data_t& get( const std::array<coordinate_t, N>& vX, const Entry<N>& rInfo ) const
    {
        auto uiIdx = indexOf<N, SANITY_UNINITIALIZED>( vX, rInfo );
        if constexpr(SANITY_OUT_OF_BOUND)
            if( uiIdx == std::numeric_limits<coordinate_t>::max( ) )
                return uiZero;
        assert(uiIdx != std::numeric_limits<coordinate_t>::max( ));
        return vData[ uiIdx ];
    }

    template <size_t N>
    inline __attribute__( ( always_inline ) ) data_t& get( const std::array<coordinate_t, N>& vX,
                                                           const Entry<N>& rInfo )
    {
        auto uiIdx = indexOf<N>( vX, rInfo );
        assert( uiIdx != std::numeric_limits<coordinate_t>::max( ) );
        return vData[ uiIdx ];
    }

    template <size_t N, bool ZERO_INIT = true, bool CAPACITY_INC_ALLOWED = false>
    Entry<N> add( const std::array<coordinate_t, N>& vAxisSizes )
    {
        Entry<N> xRet{ };
        for( size_t uiI = 0; uiI < N; uiI++ )
            xRet.vAxisSizes[ uiI ] = vAxisSizes[ uiI ];

        xRet.uiStartIndex = 0; // required to get result out of sizeOf function
        size_t uiToAdd = sizeOf( xRet );
        {
            // resize the vector in a locked fashion (this just increases the size variable, no allocation happens)
            // hence it is threadsave to query and or write the vector at the same time
            std::lock_guard<std::mutex> xGuard( xResizeLock );
            // make sure no reallocation occurs on the vector
            assert( CAPACITY_INC_ALLOWED || vData.capacity( ) >= uiToAdd + vData.size( ) );
            xRet.uiStartIndex = vData.size( );
            // make space for the new grid
            vData.resize( uiToAdd + vData.size( ) );
        } // scope for xGuard

        // stxxl does not zero initialize the resized elements :'(
        if constexpr( ZERO_INIT )
            for( size_t uiI = 0; uiI < uiToAdd; uiI++ )
                vData[ xRet.uiStartIndex + uiI ] = data_t{ };

        return xRet;
    }

    template <size_t N, bool CAPACITY_INC_ALLOWED = false>
    Entry<N> add( const std::array<coordinate_t, N>& vAxisSizes, const std::vector<data_t>& vTmp )
    {
        Entry<N> xRet{ };
        for( size_t uiI = 0; uiI < N; uiI++ )
            xRet.vAxisSizes[ uiI ] = vAxisSizes[ uiI ];

        xRet.uiStartIndex = 0; // required to get result out of sizeOf function
        size_t uiToAdd = sizeOf( xRet );
        assert( uiToAdd == vTmp.size( ) );

        if constexpr( THREADSAVE )
        {
            {
                // resize the vector in a locked fashion (this just increases the size variable, no allocation happens)
                // hence it is threadsave to query and or write the vector at the same time
                std::lock_guard<std::mutex> xGuard( xResizeLock );
                // make sure no reallocation occurs on the vector
                assert( CAPACITY_INC_ALLOWED || vData.capacity( ) >= uiToAdd + vData.size( ) );
                xRet.uiStartIndex = vData.size( );
                // make space for the new grid
                vData.resize( uiToAdd + vData.size( ) );
            } // scope for xGuard
            for( size_t uiI = 0; uiI < uiToAdd; uiI++ )
                vData[ xRet.uiStartIndex + uiI ] = vTmp[ uiI ];
        }
        else
        {
            std::lock_guard<std::mutex> xGuard( xRWLock );
            xRet.uiStartIndex = vData.extend( vTmp );
        }

        return xRet;
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

    template <size_t N, typename fSuccessor_t, typename... args_successor_t> class ParallelIterator
    {
        const Entry<N> xEntry;
        std::vector<int8_t> vNumPredComputed;
        size_t uiNumDone = 0;
        std::queue<size_t> vNext;
        std::mutex xLock;
        std::condition_variable xVar;
        fSuccessor_t fSuccessors;
        // std::function<std::vector<size_t>( size_t, size_t, const Entry<N>& )> fSuccessors;
        std::tuple<args_successor_t&...> xStoredArgs;

      public:
        ParallelIterator( const Entry<N> xEntry,
                          // std::function<std::vector<size_t>( size_t, size_t, const Entry<N>& )> fSuccessors =
                          fSuccessor_t fSuccessors,
                          args_successor_t&... vArgs )
            : xEntry( xEntry ),
              vNumPredComputed( sizeOf( xEntry ) ),
              vNext( ),
              xLock( ),
              xVar( ),
              fSuccessors( fSuccessors ),
              xStoredArgs( vArgs... )
        {
            vNext.push( xEntry.uiStartIndex );
        }

        size_t getNext( )
        {
            std::unique_lock xGuard( xLock );
            while( vNext.size( ) == 0 )
                xVar.wait( xGuard );
            assert( vNext.size( ) > 0 );
            size_t uiRet = vNext.front( );
            if( uiRet != std::numeric_limits<size_t>::max( ) )
                vNext.pop( );
            return uiRet;
        }


        void doneWith( size_t uiIdx )
        {
            for( size_t uiD = 0; uiD < N; uiD++ )
            {
                for( size_t uiIdx : std::apply(
                         [ & ]( args_successor_t&... vArgs ) { return fSuccessors( uiIdx, uiD, xEntry, vArgs... ); },
                         xStoredArgs ) )
                {
                    if( uiIdx == std::numeric_limits<size_t>::max( ) ) // poison
                    {
                        {
                            std::unique_lock xGuard( xLock );
                            vNext.push( std::numeric_limits<size_t>::max( ) );
                        }
                        xVar.notify_all( );
                        return;
                    }

                    size_t uiNumReq = 0;
                    auto vPos = NDGrid::posOf( uiIdx, xEntry );
                    for( size_t uiI = 0; uiI < N; uiI++ )
                        if( vPos[ uiI ] > 0 )
                            uiNumReq++;

                    size_t uiX;
                    {
                        std::unique_lock xGuard( xLock );
                        vNumPredComputed[ uiIdx - xEntry.uiStartIndex ]++;
                        uiX = vNumPredComputed[ uiIdx - xEntry.uiStartIndex ];
                        assert( uiX <= uiNumReq );
                        if( uiX == uiNumReq )
                            vNext.push( uiIdx );
                    }
                    if( uiX == uiNumReq )
                        xVar.notify_one( );
                }
            }
            {
                std::unique_lock xGuard( xLock );
                ++uiNumDone;
            }
        }

        void process( size_t uiNumThreads, progress_stream_t& xProg, std::function<void( size_t, size_t )> fDo )
        {
            auto fTask = [ & ]( ParallelIterator* xIt, size_t uiTid, std::function<void( size_t, size_t )> fDo ) {
                while( true )
                {
                    size_t uiIdx = xIt->getNext( );
                    if( uiIdx == std::numeric_limits<size_t>::max( ) )
                        break;

                    fDo( uiTid, uiIdx );

                    doneWith( uiIdx );

                    if( uiTid == 0 && xProg.printAgain( ) )
                        xProg << Verbosity( 0 ) << "processed " << uiNumDone << " out of " << vNumPredComputed.size( )
                              << " overlays, thats " << 100.0 * ( (double)uiNumDone / (double)vNumPredComputed.size( ) )
                              << "%.\n";
                }
            };
            if( uiNumThreads == 0 || vNumPredComputed.size( ) == 0 )
                fTask( this, 0, fDo );
            else
            {
                std::vector<std::thread> vThreads;
                for( size_t uiI = 0; uiI < uiNumThreads && uiI < vNumPredComputed.size( ); uiI++ )
                    vThreads.emplace_back( std::bind( fTask, this, uiI, fDo ) );

                for( auto& xThread : vThreads )
                    xThread.join( );
            }
        }
    }; // class

    template <size_t N> using default_func_t = std::function<std::vector<size_t>( size_t, size_t, const Entry<N>& )>;

    template <size_t N> ParallelIterator<N, default_func_t<N>> genIterator( const Entry<N> xEntry ) const
    {
        return ParallelIterator<N, default_func_t<N>>(
            xEntry, []( size_t uiIdx, size_t uiDim, const Entry<N>& xEntry ) {
                auto vPos = NDGrid::posOf( uiIdx, xEntry );
                ++vPos[ uiDim ];
                std::vector<size_t> vRet;
                if( vPos[ uiDim ] < xEntry.vAxisSizes[ uiDim ] )
                    vRet.push_back( NDGrid::indexOf( vPos, xEntry ) );
                if( uiIdx + 1 == NDGrid::sizeOf( xEntry ) + xEntry.uiStartIndex )
                    vRet.push_back( std::numeric_limits<size_t>::max( ) ); // poison
                return vRet;
            } );
    }

    template <size_t N, typename fSuccessor_t, typename... args_successor_t>
    ParallelIterator<N, fSuccessor_t, args_successor_t...>
    genIterator( const Entry<N> xEntry, fSuccessor_t fSuccessors, args_successor_t&... vArgs ) const
    {
        return ParallelIterator<N, fSuccessor_t, args_successor_t...>( xEntry, fSuccessors, vArgs... );
    }
};

template <typename type_defs, typename data_t, template <typename> typename data_tmpl_vec_generator_t>
const data_t NDGrid<type_defs, data_t, data_tmpl_vec_generator_t>::uiZero{ };


} // namespace sps
