#pragma once

#include <array>
#include <chrono>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <stxxl/vector>
#include <thread>
#include <utility>
#include <vector>

// https://stackoverflow.com/questions/24110928/overload-of-operator-not-found-when-called-from-stdostream-iterator

// has to be placed in namespace std

namespace std
{

template <typename T, size_t N> ostream& operator<<( ostream& out, const std::array<T, N>& array )
{
    out << "[";
    if( N > 0 )
        out << array[ 0 ];
    for( size_t uiI = 1; uiI < N; uiI++ )
        out << ", " << array[ uiI ];
    out << "]";
    return out;
}


template <size_t I, typename... TS> void tupleHelp( ostream& out, const std::tuple<TS...>& tuple )
{
    if constexpr( I == 0 )
        out << "(";
    else
        out << ", ";
    out << std::get<I>( tuple );
    if constexpr( I + 1 < sizeof...( TS ) )
        tupleHelp<I + 1, TS...>( out, tuple );
    else
        out << ")";
}

template <typename... TS> ostream& operator<<( ostream& out, const std::tuple<TS...>& tuple )
{
    tupleHelp<0, TS...>( out, tuple );
    return out;
}

template <typename T1, typename T2> ostream& operator<<( ostream& out, const std::pair<T1, T2>& pair )
{
    out << "(" << pair.first << ", " << pair.second << ")";
    return out;
}

template <typename T> ostream& operator<<( ostream& out, const std::vector<T>& vector )
{
    out << "{";
    if( vector.size( ) > 0 )
        out << vector[ 0 ];
    for( size_t uiI = 1; uiI < vector.size( ); uiI++ )
        out << ", " << vector[ uiI ];
    out << "}";
    return out;
}

template <typename T> ostream& stream( ostream& out, const std::vector<T>& vector, size_t uiFrom, size_t uiTo )
{
    out << "{";
    if( vector.size( ) > uiFrom )
        out << vector[ uiFrom ];
    for( size_t uiI = uiFrom + 1; uiI < uiTo; uiI++ )
        out << ", " << vector[ uiI ];
    out << "}";
    return out;
}


template <typename ValueType, unsigned PageSize, typename PagerType, unsigned BlockSize, typename AllocStr,
          typename SizeType>
ostream& operator<<( ostream& out,
                     const stxxl::vector<ValueType, PageSize, PagerType, BlockSize, AllocStr, SizeType>& vector )
{
    out << "{";
    if( vector.size( ) > 0 )
        out << vector[ 0 ];
    for( size_t uiI = 1; uiI < vector.size( ); uiI++ )
        out << ", " << vector[ uiI ];
    out << "}";
    return out;
}

template <typename ValueType, unsigned PageSize, typename PagerType, unsigned BlockSize, typename AllocStr,
          typename SizeType>
ostream& stream( ostream& out,
                 const stxxl::vector<ValueType, PageSize, PagerType, BlockSize, AllocStr, SizeType>& vector,
                 size_t uiFrom, size_t uiTo )
{
    out << "{";
    if( vector.size( ) > uiFrom )
        out << vector[ uiFrom ];
    for( size_t uiI = uiFrom + 1; uiI < uiTo; uiI++ )
        out << ", " << vector[ uiI ];
    out << "}";
    return out;
}


template <typename stream_t, typename T> optional<stream_t>& operator<<( optional<stream_t>& out, const T& rT )
{
    if( out )
        *out << rT;
    return out;
}

} // namespace std

namespace sps
{

// taken from https://www.techiedelight.com/round-next-highest-power-2/
size_t constexpr nextPower2( size_t n )
{
    // decrement `n` (to handle the case when `n` itself is a power of 2)
    n--;

    // set all bits after the last set bit
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;

    // increment `n` and return
    return ++n;
}

// align C_T to ALIGN_TO, by padding elements
template <typename C_T, size_t ALIGN_TO> class AlignTo : public C_T
{
  private:
    static_assert( sizeof( C_T ) <= ALIGN_TO );
    // make sure this is aligned to the next power of two (block load optimization of stxxl::vector)

    //#pragma GCC diagnostic push
    //#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
    // hmmm would have expected this to give and uninitialized warning... o.O
    std::array<char, ALIGN_TO - sizeof( C_T )> __buffer;
    //#pragma GCC diagnostic pop

  public:
    using C_T::C_T;

}; // class


// conditional inheritance required to force minimal memory usage (memory of disabled types is still allocated)
#define POWER_2_COND( C_T )                                                                                            \
    std::conditional<( sizeof( C_T ) < nextPower2( sizeof( C_T ) ) ), AlignTo<C_T, nextPower2( sizeof( C_T ) )>,       \
                     C_T>::type

// if the size of C_T is not an even power of 2 align it to the next even power of two
template <typename C_T> class AlignedPower2 : public POWER_2_COND( C_T )
{
    using X = typename POWER_2_COND( C_T );

  public:
    using X::X;
};


const std::string CLRLN = "\r\033[K";

#define DO_PROFILE 1

#if DO_PROFILE == 1
struct Profiler
{
    std::chrono::time_point<std::chrono::high_resolution_clock> xLastTimePoint;
    std::string sLastLabel;
    std::map<std::string, double> xTimes;
    Profiler( std::string sLabel ) : xLastTimePoint( std::chrono::high_resolution_clock::now( ) ), sLastLabel( sLabel )
    {}

    Profiler( ) : Profiler( "init" )
    {}
    ~Profiler( )
    {
        print( "" );
    }
    void step( std::string sLabel )
    {
        if( xTimes.count( sLastLabel ) == 0 )
            xTimes[ sLastLabel ] = 0;
        xTimes[ sLastLabel ] = xTimes[ sLastLabel ] + std::chrono::duration<double, std::milli>(
                                                          std::chrono::high_resolution_clock::now( ) - xLastTimePoint )
                                                          .count( );
        sLastLabel = sLabel;
        xLastTimePoint = std::chrono::high_resolution_clock::now( );
    }
    void print( std::string sLabel )
    {
        step( sLabel );
        std::cerr << std::endl;
        double dTotal = 0;
        for( auto xEntry : xTimes )
            dTotal += xEntry.second;
        for( auto xEntry : xTimes )
            std::cerr << 100 * xEntry.second / dTotal << "%\t" << xEntry.second << " ms\t" << xEntry.first << std::endl;
        xTimes.clear( );
    }
};
#else
struct Profiler
{
    Profiler( std::string )
    {}

    Profiler( )
    {}
    ~Profiler( )
    {
        print( "" );
    }
    void print( std::string )
    {}
    void step( std::string )
    {}
};
#endif

class Lockable
{
    size_t uiMaxPartialLocks;
    size_t uiNumPartialLocks = 0;
    bool bFullLock = false;
    size_t bFullLockWaiting = 0;
    std::mutex xLock;
    std::condition_variable xCv;

  public:
    Lockable( size_t uiMaxPartialLocks ) : uiMaxPartialLocks( uiMaxPartialLocks )
    {}

    class FullLock
    {
        Lockable& rX;

      public:
        FullLock( Lockable& rX ) : rX( rX )
        {
            std::unique_lock xGuard( rX.xLock );
            ++rX.bFullLockWaiting;
            while( rX.uiNumPartialLocks > 0 || rX.bFullLock )
                rX.xCv.wait( xGuard );
            assert( rX.uiNumPartialLocks == 0 && !rX.bFullLock );
            --rX.bFullLockWaiting;
            rX.bFullLock = true;
        }

        ~FullLock( )
        {
            std::unique_lock xGuard( rX.xLock );
            rX.bFullLock = false;
            rX.xCv.notify_all( );
        }
    };
    class PartialLock
    {
        Lockable& rX;

      public:
        PartialLock( Lockable& rX ) : rX( rX )
        {
            std::unique_lock xGuard( rX.xLock );
            while( rX.bFullLock || rX.uiNumPartialLocks >= rX.uiMaxPartialLocks || rX.bFullLockWaiting > 0 )
                rX.xCv.wait( xGuard );
            assert( !rX.bFullLock && rX.uiNumPartialLocks < rX.uiMaxPartialLocks && rX.bFullLockWaiting == 0 );
            ++rX.uiNumPartialLocks;
        }

        ~PartialLock( )
        {
            std::unique_lock xGuard( rX.xLock );
            assert( rX.uiNumPartialLocks > 0 );
            --rX.uiNumPartialLocks;
            if( rX.uiNumPartialLocks == 0 || rX.uiNumPartialLocks + 1 == rX.uiMaxPartialLocks )
                rX.xCv.notify_all( );
        }
    };

    FullLock fullLock( )
    {
        return FullLock( *this );
    }

    PartialLock partialLock( )
    {
        return PartialLock( *this );
    }
};

} // namespace sps