#pragma once

#include <array>
#include <iostream>
#include <stxxl/vector>
#include <utility>
#include <vector>
#include <chrono>

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

template <typename T>
ostream& stream( ostream& out, const std::vector<T>& vector, size_t uiFrom, size_t uiTo )
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
    if(out)
        *out << rT;
    return out;
}

} // namespace std

namespace sps{

// taken from https://www.techiedelight.com/round-next-highest-power-2/
size_t constexpr nextPower2(size_t n)
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

const std::string CLRLN = "\r\033[K";

#define DO_PROFILE 1

#if DO_PROFILE == 1
struct Profiler
{
    std::chrono::time_point<std::chrono::high_resolution_clock> xLastTimePoint;
    std::string sLastLabel;
    std::map<std::string, double> xTimes;
    Profiler(std::string sLabel) : xLastTimePoint(std::chrono::high_resolution_clock::now()), sLastLabel(sLabel) {}
    void step(std::string sLabel)
    {
        if(xTimes.count(sLastLabel) == 0)
            xTimes[sLastLabel] = 0;
        xTimes[sLastLabel] = xTimes[sLastLabel] + std::chrono::duration<double, std::milli>(
                        std::chrono::high_resolution_clock::now() - xLastTimePoint).count();
        sLastLabel = sLabel;
        xLastTimePoint = std::chrono::high_resolution_clock::now();
    }
    void print(std::string sLabel) {
        step(sLabel);
        std::cerr << std::endl;
        for(auto xEntry : xTimes)
            std::cerr << xEntry.first << ": " << xEntry.second << "ms" << std::endl;
        xTimes.clear();
    }
    ~Profiler() { step(""); }
};
#else
struct Profiler
{
    Profiler(std::string) {}
    void step(std::string) {}
    ~Profiler() { step(""); }
    void print(std::string) {}
};
#endif

}