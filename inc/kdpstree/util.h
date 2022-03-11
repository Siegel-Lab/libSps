#pragma once

#include <array>
#include <vector>
#include <iostream>
#include <utility>


template<typename T, size_t N>
std::ostream& operator<<(std::ostream &out, const std::array<T, N>& array)
{
    out << "[";
    if (N > 0)
        out << array[0];
    for(size_t uiI = 1; uiI < N; uiI++)
        out << ", " << array[uiI];
    out << "]";
    return out;
}


template<size_t I, typename... TS>
void tupleHelp(std::ostream &out, const std::tuple<TS...>& tuple)
{
    if constexpr(I == 0)
        out << "(";
    else
        out << ", ";
    out << std::get<I>(tuple);
    if constexpr(I+1 < sizeof...(TS))
        tupleHelp<I+1, TS...>(out, tuple);
    else
        out << ")";
}

template<typename... TS>
std::ostream& operator<<(std::ostream &out, const std::tuple<TS...>& tuple)
{
    tupleHelp<0, TS...>(out, tuple);
    return out;
}

template<typename T>
std::ostream& operator<<(std::ostream &out, const std::vector<T>& vector)
{
    out << "{";
    if (vector.size() > 0)
        out << vector[0];
    for(size_t uiI = 1; uiI < vector.size(); uiI++)
        out << ", " << vector[uiI];
    out << "}";
    return out;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream &out, const std::pair<T1, T2>& pair)
{
    out << "(" << pair.first << ", " << pair.second << ")";
    return out;
}