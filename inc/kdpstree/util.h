#pragma once

#include <array>
#include <vector>


template<typename T, size_t N>
std::ostream& operator<<(std::ostream &out, const std::array<T, N>& array)
{
    out << "a[";
    if (N > 0)
        out << array[0];
    for(size_t uiI = 1; uiI < N; uiI++)
        out << ", " << array[uiI];
    out << "]";
    return out;
}

template<typename T, size_t N, size_t I>
void tupleHelp(std::ostream &out, const std::tuple<T, N>& tuple)
{
    if constexpr(I == 0)
        out << "t[";
    else
        out << ", ";
    out << std::get<I>(tuple);
    if constexpr(I+1 < N)
        tupleHelp<T, N, I+1>(out, tuple);
    else
        out << "]";
}

template<typename T, size_t N>
std::ostream& operator<<(std::ostream &out, const std::tuple<T, N>& tuple)
{
    tupleHelp<T, N, 0>(out, tuple);
    return out;
}

template<typename T>
std::ostream& operator<<(std::ostream &out, const std::vector<T>& vector)
{
    out << "v[";
    if (vector.size() > 0)
        out << vector[0];
    for(size_t uiI = 1; uiI < vector.size(); uiI++)
        out << ", " << vector[uiI];
    out << "]";
    return out;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream &out, const std::pair<T1, T2>& pair)
{
    out << "p[" << pair.first << ", " << pair.second << "]";
    return out;
}