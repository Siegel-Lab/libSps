#include "kdtree.h"
#include <pybind11/pybind11.h>

PYBIND11_MODULE( libkdtree, m )
{
    pybind11::class_<kdtree::CacheVector>( m, "kdtree" )
        .def( pybind11::init<kdtree::CacheVector::Points, std::vector<kdtree::CacheVector::Point>, size_t, size_t>( ) ) // constructor
        .def( "count", &kdtree::CacheVector::count )
        .def( "__str__", &kdtree::CacheVector::print );
}