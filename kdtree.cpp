#include "kdtree.h"
#include <pybind11/pybind11.h>

PYBIND11_MODULE( libkdtree, m )
{
    pybind11::class_<kdtree::KDTree>( m, "kdtree" )
        .def( pybind11::init<size_t, size_t, kdtree::KDTree::Points>( ) ) // constructor
        .def( "count", &kdtree::KDTree::count );
}