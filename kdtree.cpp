#include "kdtree.h"


PYBIND11_MODULE( libkdtree, m )
{
    export_<3>(m, "kdtree_3");
    export_<2>(m, "kdtree_2");
    export_<1>(m, "kdtree_1");
}