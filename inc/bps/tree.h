#pragma once

#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif
#include <string>


namespace bps
{

template <typename type_defs> class Tree
{

    Desc<type_defs> vDesc;

  public:
    Tree( std::string sPrefix ) : vDesc( sPrefix )
    {}

}; // class

} // namespace bps