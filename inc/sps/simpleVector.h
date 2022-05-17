#pragma once

#include "sps/point.h"
#include "sps/util.h"
#include "sps/type_defs.h"

#include <cassert>
#include <functional>
#include <string>

#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif


namespace sps
{

template <typename type_defs> class SimpleVector
{
    EXTRACT_TYPE_DEFS; // macro call

    using base_point_t = BasePoint<type_defs>;
    using ortho_t = AlignedPower2<std::array<base_point_t, 2>>;


  public:
    EXTRACT_VEC_GENERATOR( ortho, ortho_t ); // macro call
    ortho_file_t xFile;
    ortho_vec_t vData;

    
    SimpleVector( std::string sPrefix, bool bWrite )
        : xFile( ortho_vec_generator.file( sPrefix + ".orthos", bWrite ) ), vData( ortho_vec_generator.vec( xFile ) )
    {}
    

    void add( pos_t vStart, pos_t vEnd )
    {
        vData.push_back( ortho_t{ base_point_t(vStart, 0 ), base_point_t(vEnd, 0 ) } );
    }
    

    val_t count( pos_t vStart, pos_t vEnd ) const
    {
        val_t uiRet = 0;
        for(const ortho_t& rP : vData)
        {
            bool bContinue = true;
            for(size_t uiI = 0; uiI < D; uiI++)
                bContinue = bContinue && vStart[uiI] <= rP[uiI][0] && vEnd[uiI] >= rP[uiI][1];
            if(bContinue)
                uiRet += 1;
        }
        return uiRet;
    }

};


} // namespace sps


#if WITH_PYTHON
template <typename type_defs> std::string exportSimpleVector( pybind11::module& m, std::string sName)
{
    using SV = typename sps::SimpleVector<type_defs>;

    pybind11::class_<SV>( m, sName.c_str( ) )
        .def( "add", &SV::add )
        .def( "count", &OI::count )

        ;

    return "    " + sName + "\n";
}
#endif
