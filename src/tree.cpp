#include "kdpstree/default.h"

struct Settings
{
    static const size_t uiNumLayers = 3;
};

const size_t Settings::uiNumLayers;

PYBIND11_MODULE( libKdpsTree, m )
{
    pybind11::class_<Settings>( m, "SETTINGS" ).def_readonly_static( "NUM_LAYERS", &Settings::uiNumLayers );

    exportTree<OnDiskTypeDef<3>>( m, "KdpsTree" );

    exportArray<OnDiskTypeDef<3>>( m, "PsArray" );
}