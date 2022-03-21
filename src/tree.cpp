#include "kdpstree/default.h"

struct Settings
{
    static const size_t uiNumLayers;
};

const size_t Settings::uiNumLayers = 3;

PYBIND11_MODULE( libKdpsTree, m )
{
    pybind11::class_<Settings>( m, "SETTINGS" ).def_readonly_static( "NUM_LAYERS", &Settings::uiNumLayers );

    exportTree<OnDiskTypeDef<Settings::uiNumLayers>>( m, "KdpsTree" );

    exportArray<OnDiskTypeDef<Settings::uiNumLayers>>( m, "PsArray" );
}