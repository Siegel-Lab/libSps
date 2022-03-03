#pragma once

#include "kdpstree/overlay_entries.h"
#include "kdpstree/overlay_meta.h"
#include "kdpstree/point.h"
#include "kdpstree/points.h"
#include "kdpstree/type_defs.h"


namespace kdpstree
{


template <typename type_defs> class OverlayKdTree
{
    EXTRACT_TYPE_DEFS; // macro call

    struct OverlayKdTreeNode
    {
        uint8_t uiSplitDimension;
        // split coordinate, pointer, points to leave
        std::array<std::tuple<coordinate_t, offset_t, bool>, b> vChildren;
    }; // struct

    vec_generator_t<OverlayKdTreeNode> kd_tree_vec_generator = vec_generator_t<OverlayKdTreeNode>( );
    using kd_tree_vec_t = typeof( kd_tree_vec_generator( "" ) );

    vec_generator_t<std::tuple<class_key_t, offset_t, bool>> root_vec_generator =
        vec_generator_t<std::tuple<class_key_t, offset_t, bool>>( );
    using root_vec_t = typeof( root_vec_generator( "" ) );

    using overlay_meta_t = OverlayMeta<type_defs>;
    vec_generator_t<overlay_meta_t> overlay_meta_vec_generator = vec_generator_t<overlay_meta_t>( );
    using overlay_meta_vec_t = typeof( overlay_meta_vec_generator( "" ) );

  public:
    kd_tree_vec_t vTree;
    root_vec_t vRoots;
    overlay_meta_vec_t vLeaves;

    OverlayKdTree( std::string sPrefix )
        : vTree( kd_tree_vec_generator( sPrefix + ".overlay_tree" ) ), //
          vRoots( root_vec_generator( sPrefix + ".overlay_roots" ) ), //
          vLeaves( overlay_meta_vec_generator( sPrefix + ".overlay_leaves" ) ) //
    {}

    std::string print( ) const
    {
        std::string sRet = "Roots:\n";
        for( size_t uiI = 0; uiI < vRoots.size( ); uiI++ )
            sRet += std::to_string( uiI ) + ": c" + std::to_string( std::get<0>( vRoots[ uiI ] ) ) + " -> o" +
                    std::to_string( std::get<1>( vRoots[ uiI ] ) ) +
                    ( std::get<2>( vRoots[ uiI ] ) ? " leaf\n" : " branch\n" );
        sRet += "Branches:\n";
        for( size_t uiI = 0; uiI < vTree.size( ); uiI++ )
        {
            sRet += std::to_string( uiI ) + ": d" + std::to_string( vTree[ uiI ].uiSplitDimension );
            for( size_t uiA = 0; uiA < b; uiA++ )
                if( std::get<1>( vTree[ uiI ].vChildren[ uiA ] ) != std::numeric_limits<offset_t>::max( ) )
                    sRet += " (p" + std::to_string( std::get<0>( vTree[ uiI ].vChildren[ uiA ] ) ) + ", o" +
                            std::to_string( std::get<1>( vTree[ uiI ].vChildren[ uiA ] ) ) +
                            ( std::get<2>( vTree[ uiI ].vChildren[ uiA ] ) ? ", leaf)" : ", branch)" );
            sRet += "\n";
        }
        sRet += "Leaves:\n";
        for( size_t uiI = 0; uiI < vLeaves.size( ); uiI++ )
            sRet += std::to_string( uiI ) + ": " + vLeaves[ uiI ].print( ) + "\n";
        return sRet;
    }

    void clear()
    {
        vTree.clear();
        vRoots.clear();
        vLeaves.clear();
    }
};

} // namespace kdpstree