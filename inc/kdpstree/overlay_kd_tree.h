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

    void clear()
    {
        vTree.clear();
        vRoots.clear();
        vLeaves.clear();
    }
};

template <typename type_defs>
std::ostream& operator<<(std::ostream& os, const OverlayKdTree<type_defs>& rTree)
{
    os << "Roots:" << std::endl;
    for( size_t uiI = 0; uiI < rTree.vRoots.size( ); uiI++ )
        os << uiI << ": c" << std::get<0>( rTree.vRoots[uiI] ) << " -> o" 
           << std::get<1>( rTree.vRoots[uiI] ) << ( std::get<2>( rTree.vRoots[uiI] ) ? " leaf" : " branch" ) 
           << std::endl;
    os << "Branches:" << std::endl;
    for( size_t uiI = 0; uiI < rTree.vTree.size( ); uiI++ )
    {
        os << uiI << ": d" << rTree.vTree[ uiI ].uiSplitDimension;
        for( size_t uiA = 0; uiA < type_defs::b; uiA++ )
            if( std::get<1>( rTree.vTree[ uiI ].vChildren[ uiA ] ) != std::numeric_limits<typename type_defs::offset_t>::max( ) )
                os << uiA << ": (p" << std::get<0>(rTree.vTree[ uiI ].vChildren[ uiA ]) << ", o" 
                   << std::get<1>( rTree.vTree[ uiI ].vChildren[ uiA ] ) 
                   << ( std::get<2>( rTree.vTree[uiI].vChildren[ uiA ] ) ? " leaf) " : " branch) " );
        os << std::endl;
    }
    os << "Leaves:" << std::endl;
    for( size_t uiI = 0; uiI < rTree.vLeaves.size( ); uiI++ )
        os << uiI << ": " << rTree.vLeaves[ uiI ] << std::endl;
    return os;
}


} // namespace kdpstree