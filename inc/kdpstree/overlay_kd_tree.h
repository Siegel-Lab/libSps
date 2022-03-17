#pragma once

#include "kdpstree/overlay_entries.h"
#include "kdpstree/overlay_meta.h"
#include "kdpstree/point.h"
#include "kdpstree/points.h"
#include "kdpstree/type_defs.h"


namespace kdpstree
{

template <size_t A, size_t B> struct AssertDiv
{
    static_assert( A % B == 0 );
    static constexpr size_t result = A / B;
};


template <typename type_defs> class OverlayKdTree
{
    EXTRACT_TYPE_DEFS; // macro call

  public:
    class UnalignedOverlayKdTreeNode
    {
      public:
        uint8_t uiSplitDimension;
        // split coordinate, pointer, points to leave
        std::array<std::tuple<coordinate_t, offset_t, bool>, b> vChildren;

        UnalignedOverlayKdTreeNode( uint8_t uiSplitDimension ) : uiSplitDimension( uiSplitDimension )
        {}

        UnalignedOverlayKdTreeNode( ) : UnalignedOverlayKdTreeNode( 0 )
        {}

    }; // struct
    class OverlayKdTreeNode : public UnalignedOverlayKdTreeNode
    {
      public:
        using UnalignedOverlayKdTreeNode::UnalignedOverlayKdTreeNode;

      private:
        static constexpr size_t ALIGN_TO = 4096;

        static_assert( sizeof( UnalignedOverlayKdTreeNode ) <= ALIGN_TO );
        // make sure sizeof(this) % 4096 == 0
        std::array<char, ALIGN_TO % sizeof( UnalignedOverlayKdTreeNode )> __buffer;
    }; // struct

    static constexpr size_t uiOverlayKdTreeNodeEle = AssertDiv<4096, sizeof( OverlayKdTreeNode )>::result;
    EXTRACT_VEC_GENERATOR_ELE( branch, OverlayKdTreeNode, uiOverlayKdTreeNodeEle ); // macro call

    using root_t = std::tuple<class_key_t, offset_t, bool>;
    EXTRACT_VEC_GENERATOR( root, root_t ); // macro call

    using overlay_meta_t = OverlayMeta<type_defs>;
    static constexpr size_t uiOverlayMetaEle = AssertDiv<4096, sizeof( overlay_meta_t )>::result;
    EXTRACT_VEC_GENERATOR_ELE( leaf, overlay_meta_t, uiOverlayMetaEle ); // macro call

    branch_file_t xTreeFile;
    branch_vec_t vTree;
    root_file_t xRootsFile;
    root_vec_t vRoots;
    leaf_file_t xLeavesFile;
    leaf_vec_t vLeaves;

    OverlayKdTree( std::string sPrefix )
        : xTreeFile( branch_vec_generator.file( sPrefix + ".overlay_branches" ) ), //
          vTree( branch_vec_generator.vec( xTreeFile ) ), //
          xRootsFile( root_vec_generator.file( sPrefix + ".overlay_roots" ) ), //
          vRoots( root_vec_generator.vec( xRootsFile ) ), //
          xLeavesFile( leaf_vec_generator.file( sPrefix + ".overlay_leaves" ) ), //
          vLeaves( leaf_vec_generator.vec( xLeavesFile ) ) //
    {}

    void clear( )
    {
        vTree.clear( );
        vRoots.clear( );
        vLeaves.clear( );
    }
};


} // namespace kdpstree

namespace std
{
template <typename type_defs>
ostream& operator<<( ostream& os, const typename kdpstree::OverlayKdTree<type_defs>& rTree )
{
    os << "Roots:" << std::endl;
    for( size_t uiI = 0; uiI < rTree.vRoots.size( ); uiI++ )
        os << uiI << ": c" << std::get<0>( rTree.vRoots[ uiI ] ) << " -> o" << std::get<1>( rTree.vRoots[ uiI ] )
           << ( std::get<2>( rTree.vRoots[ uiI ] ) ? " leaf" : " branch" ) << std::endl;
    os << "Branches:" << std::endl;
    for( size_t uiI = 0; uiI < rTree.vTree.size( ); uiI++ )
    {
        os << uiI << ": d" << (size_t)rTree.vTree[ uiI ].uiSplitDimension << " -> ";
        for( size_t uiA = 0; uiA < type_defs::b; uiA++ )
            if( std::get<1>( rTree.vTree[ uiI ].vChildren[ uiA ] ) !=
                std::numeric_limits<typename type_defs::offset_t>::max( ) )
                os << uiA << ": (p" << std::get<0>( rTree.vTree[ uiI ].vChildren[ uiA ] ) << ", o"
                   << std::get<1>( rTree.vTree[ uiI ].vChildren[ uiA ] )
                   << ( std::get<2>( rTree.vTree[ uiI ].vChildren[ uiA ] ) ? " leaf) " : " branch) " );
        os << std::endl;
    }
    os << "Leaves:" << std::endl;
    for( size_t uiI = 0; uiI < rTree.vLeaves.size( ); uiI++ )
        os << uiI << ": " << rTree.vLeaves[ uiI ] << std::endl;
    return os;
}
} // namespace std