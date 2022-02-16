#pragma once

#include "cstree/cont_sum.h"
#include "cstree/data_point.h"
#include "cstree/subtree_info.h"
#include "cstree/thread_pool.h"
#include "cstree/type_defs.h"
#include <string>


namespace cstree
{
template <typename type_defs> class Tree
{
  private:
    std::string sPrefix;

    template <typename data_point_t> using data_points_vec_t = typename type_defs::data_points_vec_t;
    using data_points_t = data_points_vec_t<DataPoint<type_defs>>;
    data_points_t vPoints;

    typename type_defs::bin_coords_vec_t vBinCoords;

    template <typename data_point_t> using cont_sums_vec_t = typename type_defs::cont_sums_vec_t;
    using cont_sums_t = cont_sums_vec_t<ContSum<type_defs>>;
    cont_sums_t vContSums;

    typename type_defs::point_desc_vec_t vDesc;

    SubtreeInfo<type_defs> xMainTree;

    size_t uiThreads;

    char cDelim;

    using data_points_generator = type_defs::vec_generator<data_points_t>;

    std::string pointsFileName( ) const
    {
        return sPrefix + ".points";
    }
    std::string binsFileName( ) const
    {
        return sPrefix + ".bins";
    }
    std::string sumsFileName( ) const
    {
        return sPrefix + ".sums";
    }
    std::string descFileName( ) const
    {
        return sPrefix + ".desc";
    }

  public:
    Tree( std::string sPrefix, size_t uiThreads, char cDelim )
        : sPrefix( sPrefix ),
          vPoints( data_points_generator( this->pointsFileName( ) ) ),
          vBinCoords( type_defs::vec_generator<typename type_defs::bin_coords_vec_t>( this->binsFileName( ) ) ),
          vContSums( type_defs::vec_generator<typename type_defs::cont_sums_t>( this->sumsFileName( ) ) ),
          vDesc( type_defs::vec_generator<typename type_defs::point_desc_vec_t>( this->descFileName( ) ) ),
          xMainTree( vPoints.begin( ), //
                     vPoints.end( ), //
                     vContSums.begin( ) ), //
          uiThreads( uiThreads ), //
          cDelim( cDelim ) //
    {
        // make sure the initial bounds are set up for the root node
        if( vBinCoords.size( ) == 0 && vContSums.size( ) == 0 )
        {
            for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            {
                xMainTree.vuiSubtreeBinCordsBegin[ uiD ] = 2 * uiD;
                xMainTree.vuiSubtreeBinCordsEnd[ uiD ] = 2 * uiD + 2;
                vBinCoords.push_back( 0 );
                vBinCoords.push_back( 0 );
            }
            vContSums.emplace_back( 0 );
        }
    }


    Tree( std::string sPrefix, size_t uiThreads ) : Tree( sPrefix, uiThreads, '\0' )
    {}

    Tree( std::string sPrefix ) : Tree( sPrefix, 0 )
    {}

    void addPoint( typename type_defs::point_t vPos, std::string sDesc )
    {
        typename type_defs::points_vec_offset_t uiDescOffset = vDesc.size( );
        for( char c : sDesc )
            vDesc.push_back( c );
        vDesc.push_back( cDelim );
        vPoints.emplace_back( vPos, uiDescOffset );

        for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
        {
            vBinCoords[ xMainTree.vuiSubtreeBinCordsBegin[ uiD ] ] =
                std::min( vBinCoords[ xMainTree.vuiSubtreeBinCordsBegin[ uiD ] ], vPos[ uiD ] );
            vBinCoords[ xMainTree.vuiSubtreeBinCordsEnd[ uiD ] - 1 ] =
                std::max( vBinCoords[ xMainTree.vuiSubtreeBinCordsEnd[ uiD ] - 1 ], vPos[ uiD ] );
        }

        xMainTree.uiPointsEnd++;

        vContSums.clear( );
    }

    typename type_defs::coordinate_t axisIndexFromPos( SubtreeInfo<type_defs>& rSubtreeInfo,
                                                       typename type_defs::point_t vPos,
                                                       typename type_defs::coordinate_t uiD ) const
    {
        auto iIt = std::lower_bound( vBinCoords.begin( ) + rSubtreeInfo.vuiSubtreeBinCordsBegin[ uiD ],
                                     vBinCoords.begin( ) + rSubtreeInfo.vuiSubtreeBinCordsEnd[ uiD ],
                                     vPos[ uiD ] );
        if( iIt == vBinCoords.begin( ) + rSubtreeInfo.vuiSubtreeBinCordsEnd[ uiD ] )
            return std::numeric_limits<typename type_defs::coordinate_t>::max( );
        return (typename type_defs::coordinate_t)( iIt - vBinCoords.begin( ) );
    }


    std::array<typename type_defs::coordinate_t, type_defs::d>
    axisIndicesFromPos( SubtreeInfo<type_defs>& rSubtreeInfo, typename type_defs::point_t vPos ) const
    {
        std::array<typename type_defs::coordinate_t, type_defs::d> vRet;
        for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            vRet[ uiD ] = axisIndexFromPos( rSubtreeInfo, vPos, uiD );
        return vRet;
    }

    // @todo points can be sorted into the perfect order from the getgo!!!!
    class PointSorter
    {
      private:
        const SubtreeInfo<type_defs>& rSubtreeInfo;
        const Tree* pTree;

      public:
        PointSorter( const SubtreeInfo<type_defs>& rSubtreeInfo, const Tree* pTree )
            : rSubtreeInfo( rSubtreeInfo ), pTree( pTree )
        {}

        bool operator( )( const typename type_defs::point_t& a, const typename type_defs::point_t& b ) const
        {
            typename type_defs::sums_vec_offset_t uiA =
                rSubtreeInfo.contSumIndexFromAxisIndices( pTree->axisIndicesFromPos( rSubtreeInfo, a ) );
            typename type_defs::sums_vec_offset_t uiB =
                rSubtreeInfo.contSumIndexFromAxisIndices( pTree->axisIndicesFromPos( rSubtreeInfo, b ) );
            return uiA < uiB;
        }
        typename type_defs::sums_vec_offset_t min_value( ) const
        {
            return 0;
        }
        typename type_defs::sums_vec_offset_t max_value( ) const
        {
            return pTree->vContSums.size( );
        }
    };

    template <typename type_defs::coordinate_t D>
    void integrateHelper( ThreadPool& rPool, SubtreeInfo<type_defs>& rSubtreeInfo,
                          std::array<typename type_defs::coordinate_t, type_defs::d>& vAxisIndices,
                          typename type_defs::coordinate_t uiFixedDim )
    {
        if( D == uiFixedDim )
            integrateHelper<D + 1>( rSubtreeInfo, vAxisIndices, uiFixedDim );
        else
            for( typename type_defs::points_vec_offset_t uiI = 0; uiI < rSubtreeInfo.axisSize( D ); uiX++ )
            {
                vAxisIndices[ D ] = uiI;
                integrateHelper<D + 1>( rPool, rSubtreeInfo, vAxisIndices, uiFixedDim );
            }
    }


    void subdivideBin( SubtreeInfo<type_defs>& rSubtreeInfo,
                       std::array<typename type_defs::coordinate_t, type_defs::d>
                           vBinCordsBegin,
                       std::array<typename type_defs::coordinate_t, type_defs::d>
                           vBinCordsEnd,
                       typename type_defs::coordinate_t uiLayer )
    {
        // 1: setup bin coordinates
        auto vuiAxisIndices = rSubtreeInfo.axisIndicesFromContSumIndex( uiContSumIdx );
        for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
        {
            rSubtreeInfo.vuiSubtreeBinCordsBegin[ uiD ] = vBinCoords.size( );
            auto vBinCoords = bin_generator( vBinCordsBegin[ uiD ], vBinCordsEnd[ uiD ], uiD, uiLayer );
            for( typename type_defs::coordinate_t uiBinCord : vBinCoords )
                vBinCoords.push_back( uiBinCord );
            rSubtreeInfo.vuiSubtreeBinCordsEnd[ uiD ] = vBinCoords.size( );
        }

        // 2: sort points to into their bins
        // @todo this is not necessary since points cna be sorted into the perfect order from the getgo
        type_defs::sort_func_t( vPoints.begin( ) + rSubtreeInfo.uiPointsBegin,
                                vPoints.begin( ) + rSubtreeInfo.uiPointsEnd,
                                PointSorter( rSubtreeInfo, this ) );

        // 3: setup cont sum
        rSubtreeInfo.uiSubtreeOffset = vContSums.size( );
        vContSums.resize( vContSums.size( ) + xSubcontSum.numContSumBins( ) );

        // 4: adjust point begins and ends of cont_sums & set uiVal to num
        typename type_defs::points_vec_offset_t uiCurrPointIdx = rSubtreeInfo.uiPointsBegin;
        for( typename type_defs::sums_vec_offset_t uiContSumIdx = rSubtreeInfo.uiSubtreeOffset;
             uiContSumIdx < rSubtreeInfo.countSumBinsEnd( ) + rSubtreeInfo.uiSubtreeOffset;
             uiContSumIdx++ )
        {
            vContSums[ uiContSumIdx ].xSubtree.uiPointsBegin = uiCurrPointIdx;
            while( uiCurrPointIdx < rSubtreeInfo.uiPointsEnd )
            {
                typename type_defs::sums_vec_offset_t uiX = rSubtreeInfo.contSumIndexFromAxisIndices(
                    axisIndicesFromPos( rSubtreeInfo, vPoints[ uiCurrPointIdx ] ) );
                if( uiX != uiContSumIdx )
                    break;
                uiCurrPointIdx++;
                vContSums[ uiContSumIdx ].xSubtree.uiPointsEnd = uiCurrPointIdx;
            }
            vContSums[ uiContSumIdx ].uiVal =
                vContSums[ uiContSumIdx ].xSubtree.uiPointsEnd - vContSums[ uiContSumIdx ].xSubtree.uiPointsBegin;
        }

        // 5: integrate cont_sums
        for( typename type_defs::coordinate_t uiFixedDim = 0; uiFixedDim < type_defs::d; uiFixedDim++ )
        {
            ThreadPool xPool( uiThreads );
            std::array<typename type_defs::coordinate_t, type_defs::d> vAxisIndices;
            integrateHelper<0>( xPool, rSubtreeInfo, vAxisIndices, uiFixedDim );
        }
    }

    void recursiveSubdivideHelper( SubtreeInfo<type_defs>& rSubtreeInfo, typename type_defs::coordinate_t uiLayer,
                                   typename type_defs::points_vec_offset_t uiMinPoints )
    {
        for( typename type_defs::sums_vec_offset_t uiContSumIdx = rSubtreeInfo.uiSubtreeOffset;
             uiContSumIdx < rSubtreeInfo.countSumBinsEnd( ) + rSubtreeInfo.uiSubtreeOffset;
             uiContSumIdx++ )
        {
            if( vContSums[ uiContSumIdx ].xSubtree.numPoints >= uiMinPoints )
            {
                std::array<typename type_defs::coordinate_t, type_defs::d> vBinCordsBegin =
                    rSubtreeInfo.axisIndicesFromContSumIndex( uiContSumIdx );
                std::array<typename type_defs::coordinate_t, type_defs::d> vBinCordsEnd =
                    rSubtreeInfo.axisIndicesFromContSumIndex( uiContSumIdx, 1 );
                subdivideBin( vContSums[ uiContSumIdx ].xSubtree, vBinCordsBegin, vBinCordsEnd, uiLayer );
                recursiveSubdivideHelper( vContSums[ uiContSumIdx ].xSubtree, uiLayer + 1, uiMinPoints );
            }
        }
    }

    void recursiveSubdivide( typename type_defs::points_vec_offset_t uiMinPoints )
    {
        recursiveSubdivideHelper( xMainTree, 0, uiMinPoints );
    }

    template <typename type_defs::coordinate_t D>
    void printHelper( SubtreeInfo<type_defs>& rSubtreeInfo,
                      std::array<typename type_defs::coordinate_t, type_defs::d>& vAxisIndices,
                      typename type_defs::coordinate_t uiLayer, std::string sIndet )
    {
        for( typename type_defs::points_vec_offset_t uiI = 0; uiI < rSubtreeInfo.axisSize( D ); uiX++ )
        {
            vAxisIndices[ D ] = uiI;
            printHelper<D + 1>( rSubtreeInfo, vAxisIndices, uiLayer, sIndet );
        }
    }

    std::string print( )
    {
        std::array<typename type_defs::coordinate_t, type_defs::d> vAxisIndices;
        return printHelper<0>( xMainTree, vAxisIndices, 0, "  " );
    }
};

template <typename type_defs>
void Tree::integrateHelper<type_defs::d>( ThreadPool& rPool,
                                          SubtreeInfo<type_defs>& rSubtreeInfo,
                                          std::array<typename type_defs::coordinate_t, type_defs::d>& vAxisIndices,
                                          typename type_defs::coordinate_t uiFixedDim )
{
    xPool.enqueue(
        [ & ]( size_t uiTid, std::array<typename type_defs::coordinate_t, type_defs::d> vAxisIndices ) {
            typename type_defs::cont_sum_val_t uiContSum = 0;
            for( typename type_defs::points_vec_offset_t uiI = 0; uiI < rSubtreeInfo.axisSize( uiFixedDim ); uiX++ )
            {
                vAxisIndices[ uiFixedDim ] = uiI;
                typename type_defs::sums_vec_offset_t uiContSumIdx =
                    rSubtreeInfo.contSumIndexFromAxisIndices( vAxisIndices );
                uiContSum += vContSums[ uiContSumIdx ].uiVal;
                vContSums[ uiContSumIdx ].uiVal = uiContSum;
            }
        },
        vAxisIndices );
}


template <typename type_defs>
std::string void
Tree::printHelper<type_defs::d>( SubtreeInfo<type_defs>& rSubtreeInfo,
                                 std::array<typename type_defs::coordinate_t, type_defs::d>& vAxisIndices,
                                 typename type_defs::coordinate_t uiLayer, std::string sIndet )
{
    std::string sRet = "";
    typename type_defs::sums_vec_offset_t uiContSumIdx = rSubtreeInfo.contSumIndexFromAxisIndices( vAxisIndices );
    for( typename type_defs::coordinate_t uiX = 0; uiX < uiLayer; uiX++ )
        sRet += sIndet;
    sRet += "Node: " + std::to_string( uiContSumIdx ) +
            " Cont. Sum: " + std::to_string( vContSums[ uiContSumIdx ].uiVal ) +
            " Num. Points: " + std::to_string( vContSums[ uiContSumIdx ].xSubtree.numPoints( ) ) +
            " Num. Children: " + std::to_string( vContSums[ uiContSumIdx ].xSubtree.numCountSumBins( ) );
    sRet += "\n";
    if( vContSums[ uiContSumIdx ].xSubtree.uiSubtreeOffset == 0 && uiLayer != 0 )
        for( typename type_defs::points_vec_offset_t uiP = vContSums[ uiContSumIdx ].xSubtree.uiPointsBegin;
             uiP < vContSums[ uiContSumIdx ].xSubtree.uiPointsBegin;
             uiP++ )
        {
            for( typename type_defs::coordinate_t uiX = 0; uiX < uiLayer + 1; uiX++ )
                sRet += sIndet;
            sRet += "(";
            for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            {
                if( uiD > 0 )
                    sRet += ", ";
                sRet += vPoints[ uiP ][ uiD ];
            }
            sRet += ")\n";
        }
    else
    {
        std::array<typename type_defs::coordinate_t, type_defs::d> vAxisIndices;
        sRet += printHelper<0>( vContSums[ uiContSumIdx ].xSubtree, vAxisIndices, uiLayer + 1, sIndet );
    }
    return sRet;
}

} // namespace cstree


template <typename type_defs> void exportTree( pybind11::module& m, std::string sName )
{
    pybind11::class_<cstree::Tree<type_defs>>( m, sName.c_str( ) )
        .def( pybind11::init<std::string>( ) ) // constructor
        .def( pybind11::init<std::string, size_t>( ) ) // constructor
        .def( "add_point", &cstree::Tree<type_defs>::addPoint )
        .def( "__str__", &cstree::Tree<type_defs>::print )
        .def( "subdivide", &cstree::Tree<type_defs>::recursiveSubdivide )

        ;
}
