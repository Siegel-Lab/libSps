#pragma once

#include "cstree/cont_sum.h"
#include "cstree/data_point.h"
#include "cstree/subtree_info.h"
#include "cstree/thread_pool.h"
#include "cstree/type_defs.h"
#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif
#include <string>


namespace cstree
{
template <typename type_defs> class Tree
{
  private:
    std::string sPrefix;

    const typename type_defs::bin_cords_generator xBinCoordsGen;

    using data_points_t = typename type_defs::template data_points_vec_t<DataPoint<type_defs>>;
    data_points_t vPoints;


    using cont_sums_t = typename type_defs::template cont_sums_vec_t<ContSum<type_defs>>;
    cont_sums_t vContSums;

    typename type_defs::point_desc_vec_t vDesc;

    AnnotatedSubtreeInfo<type_defs> xMainTree;

    size_t uiThreads;

    char cDelim;

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

    std::array<typename type_defs::coordinate_t, type_defs::d> mainBinCoordsBegin( ) const
    {
        std::array<typename type_defs::coordinate_t, type_defs::d> vuiBinCoordsBegin;
        for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            vuiBinCoordsBegin[ uiD ] = xBinCoordsGen.minCoord( uiD );
    }

    std::array<typename type_defs::coordinate_t, type_defs::d> mainBinCoordsEnd( ) const
    {
        std::array<typename type_defs::coordinate_t, type_defs::d> vuiBinCoordsEnd;
        for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            vuiBinCoordsEnd[ uiD ] = xBinCoordsGen.maxCoord( uiD );
    }

  public:
    Tree( std::string sPrefix, const typename type_defs::bin_cords_generator xBinCoordsGen, size_t uiThreads,
          char cDelim )
        : sPrefix( sPrefix ),
          xBinCoordsGen( xBinCoordsGen ),
          vPoints( typename type_defs::template vec_generator<data_points_t>( )( this->pointsFileName( ) ) ),
          vContSums( typename type_defs::template vec_generator<cont_sums_t>( )( this->sumsFileName( ) ) ),
          vDesc( typename type_defs::template vec_generator<typename type_defs::point_desc_vec_t>( )(
              this->descFileName( ) ) ),
          xMainTree( 0, &xBinCoordsGen, &vContSums, mainBinCoordsBegin( ), mainBinCoordsEnd( ), 0 ), //
          uiThreads( uiThreads ), //
          cDelim( cDelim ) //
    {
        // make sure the initial bounds are set up for the root node
        if( vContSums.size( ) == 0 )
            vContSums.emplace_back( 0 );
    }


    Tree( std::string sPrefix, typename type_defs::bin_cords_generator xBinCoordsGen, size_t uiThreads )
        : Tree( sPrefix, xBinCoordsGen, uiThreads, '\0' )
    {}

    Tree( std::string sPrefix, typename type_defs::bin_cords_generator xBinCoordsGen )
        : Tree( sPrefix, xBinCoordsGen, 0 )
    {}

    void addPoint( typename type_defs::point_t vPos, std::string sDesc )
    {
        typename type_defs::points_vec_offset_t uiDescOffset = vDesc.size( );
        for( char c : sDesc )
            vDesc.push_back( c );
        vDesc.push_back( cDelim );
        vPoints.emplace_back( vPos, uiDescOffset );

        xMainTree.setPointsEnd( vPoints.size( ) );

        vContSums.clear( );
    }

    // @todo points can be sorted into the perfect order from the getgo!!!!
    class PointSorter
    {
      private:
        const AnnotatedSubtreeInfo<type_defs> xMainTree;

      public:
        PointSorter( const AnnotatedSubtreeInfo<type_defs> xMainTree ) : xMainTree( xMainTree )
        {}

        bool operator( )( const DataPoint<type_defs>& a, const DataPoint<type_defs>& b ) const
        {
            // filter out points with equal positions
            bool bSame = true;
            for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d && bSame; uiD++ )
                if( a.vPos[ uiD ] != b.vPos[ uiD ] )
                    bSame = false;
            if( bSame )
                return false;

            // deal with points in different positions
            // eventurally they will end up in different bins: just keep goining down with the layers
            AnnotatedSubtreeInfo<type_defs> xRoot = xMainTree;
            while( true )
            {
                // @todo recursiveley go into subtrees
                typename type_defs::sums_vec_offset_t uiA = xRoot.contSumIndexFromPos( a.vPos );
                typename type_defs::sums_vec_offset_t uiB = xRoot.contSumIndexFromPos( b.vPos );
                if( uiA != uiB )
                    return uiA < uiB;
                // uiA == uiB
                xRoot = xRoot.getSubtree( uiA );
            }
            // will never get here as only points on equal positions are always in the same bin
            // these bins are filtered out above
        }

        typename type_defs::sums_vec_offset_t min_value( ) const
        {
            return 0;
        }

        typename type_defs::sums_vec_offset_t max_value( ) const
        {
            return std::numeric_limits<typename type_defs::sums_vec_offset_t>::max( );
        }
    };

    void sortPoints( )
    {
        typename type_defs::template sort_func_t<typename data_points_t::iterator,
                                                 typeof( PointSorter( xMainTree ) )>( )(
            vPoints.begin( ), vPoints.end( ), PointSorter( xMainTree ) );
    }

    void addPoints( std::vector<std::tuple<typename type_defs::point_t, std::string>> vPoints )
    {
        for( auto xPoint : vPoints )
            addPoint( std::get<0>( xPoint ), std::get<1>( xPoint ) );
        sortPoints( );
    }


    template <typename type_defs::coordinate_t D>
    void integrateHelper( ThreadPool& rPool, AnnotatedSubtreeInfo<type_defs> xSubtreeInfo,
                          std::array<typename type_defs::coordinate_t, type_defs::d>& vAxisIndices,
                          typename type_defs::coordinate_t uiFixedDim )
    {
        // @note this constexpr is required so that the compiler does not have to evaluate an infinite loop
        if constexpr( D == type_defs::d )
        {
            rPool.enqueue(
                [ & ]( size_t uiTid, std::array<typename type_defs::coordinate_t, type_defs::d> vAxisIndices ) {
                    typename type_defs::cont_sum_val_t uiContSum = 0;
                    for( typename type_defs::points_vec_offset_t uiI = 0; uiI < xSubtreeInfo.axisSize( uiFixedDim );
                         uiI++ )
                    {
                        vAxisIndices[ uiFixedDim ] = uiI;
                        typename type_defs::sums_vec_offset_t uiContSumIdx =
                            xSubtreeInfo.contSumIndexFromAxisIndices( vAxisIndices );
                        uiContSum += vContSums[ uiContSumIdx ].uiVal;
                        vContSums[ uiContSumIdx ].uiVal = uiContSum;
                    }
                },
                vAxisIndices );
        }
        else
        {
            if( D == uiFixedDim )
                integrateHelper<D + 1>( rPool, xSubtreeInfo, vAxisIndices, uiFixedDim );
            else
                for( typename type_defs::points_vec_offset_t uiI = 0; uiI < xSubtreeInfo.axisSize( D ); uiI++ )
                {
                    vAxisIndices[ D ] = uiI;
                    integrateHelper<D + 1>( rPool, xSubtreeInfo, vAxisIndices, uiFixedDim );
                }
        }
    }

    void subdivideBin( AnnotatedSubtreeInfo<type_defs> xCurrRoot )
    {

        // 1: setup cont sum
        xCurrRoot.setSumBinsBegin( vContSums.size( ) );
        vContSums.resize( vContSums.size( ) + xCurrRoot.numContSumBins( ) );

        // 4: adjust point begins and ends of cont_sums & set uiVal to num
        typename type_defs::points_vec_offset_t uiCurrPointIdx = xCurrRoot.pointsBegin( );
        for( typename type_defs::sums_vec_offset_t uiContSumIdx = xCurrRoot.contSumBinsBegin( );
             uiContSumIdx < xCurrRoot.contSumBinsEnd( );
             uiContSumIdx++ )
        {
            auto xChild = xCurrRoot.getSubtree( uiContSumIdx );
            xChild.setPointsBegin( uiCurrPointIdx );
            while( uiCurrPointIdx < xCurrRoot.pointsEnd( ) )
            {
                typename type_defs::sums_vec_offset_t uiX =
                    xCurrRoot.contSumIndexFromPos( vPoints[ uiCurrPointIdx ].vPos );
                if( uiX != uiContSumIdx )
                    break;
                uiCurrPointIdx++;
            }
            xChild.setPointsEnd( uiCurrPointIdx );
            vContSums[ uiContSumIdx ].uiVal = xChild.numPoints( );
        }

        // 5: integrate cont_sums
        for( typename type_defs::coordinate_t uiFixedDim = 0; uiFixedDim < type_defs::d; uiFixedDim++ )
        {
            ThreadPool xPool( uiThreads );
            std::array<typename type_defs::coordinate_t, type_defs::d> vAxisIndices;
            integrateHelper<0>( xPool, xCurrRoot, vAxisIndices, uiFixedDim );
        }
    }

    void recursiveSubdivideHelper( AnnotatedSubtreeInfo<type_defs> xCurrRoot,
                                   typename type_defs::points_vec_offset_t uiMinPoints )
    {
        for( typename type_defs::sums_vec_offset_t uiContSumIdx = xCurrRoot.contSumBinsBegin( );
             uiContSumIdx < xCurrRoot.contSumBinsEnd( );
             uiContSumIdx++ )
        {
            auto xChild = xCurrRoot.getSubtree( uiContSumIdx );
            if( xChild.numPoints( ) >= uiMinPoints )
            {
                subdivideBin( xChild );
                recursiveSubdivideHelper( xChild, uiMinPoints );
            }
        }
    }

    void recursiveSubdivide( typename type_defs::points_vec_offset_t uiMinPoints )
    {
        recursiveSubdivideHelper( xMainTree, uiMinPoints );
    }

    template <typename type_defs::coordinate_t D>
    std::string printHelper( AnnotatedSubtreeInfo<type_defs> xCurrRoot,
                             std::array<typename type_defs::coordinate_t, type_defs::d>& vAxisIndices,
                             std::string sIndet )
    {
        std::string sRet = "";
        // @note this constexpr is required so that the compiler does not have to evaluate an infinite loop
        if constexpr( D == type_defs::d )
        {
            typename type_defs::sums_vec_offset_t uiContSumIdx = xCurrRoot.contSumIndexFromAxisIndices( vAxisIndices );
            auto xChild = xCurrRoot.getSubtree( uiContSumIdx );
            for( typename type_defs::coordinate_t uiX = 0; uiX < xChild.uiLayer; uiX++ )
                sRet += sIndet;
            sRet += "Node: " + std::to_string( uiContSumIdx ) +
                    " Cont. Sum: " + std::to_string( vContSums[ uiContSumIdx ].uiVal ) +
                    " Num. Points: " + std::to_string( xChild.numPoints( ) ) +
                    " Num. Children: " + std::to_string( xChild.numContSumBins( ) );
            sRet += "\n";
            if( xChild.contSumBinsBegin( ) == 0 && xChild.uiLayer != 0 )
                for( typename type_defs::points_vec_offset_t uiP = xChild.pointsBegin( ); uiP < xChild.pointsEnd( );
                     uiP++ )
                {
                    for( typename type_defs::coordinate_t uiX = 0; uiX < xChild.uiLayer + 1; uiX++ )
                        sRet += sIndet;
                    sRet += "(";
                    for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
                    {
                        if( uiD > 0 )
                            sRet += ", ";
                        sRet += vPoints[ uiP ].vPos[ uiD ];
                    }
                    sRet += ")\n";
                }
            else
            {
                std::array<typename type_defs::coordinate_t, type_defs::d> vAxisIndices;
                sRet += printHelper<0>( xChild, vAxisIndices, sIndet );
            }
        }
        else
        {
            for( typename type_defs::points_vec_offset_t uiI = 0; uiI < xCurrRoot.axisSize( D ); uiI++ )
            {
                vAxisIndices[ D ] = uiI;
                sRet += printHelper<D + 1>( xCurrRoot, vAxisIndices, sIndet );
            }
        }
        return sRet;
    }

    std::string print( )
    {
        std::array<typename type_defs::coordinate_t, type_defs::d> vAxisIndices;
        return printHelper<0>( xMainTree, vAxisIndices, "  " );
    }

    std::string printRaw( )
    {
        std::string sRet;

        sRet += "Desc (size " + std::to_string( vDesc.size( ) ) + "):\n0\t";
        size_t uiI = 0;
        for( char c : vDesc )
        {
            uiI++;
            if( c == cDelim )
                sRet += "\n" + std::to_string( uiI ) + "\t";
            else
                sRet += c;
        }
        for( size_t uiI = 0; uiI < std::to_string( uiI ).size( ) + 1; uiI++ )
            sRet.pop_back( );
        sRet += "Points (size " + std::to_string( vPoints.size( ) ) + "):\n";
        uiI = 0;
        for( auto xPoint : vPoints )
        {
            sRet += std::to_string( uiI ) + "\t(";
            for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            {
                if( uiD > 0 )
                    sRet += ", ";
                sRet += std::to_string( xPoint.vPos[ uiD ] );
            }
            sRet += ") -> ";
            sRet += std::to_string( xPoint.uiDesc );
            sRet += "\n";
            uiI++;
        }

        sRet += "Cont Sums (size " + std::to_string( vContSums.size( ) ) + "):\n";
        uiI = 0;
        for( typename type_defs::sums_vec_offset_t uiSubtreeIdx = 0; uiSubtreeIdx < vContSums.size( ); uiSubtreeIdx++ )
        {
            AnnotatedSubtreeInfo<type_defs> xSum( uiSubtreeIdx, &xBinCoordsGen, &vContSums,
                                                  std::array<typename type_defs::coordinate_t, type_defs::d>{ },
                                                  std::array<typename type_defs::coordinate_t, type_defs::d>{ }, 0 );
            sRet += std::to_string( uiI ) //
                    + "\tval: " + std::to_string( vContSums[ uiSubtreeIdx ].uiVal ) //
                    + "\n\tpoints begin: " + std::to_string( xSum.pointsBegin( ) ) //
                    + "\n\tpoints end: " + std::to_string( xSum.pointsEnd( ) ) //
                    + "\n\tsubtree offset: " + std::to_string( xSum.contSumBinsBegin( ) );
            sRet += "\n";
            uiI++;
        }

        sRet += "Root:\n";
        sRet += "points begin: " + std::to_string( xMainTree.pointsBegin( ) ) //
                + "\npoints end: " + std::to_string( xMainTree.pointsEnd( ) ) //
                + "\nsubtree offset: " + std::to_string( xMainTree.contSumBinsBegin( ) );
        sRet += "\n";


        return sRet;
    }
};


} // namespace cstree


#if WITH_PYTHON
template <typename type_defs> void exportTree( pybind11::module& m, std::string sName )
{
    pybind11::class_<cstree::Tree<type_defs>>( m, sName.c_str( ) )
        .def( pybind11::init<std::string, typename type_defs::bin_cords_generator>( ) ) // constructor
        .def( pybind11::init<std::string, typename type_defs::bin_cords_generator, size_t>( ) ) // constructor
        .def( "add_point", &cstree::Tree<type_defs>::addPoint )
        .def( "sort_points", &cstree::Tree<type_defs>::sortPoints )
        .def( "add_points", &cstree::Tree<type_defs>::addPoints )
        .def( "__str__", &cstree::Tree<type_defs>::print )
        .def( "str_raw", &cstree::Tree<type_defs>::printRaw )
        .def( "subdivide", &cstree::Tree<type_defs>::recursiveSubdivide )

        ;
}
#endif
