#pragma once

#include "cstree/cont_sum.h"
#include "cstree/type_defs.h"


namespace cstree
{

template <typename type_defs> class AnnotatedSubtreeInfo;

template <typename type_defs> class SubtreeInfo
{
    typename type_defs::points_vec_offset_t uiPointsBegin;
    typename type_defs::points_vec_offset_t uiPointsEnd;
    typename type_defs::sums_vec_offset_t uiSubtreeOffset;

  public:
    SubtreeInfo( typename type_defs::points_vec_offset_t uiPointsBegin,
                 typename type_defs::points_vec_offset_t uiPointsEnd,
                 typename type_defs::sums_vec_offset_t uiSubtreeOffset )
        : uiPointsBegin( uiPointsBegin ), uiPointsEnd( uiPointsEnd ), uiSubtreeOffset( uiSubtreeOffset )
    {}

    SubtreeInfo( ) : uiPointsBegin( 0 ), uiPointsEnd( 0 ), uiSubtreeOffset( 0 )
    {}

    friend class AnnotatedSubtreeInfo<type_defs>;
};

template <typename type_defs> class ContSum;

template <typename type_defs> class AnnotatedSubtreeInfo
{
    typename type_defs::sums_vec_offset_t uiSubtreeIdx;
    const typename type_defs::bin_cords_generator* pBinCoordsGen;

    typename type_defs::template cont_sums_vec_t<ContSum<type_defs>>* pvSubtreeArray;

    const SubtreeInfo<type_defs>& rInfo( ) const
    {
        return ( *pvSubtreeArray )[ uiSubtreeIdx ].xSubtree;
    }

    SubtreeInfo<type_defs>& rInfo( )
    {
        return ( *pvSubtreeArray )[ uiSubtreeIdx ].xSubtree;
    }

  public:
    std::array<typename type_defs::coordinate_t, type_defs::d> vuiBinCoordsBegin;
    std::array<typename type_defs::coordinate_t, type_defs::d> vuiBinCoordsEnd;
    typename type_defs::coordinate_t uiLayer;

    AnnotatedSubtreeInfo( typename type_defs::sums_vec_offset_t uiSubtreeIdx,
                          const typename type_defs::bin_cords_generator* pBinCoordsGen,
                          typename type_defs::template cont_sums_vec_t<ContSum<type_defs>>* pvSubtreeArray,
                          std::array<typename type_defs::coordinate_t, type_defs::d>
                              vuiBinCoordsBegin,
                          std::array<typename type_defs::coordinate_t, type_defs::d>
                              vuiBinCoordsEnd,
                          typename type_defs::coordinate_t uiLayer )
        : uiSubtreeIdx( uiSubtreeIdx ),
          pBinCoordsGen( pBinCoordsGen ),
          pvSubtreeArray( pvSubtreeArray ),
          vuiBinCoordsBegin( vuiBinCoordsBegin ),
          vuiBinCoordsEnd( vuiBinCoordsEnd ),
          uiLayer( uiLayer )
    {}

    typename type_defs::points_vec_offset_t pointsBegin( ) const
    {
        return rInfo().uiPointsBegin;
    }
    typename type_defs::points_vec_offset_t pointsEnd( ) const
    {
        return rInfo().uiPointsEnd;
    }

    typename type_defs::points_vec_offset_t numPoints( ) const
    {
        return pointsEnd( ) - pointsBegin( );
    }

    typename type_defs::coordinate_t axisSize( typename type_defs::coordinate_t uiD ) const
    {
        return pBinCoordsGen->axisSize( vuiBinCoordsBegin[ uiD ], vuiBinCoordsEnd[ uiD ], uiD );
    }

    typename type_defs::sums_vec_offset_t numContSumBins( ) const
    {
        typename type_defs::sums_vec_offset_t uiRet = 1;
        for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            uiRet = uiRet * axisSize( uiD );
        return uiRet;
    }

    typename type_defs::sums_vec_offset_t contSumBinsBegin( ) const
    {
        return rInfo().uiSubtreeOffset;
    }

    typename type_defs::sums_vec_offset_t contSumBinsEnd( ) const
    {
        return numContSumBins( ) + contSumBinsBegin( );
    }

    void setPointsEnd( typename type_defs::points_vec_offset_t uiX )
    {
        rInfo().uiPointsEnd = uiX;
    }
    void setPointsBegin( typename type_defs::points_vec_offset_t uiX )
    {
        rInfo().uiPointsBegin = uiX;
    }

    void setSumBinsBegin( typename type_defs::sums_vec_offset_t uiSubtreeOffset )
    {
        rInfo().uiSubtreeOffset = uiSubtreeOffset;
    }


    std::array<typename type_defs::coordinate_t, type_defs::d>
    axisIndicesFromPos( typename type_defs::point_t vPos ) const
    {
        std::array<typename type_defs::coordinate_t, type_defs::d> vRet;
        for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            vRet[ uiD ] =
                pBinCoordsGen->indexForPos( vuiBinCoordsBegin[ uiD ], vuiBinCoordsEnd[ uiD ], uiD, vPos[ uiD ] );
        return vRet;
    }

    typename type_defs::sums_vec_offset_t
    contSumIndexFromAxisIndices( std::array<typename type_defs::coordinate_t, type_defs::d> vAxisIndices ) const
    {
        typename type_defs::points_vec_offset_t uiRet = 0;
        for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            uiRet = uiRet * axisSize( uiD ) + vAxisIndices[ uiD ];
        return uiRet + rInfo().uiSubtreeOffset;
    }

    typename type_defs::sums_vec_offset_t contSumIndexFromPos( typename type_defs::point_t vPos ) const
    {
        return contSumIndexFromAxisIndices( axisIndicesFromPos( vPos ) );
    }

    std::array<typename type_defs::coordinate_t, type_defs::d>
    axisIndicesFromContSumIndex( typename type_defs::coordinate_t uiContSumIdx,
                                 typename type_defs::coordinate_t uiOffset ) const
    {
        uiContSumIdx -= rInfo().uiSubtreeOffset;
        std::array<typename type_defs::coordinate_t, type_defs::d> vRet{ };
        for( typename type_defs::coordinate_t _uiD = 0; _uiD < type_defs::d; _uiD++ )
        {
            typename type_defs::coordinate_t uiD = type_defs::d - ( _uiD + 1 );
            vRet[ uiD ] =
                pBinCoordsGen->posForIndex( vuiBinCoordsBegin[ uiD ], vuiBinCoordsEnd[ uiD ], uiD, uiContSumIdx ) +
                uiContSumIdx % axisSize( uiD ) + uiOffset;
            uiContSumIdx = uiContSumIdx / axisSize( uiD );
        }
        return vRet;
    }
    std::array<typename type_defs::coordinate_t, type_defs::d>
    axisIndicesFromContSumIndex( typename type_defs::coordinate_t uiContSumIdx ) const
    {
        return axisIndicesFromContSumIndex( uiContSumIdx, 0 );
    }

    AnnotatedSubtreeInfo getSubtree( typename type_defs::sums_vec_offset_t uiSubtreeIdx )
    {
        std::array<typename type_defs::coordinate_t, type_defs::d> vuiBinCoordsBegin;
        std::array<typename type_defs::coordinate_t, type_defs::d> vuiBinCoordsEnd;
        auto vuiAxisIndices = axisIndicesFromContSumIndex( uiSubtreeIdx );
        for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
        {
            vuiBinCoordsBegin[ uiD ] = pBinCoordsGen->posForIndex(
                this->vuiBinCoordsBegin[ uiD ], this->vuiBinCoordsEnd[ uiD ], uiD, vuiAxisIndices[ uiD ] );
            vuiBinCoordsEnd[ uiD ] = pBinCoordsGen->posForIndex(
                this->vuiBinCoordsBegin[ uiD ], this->vuiBinCoordsEnd[ uiD ], uiD, vuiAxisIndices[ uiD ] + 1 );
        }
        return AnnotatedSubtreeInfo( uiSubtreeIdx, pBinCoordsGen, pvSubtreeArray, vuiBinCoordsBegin, vuiBinCoordsEnd,
                                     uiLayer + 1 );
    }

}; // class

} // namespace cstree