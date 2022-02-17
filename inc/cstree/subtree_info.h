#pragma once

#include "cstree/cont_sum.h"
#include "cstree/type_defs.h"


namespace cstree
{
template <typename type_defs> class SubtreeInfo
{
  public:
    typename type_defs::points_vec_offset_t uiPointsBegin;
    typename type_defs::points_vec_offset_t uiPointsEnd;
    typename type_defs::sums_vec_offset_t uiSubtreeOffset;
    // @todo bin cordinated need not be stored they can be dynamically computed as needed
    std::array<typename type_defs::points_vec_offset_t, type_defs::d> vuiSubtreeBinCordsBegin;
    std::array<typename type_defs::points_vec_offset_t, type_defs::d> vuiSubtreeBinCordsEnd;

    SubtreeInfo( typename type_defs::points_vec_offset_t uiPointsBegin,
                 typename type_defs::points_vec_offset_t uiPointsEnd,
                 typename type_defs::sums_vec_offset_t uiSubtreeOffset )
        : uiPointsBegin( uiPointsBegin ),
          uiPointsEnd( uiPointsEnd ),
          uiSubtreeOffset( uiSubtreeOffset ),
          vuiSubtreeBinCordsBegin{ },
          vuiSubtreeBinCordsEnd{ }
    {}

    SubtreeInfo( )
        : uiPointsBegin( 0 ),
          uiPointsEnd( 0 ),
          uiSubtreeOffset( 0 ),
          vuiSubtreeBinCordsBegin{ },
          vuiSubtreeBinCordsEnd{ }
    {}

    typename type_defs::points_vec_offset_t numPoints( ) const
    {
        return uiPointsEnd - uiPointsBegin;
    }

    typename type_defs::coordinate_t axisSize( typename type_defs::coordinate_t uiD ) const
    {
        return vuiSubtreeBinCordsEnd[ uiD ] - vuiSubtreeBinCordsBegin[ uiD ];
    }

    typename type_defs::sums_vec_offset_t numContSumBins( ) const
    {
        typename type_defs::sums_vec_offset_t uiRet = 1;
        for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            uiRet = uiRet * axisSize( uiD );
        return uiRet;
    }
    typename type_defs::sums_vec_offset_t contSumBinsEnd( ) const
    {
        return numContSumBins( ) + uiSubtreeOffset;
    }

    typename type_defs::sums_vec_offset_t
    contSumIndexFromAxisIndices( std::array<typename type_defs::coordinate_t, type_defs::d> vAxisIndices ) const
    {
        typename type_defs::points_vec_offset_t uiRet = 0;
        for( typename type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            uiRet = uiRet * axisSize( uiD ) + vAxisIndices[ uiD ];
        return uiRet + uiSubtreeOffset;
    }

    std::array<typename type_defs::coordinate_t, type_defs::d>
    axisIndicesFromContSumIndex( typename type_defs::coordinate_t uiContSumIdx,
                                 typename type_defs::coordinate_t uiOffset ) const
    {
        uiContSumIdx -= uiSubtreeOffset;
        std::array<typename type_defs::coordinate_t, type_defs::d> vRet{ };
        for( typename type_defs::coordinate_t _uiD = 0; _uiD < type_defs::d; _uiD++ )
        {
            typename type_defs::coordinate_t uiD = type_defs::d - ( _uiD + 1 );
            vRet[ uiD ] = vuiSubtreeBinCordsBegin[ uiD ] + uiContSumIdx % axisSize( uiD ) + uiOffset;
            uiContSumIdx = uiContSumIdx / axisSize( uiD );
        }
        return vRet;
    }
    std::array<typename type_defs::coordinate_t, type_defs::d>
    axisIndicesFromContSumIndex( typename type_defs::coordinate_t uiContSumIdx ) const
    {
        return axisIndicesFromContSumIndex( uiContSumIdx, 0 );
    }
};
} // namespace cstree