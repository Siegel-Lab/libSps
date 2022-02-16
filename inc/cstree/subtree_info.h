#pragma once

#include "cstree/cont_sum.h"
#include "cstree/type_defs.h"


namespace cstree
{
template <typename type_defs> class SubtreeInfo
{
  private:
    type_defs::points_vec_offset_t uiPointsBegin;
    type_defs::points_vec_offset_t uiPointsEnd;
    type_defs::sums_vec_offset_t uiSubtreeOffset;
    // @todo bin cordinated need not be stored they can be dynamically computed as needed
    std::array<type_defs::points_vec_offset_t, type_defs::d> vuiSubtreeBinCordsBegin;
    std::array<type_defs::points_vec_offset_t, type_defs::d> vuiSubtreeBinCordsEnd;

  public:
    SubtreeInfo( typename type_defs::points_vec_offset_t uiPointsBegin,
                 typename type_defs::points_vec_offset_t uiPointsEnd,
                 typename type_defs::sums_vec_offset_t itSubtreeOffset )
        : uiPointsBegin( uiPointsBegin ),
          uiPointsEnd( uiPointsEnd ),
          itSubtreeOffset( itSubtreeOffset ),
          vuiSubtreeBinCordsBegin{ },
          vuiSubtreeBinCordsEnd{ }
    {}

    SubtreeInfo( )
        : uiPointsBegin( 0 ),
          uiPointsEnd( 0 ),
          itSubtreeOffset( 0 ),
          vuiSubtreeBinCordsBegin{ },
          vuiSubtreeBinCordsEnd{ }
    {}

    type_defs::points_vec_offset_t numPoints( ) const
    {
        return uiPointsEnd - uiPointsBegin;
    }

    type_defs::coordinate_t axisSize( type_defs::coordinate_t uiD ) const
    {
        return vuiSubtreeBinCordsEnd[ uiD ] - vuiSubtreeBinCordsBegin[ uiD ];
    }

    type_defs::sums_vec_offset_t numCountSumBins( ) const
    {
        type_defs::sums_vec_offset_t uiRet = 1;
        for( type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            uiRet = uiRet * axisSize( uiD );
        return uiRet;
    }
    type_defs::sums_vec_offset_t countSumBinsEnd( ) const
    {
        return numCountSumBins( ) + uiSubtreeOffset;
    }

    type_defs::sums_vec_offset_t
    countSumIndexFromAxisIndices( std::array<type_defs::coordinate_t, type_defs::d> vAxisIndices ) const
    {
        type_defs::points_vec_offset_t uiRet = 0;
        for( type_defs::coordinate_t uiD = 0; uiD < type_defs::d; uiD++ )
            uiRet = uiRet * axisSize( uiD ) + vAxisIndices[ uiD ];
        return uiRet + uiSubtreeOffset;
    }

    std::array<type_defs::coordinate_t, type_defs::d>
    axisIndicesFromContSumIndex( type_defs::coordinate_t uiContSumIdx, type_defs::coordinate_t uiOffset ) const
    {
        uiContSumIdx -= uiSubtreeOffset;
        std::array<type_defs::coordinate_t, type_defs::d> vRet{ };
        for( type_defs::coordinate_t _uiD = 0; _uiD < type_defs::d; _uiD++ )
        {
            type_defs::coordinate_t uiD = type_defs::d - ( _uiD + 1 );
            vRet[ uiD ] = vuiSubtreeBinCordsBegin[ uiD ] + uiContSumIdx % axisSize( uiD ) + uiOffset;
            uiContSumIdx = uiContSumIdx / axisSize( uiD );
        }
        assert( uiContSumIdx == 0 );
        return vRet;
    }
    std::array<type_defs::coordinate_t, type_defs::d>
    axisIndicesFromContSumIndex( type_defs::coordinate_t uiContSumIdx ) const
    {
        return axisIndicesFromContSumIndex( uiContSumIdx, 0 );
    }
}
} // namespace cstree