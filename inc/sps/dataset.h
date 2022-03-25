#pragma once

#include "sps/desc.h"
#include "sps/point.h"
#include "sps/type_defs.h"
#include "sps/sparse_coordinate.h"
#include "sps/overlay.h"
#include <cassert>
#include <functional>
#include <string>


namespace sps
{

template <typename type_defs> class Dataset
{
    EXTRACT_TYPE_DEFS; // macro call

    using sparse_coord_t = SparseCoord<type_defs>;
    using overlay_t = Overlay<type_defs>;
    using overlay_grid_t = NDGrid<type_defs, overlay_t>;
    using prefix_sum_grid_t = NDGrid<type_defs, val_t>;
    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;

    std::array<typename sparse_coord_t::Entry, D> vSparseCoords;
    class_key_t xDatasetId;
    typename overlay_grid_t::template Entry<D> xOverlays;

    class DivIterator
    {
        typename points_t::EntryIterator xIt;
        size_t uiDivBy;
        size_t uiDimension;
    public:
        DivIterator(typename points_t::EntryIterator xIt, size_t uiDivBy, size_t uiDimension) :
            xIt(xIt),
            uiDivBy(uiDivBy),
            uiDimension(uiDimension)
        {}

        void operator++()
        {
            ++xIt;
        }

        const coordinate_t operator*() const
        {
            return xIt->vPos[uiDimension] / uiDivBy;
        }

        bool operator!=(const DivIterator& rOther) const
        {
            return xIt != rOther.xIt;
        }
    };

    using points_it_t = typename points_t::points_vec_t::iterator;
    struct PointsComperator
    { 
        const Dataset& rDataset;
        overlay_grid_t& rOverlays;
        sparse_coord_t& rSparseCoords;

        PointsComperator( const Dataset& rDataset, overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords ) : 
            rDataset( rDataset ), rOverlays(rOverlays), rSparseCoords(rSparseCoords)
        {}

        bool operator( )( const point_t& a, const point_t& b ) const
        {
            return rDataset.overlayIndex(rOverlays, rSparseCoords, a.vPos) < 
                    rDataset.overlayIndex(rOverlays, rSparseCoords, b.vPos);
        }

        point_t min_value( ) const
        {
            return point_t(rOverlays.posOf(rDataset.xOverlays.uiStartIndex, rDataset.xOverlays), 0);
        };

        point_t max_value( ) const
        {
            return point_t(rOverlays.posOf(rDataset.xOverlays.uiStartIndex + rOverlays.sizeOf(rDataset.xOverlays),
                                           rDataset.xOverlays), 
                           0);
        };
    };
    sort_func_t<points_it_t, PointsComperator> sort_points = sort_func_t<points_it_t, PointsComperator>( );

  public:
    Dataset() 
    {}

    Dataset(overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords, prefix_sum_grid_t& rPrefixSums, 
            points_t& vPoints, class_key_t xDatasetId, typename points_t::Entry xPoints) :
            xDatasetId(xDatasetId)
    {
        // generate the overall sparse coordinates
        size_t uiNumPerDimension = (size_t)std::pow(xPoints.uiEndIndex - xPoints.uiStartIndex, 1.0/(D*D));
        pos_t vMaxCoord{};
        vPoints.iterate([&](const point_t& xPoint){
            for(size_t uiI = 0; uiI < D; uiI++)
                vMaxCoord[uiI] = std::max(vMaxCoord[uiI], xPoint.vPos[uiI]);
        }, xPoints);

        for(size_t uiI = 0; uiI < D; uiI++)
        {
            size_t uiBlockSize = vMaxCoord[uiI] / uiNumPerDimension;
            vPoints.sortByDim(uiI, xPoints);
            vSparseCoords[uiI] = rSparseCoords.add(
                DivIterator(vPoints.cbegin(xPoints), uiBlockSize, uiI), 
                DivIterator(vPoints.cend(xPoints), uiBlockSize, uiI));
        }
        // generate overlay grid
        xOverlays = rOverlays.add(rSparseCoords.axisSizes(vSparseCoords));
        
        // sort points so that they match the overlay grid order
        sort_points(vPoints.vData.begin( ) + xPoints.uiStartIndex, vPoints.vData.begin( ) + xPoints.uiEndIndex, 
                    PointsComperator(*this, rOverlays, rSparseCoords ));

        // generate all overlays
        typename points_t::Entry xCurrPoints{};
        for(coordinate_t uiI = 0; uiI < rOverlays.sizeOf(xOverlays); uiI++)
        {
            // collect points for overlay uiI
            while(xCurrPoints.uiEndIndex < vPoints.size() && 
                    overlayIndex(rOverlays, rSparseCoords, vPoints.vData[xCurrPoints.uiEndIndex].vPos) == uiI)
                ++xCurrPoints.uiEndIndex;

            // get bottom left position (compressed)
            pos_t vPos = rOverlays.posOf(uiI + xOverlays.uiStartIndex, xOverlays);

            // collect direct predecessor overlays for each dimension
            std::array<overlay_t*, D> vPredecessors;
            bool bAllLargerZero = true;
            for(size_t uiI = 0; uiI < D; uiI++)
                if(vPos[uiI] > 0)
                {
                    --vPos[uiI];
                    vPredecessors[uiI] = &rOverlays.get(vPos, xOverlays);
                    ++vPos[uiI];
                }
                else
                {
                    vPredecessors[uiI] = nullptr;
                    bAllLargerZero = false;
                }
            
            // compute the prefix sum of the overlay to the bottom-left-front-...
            val_t uiBottomLeft = 0;
            if (bAllLargerZero)
            {
                for(size_t uiI = 0; uiI < D; uiI++)
                    --vPos[uiI];
                auto& rPreEntries = rOverlays.get(vPos, xOverlays).xInternalEntires;
                uiBottomLeft = rPrefixSums.vData[
                    rPreEntries.uiStartIndex + rPrefixSums.sizeOf(rPreEntries) - 1
                ];
                for(size_t uiI = 0; uiI < D; uiI++)
                    ++vPos[uiI];
            }

            // generate the overlay
            rOverlays.vData[xOverlays.uiStartIndex].generate(
                rSparseCoords,
                rPrefixSums,
                vPoints,
                xCurrPoints,
                vPredecessors,
                vPos,
                uiBottomLeft,
                vSparseCoords
            );

            // prepare for collecting the next set of points
            xCurrPoints.uiStartIndex = xCurrPoints.uiEndIndex;
        }
    }

    coordinate_t overlayIndex(overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords, const pos_t& vPos) const
    {
        return rOverlays.indexOf(rSparseCoords.sparse(vPos, vSparseCoords), xOverlays);
    }
    
    val_t get(const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords, 
              const prefix_sum_grid_t& rPrefixSums, const pos_t& vPos) const
    {
        auto vSparsePos = rSparseCoords.sparse(vPos, vSparseCoords);
        return rOverlays.get(vSparsePos, xOverlays).get( rSparseCoords, rPrefixSums, vSparsePos );
    }
};


} // namespace sps

namespace std
{


} // namespace std