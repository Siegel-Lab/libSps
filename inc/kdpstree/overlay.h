#pragma once

#include "kdpstree/points.h"
#include "kdpstree/nd_grid.h"
#include "kdpstree/sparse_coordinate.h"
#include <cassert>
#include <functional>
#include <string>


namespace kdpstree
{

template <typename type_defs> class Overlay
{
    EXTRACT_TYPE_DEFS; // macro call

    using sparse_coord_t = SparseCoord<type_defs>;
    using prefix_sum_grid_t = NDGrid<type_defs, val_t>;


    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;

    using cord_it_t = typename sparse_coord_t::EntryIterator;
    using point_it_t = typename points_t::EntryIterator;

    using red_pos_t = std::array<coordinate_t, D-1>;

    using entry_arr_t = std::array<typename sparse_coord_t::Entry, D>;
    using red_entry_arr_t = std::array<typename sparse_coord_t::Entry, D-1>;

    class MergeIterator
    {
        cord_it_t xA, xB;
        cord_it_t xAE, xBE;
    public:
        MergeIterator(cord_it_t xA, cord_it_t xB, cord_it_t xAE, cord_it_t xBE) :
            xA(xA),
            xB(xB),
            xAE(xAE),
            xBE(xBE)
        {}

        void operator++()
        {
            if(xA != xAE && xB != xBE && xA->first == xB->first)
            {
                ++xA;
                ++xB;
            }
            else if(xA != xAE && (!(xB != xBE) || xA->first < xB->first))
                ++xA;
            else if(xB != xBE && (!(xA != xAE) || xB->first < xA->first))
                ++xB;
        }

        const coordinate_t& operator*() const
        {
            if(xA != xAE && xA->first < xB->first)
                return xA->first;
            return xB->first;
        }

        bool operator!=(const MergeIterator& rOther) const
        {
            return xA != rOther.xA || xB != rOther.xB;
        }
    };
    
    
    class PointIterator
    {
        const sparse_coord_t& rSparseCoords;
        const entry_arr_t& rGlobalSparseCords;
        const Overlay& rOverlay;
        point_it_t xIt;
        size_t uiDim;
    public:
        PointIterator(const sparse_coord_t& rSparseCoords, 
                      const entry_arr_t& rGlobalSparseCords, const Overlay& rOverlay,
                      point_it_t xIt, size_t uiDim) :
            rSparseCoords(rSparseCoords),
            rGlobalSparseCords(rGlobalSparseCords),
            rOverlay(rOverlay),
            xIt(xIt),
            uiDim(uiDim)
        {}

        void operator++()
        {
            ++xIt;
        }

        const coordinate_t operator*() const
        {
            return rSparseCoords.replace(rSparseCoords.replace(xIt->vPos[uiDim], rGlobalSparseCords[uiDim]),
                                         rOverlay.vSparseCoordsInternal[uiDim]);
        }

        bool operator!=(const PointIterator& rOther) const
        {
            return xIt != rOther.xIt;
        }
    };


  public:
    entry_arr_t vSparseCoordsOverlay;
    entry_arr_t vSparseCoordsInternal;
    std::array<typename prefix_sum_grid_t::template Entry<D-1>, D> vOverlayEntries;
    typename prefix_sum_grid_t::template Entry<D> xInternalEntires;
    typename points_t::Entry xPoints;
    val_t uiBottomLeft;

    Overlay() :
        vSparseCoordsOverlay{}, vSparseCoordsInternal{}, vOverlayEntries{}, xInternalEntires{}, uiBottomLeft{}
    {}

    void generate( sparse_coord_t& rSparseCoords, prefix_sum_grid_t& rPrefixSums, points_t& vPoints, 
            typename points_t::Entry xPoints, std::array<Overlay*, D> vPredecessors, std::array<coordinate_t, D> vPos, 
             val_t uiBottomLeft, const entry_arr_t& rGlobalSparseCords )
    {
        this->xPoints = xPoints;
        this->uiBottomLeft = uiBottomLeft;
        // construct sparse coordinates for each dimension
        for(size_t uiI = 0; uiI < D; uiI++)
            if(vPredecessors[uiI] != nullptr)
            {
                // add coordinates from previous overlay to the overlay entries
                auto xA = rSparseCoords.cbegin(vPredecessors[uiI]->vSparseCoordsOverlay[uiI]);
                auto xAE = rSparseCoords.cend(vPredecessors[uiI]->vSparseCoordsOverlay[uiI]);

                auto xB = rSparseCoords.cbegin(vPredecessors[uiI]->vSparseCoordsInternal[uiI]);
                auto xBE = rSparseCoords.cend(vPredecessors[uiI]->vSparseCoordsInternal[uiI]);

                MergeIterator xBegin( xA, xB, xAE, xBE );
                MergeIterator xEnd( xAE, xBE, xAE, xBE );

                vSparseCoordsOverlay[uiI] = rSparseCoords.add(xBegin, xEnd);
            }

        for(size_t uiI = 0; uiI < D; uiI++)
        {
            // add coordinates from previous overlay to the overlay entries 
            vPoints.sortByDim(uiI, xPoints);
            vSparseCoordsInternal[uiI] = rSparseCoords.add(
                    PointIterator(rSparseCoords, rGlobalSparseCords, *this, vPoints.cbegin(xPoints), uiI),
                    PointIterator(rSparseCoords, rGlobalSparseCords, *this, vPoints.cend(xPoints), uiI)
                );
        }

        // construct internal grid
        pos_t vInternalAxisSizes = rSparseCoords.axisSizes(vSparseCoordsInternal);
        xInternalEntires = rPrefixSums.add(vInternalAxisSizes);
        typename prefix_sum_grid_t::template Entry<D> xInternalEntiresHelper = rPrefixSums.add(vInternalAxisSizes);
        vPoints.iterate([&](const point_t& xPoint){
            rPrefixSums.get(
                    rSparseCoords.sparse(rSparseCoords.sparse(xPoint.vPos, rGlobalSparseCords),
                                         vSparseCoordsInternal),
                    xInternalEntiresHelper
                ) += 1;
        }, xPoints);

        // compute internal prefix sum
        for(size_t uiI = 0; uiI < D; uiI++)
        {
            red_entry_arr_t vRelevantSparseCoordsInternal = relevant(vSparseCoordsInternal, uiI);
            rSparseCoords.template iterate<D-1>([&]( const red_pos_t&, const red_pos_t& vTo )
            {
                pos_t vFullTo = expand(vTo, uiI);
                val_t uiPrefixSum = 0;
                rSparseCoords.iterate([&](coordinate_t, coordinate_t uiTo)
                {
                    vFullTo[uiI] = uiTo;
                    uiPrefixSum += rPrefixSums.get(vFullTo, xInternalEntiresHelper);
                    rPrefixSums.get(vFullTo, xInternalEntires) = uiPrefixSum;
                }, vPredecessors[uiI]->vSparseCoordsInternal[uiI]);
            }, vRelevantSparseCoordsInternal);
        }
        rPrefixSums.remove(xInternalEntiresHelper);

        // construct overlay sum grid
        for(size_t uiI = 0; uiI < D; uiI++)
            if(vPredecessors[uiI] != nullptr)
            {
                red_entry_arr_t vRelevantSparseCoordsOverlay = relevant(vSparseCoordsOverlay, uiI);

                red_pos_t vAxisSizes = rSparseCoords.axisSizes(vRelevantSparseCoordsOverlay);
                vOverlayEntries[uiI] = rPrefixSums.add(vAxisSizes);

                rSparseCoords.template iterate<D-1>([&]( const red_pos_t& vFrom, const red_pos_t& vTo )
                {
                    pos_t vFullFrom = expand(vFrom, uiI);
                    vFullFrom[uiI] = vPos[uiI];

                    rPrefixSums.get(vTo, vOverlayEntries[uiI]) = 
                        vPredecessors[uiI]->get(rSparseCoords, rPrefixSums, vFullFrom);
                }, vRelevantSparseCoordsOverlay);
            }
    }

    val_t get(const sparse_coord_t& rSparseCoords, const prefix_sum_grid_t& rPrefixSums, pos_t vCoords) const
    {
        val_t uiRet = 0;

        for(size_t uiI = 0; uiI < D; uiI++)
            uiRet += rPrefixSums.get(
                rSparseCoords.sparse(relevant(vCoords, uiI), relevant(vSparseCoordsOverlay, uiI)), 
                vOverlayEntries[uiI]);

        assert(uiRet >= (D-1) * uiBottomLeft);
        uiRet -= (D-1) * uiBottomLeft;

        uiRet += rPrefixSums.get(rSparseCoords.sparse(vCoords, vSparseCoordsInternal), xInternalEntires);

        return uiRet;
    }

    template<typename T>
    std::array<T, D-1> relevant( const std::array<T, D>& vAllEntries, size_t uiI) const
    {
        std::array<T, D-1> vRelevantEntries;
        for(size_t uiJ = 0; uiJ < uiI; uiJ++)
            vRelevantEntries[uiJ] = vAllEntries[uiJ];
        for(size_t uiJ = uiI+1; uiJ < D; uiJ++)
            vRelevantEntries[uiJ-1] = vAllEntries[uiJ];
        return vRelevantEntries;
    }

    template<typename T>
    std::array<T, D> expand( const std::array<T, D-1>& vCompressed, size_t uiI) const
    {
        std::array<T, D> vAllEntries;
        for(size_t uiJ = 0; uiJ < uiI; uiJ++)
            vAllEntries[uiJ] = vCompressed[uiJ];
        for(size_t uiJ = uiI+1; uiJ < D; uiJ++)
            vAllEntries[uiJ] = vCompressed[uiJ-1];
        return vAllEntries;
    }
};


} // namespace kdpstree


#if WITH_PYTHON
#endif
