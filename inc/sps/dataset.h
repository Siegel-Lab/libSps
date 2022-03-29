#pragma once

#include "sps/util.h"
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
    using desc_t = Desc<type_defs>;

    std::array<typename sparse_coord_t::Entry, D> vSparseCoords;
    typename overlay_grid_t::template Entry<D> xOverlays;

    class DimIterator
    {
        typename points_t::EntryIterator xIt;
        size_t uiDimension;
    public:
        DimIterator(typename points_t::EntryIterator xIt, size_t uiDimension) :
            xIt(xIt),
            uiDimension(uiDimension)
        {}

        void operator++()
        {
            ++xIt;
        }

        const coordinate_t operator*() const
        {
            return xIt->vPos[uiDimension];
        }

        bool operator!=(const DimIterator& rOther) const
        {
            return xIt != rOther.xIt;
        }
    
        friend std::ostream& operator<<( std::ostream& os, const DimIterator& rIt )
        {
            os << rIt.xIt;

            os << " d";
            os << rIt.uiDimension;

            return os;
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
            coordinate_t uiA = rDataset.overlayIndex(rOverlays, rSparseCoords, a.vPos);
            coordinate_t uiB = rDataset.overlayIndex(rOverlays, rSparseCoords, b.vPos);
            return (uiA == uiB && a.uiDescOffset < b.uiDescOffset) || uiA < uiB;
        }

        point_t min_value( ) const
        {
            // @todo undo sparsity
            return point_t(rOverlays.posOf(rDataset.xOverlays.uiStartIndex, rDataset.xOverlays), 0);
        };

        point_t max_value( ) const
        {
            // @todo undo sparsity
            auto vPos = rOverlays.posOf(rDataset.xOverlays.uiStartIndex + rOverlays.sizeOf(rDataset.xOverlays) - 1,
                                           rDataset.xOverlays);
            return point_t(vPos, std::numeric_limits<size_t>::max());
        };
    };
    sort_func_t<points_it_t, PointsComperator> sort_points = sort_func_t<points_it_t, PointsComperator>( );

  public:
    Dataset() 
    {}

    Dataset(overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords, prefix_sum_grid_t& rPrefixSums, 
            points_t& vPoints, typename points_t::Entry xPoints, std::optional<progress_stream_t> xProg)
    {
        // generate the overall sparse coordinates
        size_t uiNumPerDimension = (size_t)std::pow(xPoints.uiEndIndex - xPoints.uiStartIndex, 1.0/(D*D));
        xProg << "generating " << uiNumPerDimension << " overlays per dimension\n";
        pos_t vMaxCoord{};
        vPoints.iterate([&](const point_t& xPoint){
            for(size_t uiI = 0; uiI < D; uiI++)
                vMaxCoord[uiI] = std::max(vMaxCoord[uiI], xPoint.vPos[uiI]);
        }, xPoints);

        for(size_t uiI = 0; uiI < D; uiI++)
        {
            size_t uiBlockSize = std::max(1ul, vMaxCoord[uiI] / uiNumPerDimension);
            vPoints.sortByDim(uiI, xPoints);
            vSparseCoords[uiI] = rSparseCoords.add(
                DimIterator(vPoints.cbegin(xPoints), uiI), 
                DimIterator(vPoints.cend(xPoints), uiI), uiBlockSize);
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
            while(xCurrPoints.uiEndIndex < xPoints.uiEndIndex && 
                    overlayIndex(rOverlays, rSparseCoords, vPoints.vData[xCurrPoints.uiEndIndex].vPos) == uiI)
                ++xCurrPoints.uiEndIndex;
            xProg << "generating overlay for points " << xCurrPoints.uiStartIndex << " - " 
                  << xCurrPoints.uiEndIndex << " of " << xPoints.uiStartIndex << " - " << xPoints.uiEndIndex << "\n";

            // get bottom left position (compressed)
            pos_t vPos = rOverlays.posOf(uiI + xOverlays.uiStartIndex, xOverlays);

            xProg << "overlay anchor is " << vPos << "\n";


            // collect direct predecessor overlays for each dimension
            std::array<overlay_t*, D> vPredecessors;
            bool bAllLargerZero = true;
            for(size_t uiD = 0; uiD < D; uiD++)
                if(vPos[uiD] > 0)
                {
                    --vPos[uiD];
                    assert(rOverlays.indexOf(vPos, xOverlays) < uiI);
                    xProg << "predecessor dim " << uiD << " is " << rOverlays.indexOf(vPos, xOverlays) << "\n";
                    vPredecessors[uiD] = &rOverlays.get(vPos, xOverlays);
                    vPredecessors[uiD]->stream( std::cout, rSparseCoords,  rPrefixSums, vPoints ) << std::endl;
                    ++vPos[uiD];
                }
                else
                {
                    xProg << "predecessor dim " << uiD << " is nullptr\n";
                    vPredecessors[uiD] = nullptr;
                    bAllLargerZero = false;
                }
            
            // compute the prefix sum of the overlay to the bottom-left-front-...
            val_t uiBottomLeft = 0;
            if (bAllLargerZero)
            {
                for(size_t uiI = 0; uiI < D; uiI++)
                    --vPos[uiI];
                uiBottomLeft = rOverlays.get(vPos, xOverlays).get(rSparseCoords, rPrefixSums, vPos);
                for(size_t uiI = 0; uiI < D; uiI++)
                    ++vPos[uiI];
            }

            xProg << "generating overlay now...\n";
            // generate the overlay
            rOverlays.vData[xOverlays.uiStartIndex + uiI].generate(
                rSparseCoords,
                rPrefixSums,
                vPoints,
                xCurrPoints,
                vPredecessors,
                vPos,
                uiBottomLeft,
                vSparseCoords,
                xProg
            );

            std::cout << rOverlays.vData[xOverlays.uiStartIndex  + uiI] << std::endl;
            rOverlays.vData[xOverlays.uiStartIndex + uiI].stream( std::cout, rSparseCoords,  rPrefixSums, vPoints ) 
                << std::endl;

            // prepare for collecting the next set of points
            xCurrPoints.uiStartIndex = xCurrPoints.uiEndIndex;
        }
    }

    coordinate_t overlayIndex(overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords, const pos_t& vPos) const
    {
        return rOverlays.indexOf(rSparseCoords.sparse(vPos, vSparseCoords), xOverlays);
    }
    
    val_t get(const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords, 
              const prefix_sum_grid_t& rPrefixSums, const pos_t& vPos, 
              std::optional<progress_stream_t> xProg = { } ) const
    {
        auto vSparsePos = rSparseCoords.sparse(vPos, vSparseCoords);
        xProg << vPos << " -> " << vSparsePos << "\n";
        return rOverlays.get(vSparsePos, xOverlays).get( rSparseCoords, rPrefixSums, vPos, xProg );
    }
    
    friend std::ostream& operator<<( std::ostream& os, const Dataset& rDataset )
    {
        os << "<" << std::endl;
        os << "\tvSparseCoords: ";
        os << rDataset.vSparseCoords << std::endl;

        os << "\txOverlayGrid: ";
        os << rDataset.xOverlays << std::endl;
        os << ">";

        return os;
    }
    
    std::ostream& stream( std::ostream& os, const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
                          const prefix_sum_grid_t& rPrefixSums, const points_t& vPoints, const desc_t& vDesc ) const
    {
        os << "<" << std::endl;
        os << "\tvSparseCoords: ";
        for(size_t uiI = 0; uiI < D; uiI++)
            vSparseCoords[uiI].stream(os, rSparseCoords) << " ";
        os << std::endl;

        os << "\txOverlayGrid: ";
        xOverlays.stream(os, rOverlays, rSparseCoords, rPrefixSums, vPoints, vDesc) << std::endl;

        os << ">";

        return os;
    }
};


} // namespace sps
