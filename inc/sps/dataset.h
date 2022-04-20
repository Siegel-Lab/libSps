#pragma once

#include "sps/desc.h"
#include "sps/overlay.h"
#include "sps/point.h"
#include "sps/sparse_coordinate.h"
#include "sps/type_defs.h"
#include "sps/util.h"
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

    typename sparse_coord_t::EntryArray xSparseCoordsDependantDimension;

    typename overlay_grid_t::template Entry<D> xOverlays;

    class DimIterator
    {
        const Dataset& rDataset;
        sparse_coord_t& rSparseCoords;
        typename points_t::EntryIterator xIt, xItEnd;
        size_t uiDimension;
        size_t uiBlockSize;

      public:
        DimIterator( const Dataset& rDataset, sparse_coord_t& rSparseCoords, typename points_t::EntryIterator xIt, 
                     typename points_t::EntryIterator xItEnd, size_t uiDimension,
                     size_t uiBlockSize ) : 
                rDataset( rDataset ), rSparseCoords( rSparseCoords ), xIt( xIt ), xItEnd(xItEnd), 
                uiDimension( uiDimension ), uiBlockSize(uiBlockSize)
        {}

        void operator++( )
        {
            for(size_t uiI = 0; uiI < uiBlockSize && xIt != xItEnd; uiI++)
                ++xIt;
        }

        const coordinate_t operator*( ) const
        {
            return xIt->vPos[ uiDimension ];
        }

        bool operator!=( const DimIterator& rOther ) const
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
    struct PointsBinComperator
    {
        const Dataset& rDataset;
        overlay_grid_t& rOverlays;
        sparse_coord_t& rSparseCoords;

        PointsBinComperator( const Dataset& rDataset, overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords )
            : rDataset( rDataset ),
              rOverlays( rOverlays ),
              rSparseCoords( rSparseCoords )
        {}

        bool operator( )( const point_t& a, const point_t& b ) const
        {
            coordinate_t uiA = rDataset.overlayIndex( rOverlays, rSparseCoords, a.vPos );
            coordinate_t uiB = rDataset.overlayIndex( rOverlays, rSparseCoords, b.vPos );
            if( uiA == uiB )
                return a.uiDescOffset < b.uiDescOffset;
            return uiA < uiB;
        }

        point_t min_value( ) const
        {
            return point_t( rSparseCoords.invSparse(
                                rOverlays.posOf( rDataset.xOverlays.uiStartIndex, rDataset.xOverlays ), 
                                rDataset.vSparseCoords ),
                            0 );
        };

        point_t max_value( ) const
        {
            coordinate_t uiIndex = rDataset.xOverlays.uiStartIndex + rOverlays.sizeOf( rDataset.xOverlays ) - 1;

            pos_t vPos = rOverlays.posOf( uiIndex, rDataset.xOverlays );
            pos_t vActPos = rDataset.actualFromGridPos(rSparseCoords, vPos);

            return point_t( vActPos, std::numeric_limits<size_t>::max( ) );
        };
    };
    sort_func_t<points_it_t, PointsBinComperator> sort_points_bin = sort_func_t<points_it_t, PointsBinComperator>( );

    struct PointsComperator
    {
        const size_t uiDim;

        PointsComperator( size_t uiDim )
            : uiDim( uiDim )
        {}

        bool operator( )( const point_t& a, const point_t& b ) const
        {
            return a.vPos[uiDim] < b.vPos[uiDim];
        }

        point_t min_value( ) const
        {
            point_t xRet{ };
            return xRet;
        };

        point_t max_value( ) const
        {
            point_t xRet{ };
            xRet.vPos[ uiDim ] = std::numeric_limits<coordinate_t>::max( );
            return xRet;
        };
    };
    sort_func_t<points_it_t, PointsComperator> sort_points = sort_func_t<points_it_t, PointsComperator>( );

  public:
    Dataset( ) : xSparseCoordsDependantDimension()
    {}

    typename sparse_coord_t::Entry makeSparseCoords(
                                            sparse_coord_t& rSparseCoords, points_t& vPoints,
                                            typename points_t::Entry xPoints, size_t uiDim, progress_stream_t xProg,
                                            coordinate_t uiFixedStart = std::numeric_limits<coordinate_t>::max(), 
                                            coordinate_t uiFixedEnd = std::numeric_limits<coordinate_t>::max()) const
    {
        size_t uiNumCoords = 0;
        coordinate_t uiLast = std::numeric_limits<coordinate_t>::max();
        sort_points(vPoints.vData.begin( ) + xPoints.uiStartIndex, vPoints.vData.begin( ) + xPoints.uiEndIndex,
                    PointsComperator(uiDim));
        vPoints.iterate(
            [ & ]( const point_t& xPoint ) {
                coordinate_t uiCurr = xPoint.vPos[ uiDim ];
                if(uiCurr != uiLast)
                {
                    uiLast = uiCurr;
                    ++uiNumCoords;
                }
            },
            xPoints );
        size_t uiNumPerDimension = (size_t)std::pow( uiNumCoords, 1.0 / (float)( D ) ); // 
        xProg << "generating " << uiNumPerDimension << " overlays in dimension " << uiDim 
                                << " for " << uiNumCoords << " different coordinates. Points: " << xPoints << "\n";
        size_t uiBlockSize = std::max( 1ul, uiNumCoords / uiNumPerDimension );
        auto xStart = DimIterator( *this, rSparseCoords, vPoints.cbegin( xPoints ), 
                                    vPoints.cend( xPoints ), uiDim, uiBlockSize );
        auto xEnd = DimIterator( *this, rSparseCoords, vPoints.cend( xPoints ), 
                                 vPoints.cend( xPoints ), uiDim, uiBlockSize );
        if(uiFixedStart == std::numeric_limits<coordinate_t>::max())
            return rSparseCoords.add(xStart, xEnd);
        return rSparseCoords.addStartEnd(xStart, xEnd, uiFixedStart, uiFixedEnd);
    }

#pragma GCC diagnostic push
// vPosTopRightActual not used with DEPENDANT_DIMENSION == false
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
    std::array<std::vector<coordinate_t>, D> getPredecessor(
            const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords, 
            pos_t vGridPos, pos_t vPosBottomLeftActual, pos_t vPosTopRightActual, progress_stream_t xProg) const
    {

        auto vbHasPredecessor = hasPredecessor(rSparseCoords, vGridPos);
        std::array<std::vector<coordinate_t>, D> vPredecessors {};
        for( size_t uiD = 0; uiD < D; uiD++ )
        {
            if( vbHasPredecessor[ uiD ] )
            {
                if constexpr(DEPENDANT_DIMENSION)
                {
#ifndef NDEBUG
                    coordinate_t uiIdx = rOverlays.indexOf(vGridPos, xOverlays);
#endif
                    if(uiD != 1)
                    {
                        pos_t vItrPos = vPosBottomLeftActual;
                        assert(vItrPos[uiD] > 0);
                        // move outside the current overlay
                        --vItrPos[uiD];
                        pos_t vItrTopRight;
                        do
                        {
                            coordinate_t uiItrIndex = overlayIndex(rOverlays, rSparseCoords, vItrPos);
                            vPredecessors[ uiD ].push_back( uiItrIndex );
                            pos_t vItrGridPos = rOverlays.posOf(uiItrIndex, xOverlays);
                            vItrTopRight = actualTopRightFromGridPos(rSparseCoords, vItrGridPos);
                            
                            xProg << Verbosity( 2 ) << "uiD " << uiD << " vItrPos " << vItrPos << " uiItrIndex "
                                << uiItrIndex << " vItrGridPos " << vItrGridPos << " vItrTopRight " << vItrTopRight 
                                << " vPosBottomLeftActual " << vPosBottomLeftActual
                                << " vPosTopRightActual " << vPosTopRightActual << "\n";
                            assert(uiItrIndex != uiIdx);

                            // move upwards
                            assert(vItrTopRight[1] > vItrPos[1] || 
                                   vItrTopRight[1] == std::numeric_limits<coordinate_t>::max());
                            vItrPos[1] = vItrTopRight[1];
                        }
                        while(vItrTopRight[1] < vPosTopRightActual[1] && 
                              vItrTopRight[1] != std::numeric_limits<coordinate_t>::max());

                        continue;
                    }
                }

                --vPosBottomLeftActual[ uiD ];

                vPredecessors[ uiD ].push_back(overlayIndex(rOverlays, rSparseCoords, vPosBottomLeftActual ));

                xProg << Verbosity( 3 ) << "predecessor dim " << uiD << " pos " 
                        << vPosBottomLeftActual << " index " 
                        << vPredecessors[ uiD ].back() << "\n";

                ++vPosBottomLeftActual[ uiD ];
            }
            else
                xProg << Verbosity( 2 ) << "predecessor dim " << uiD << " nonexistant\n";
            
        }
        return vPredecessors;
    }
#pragma GCC diagnostic pop

    Dataset( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords, prefix_sum_grid_t& rPrefixSums,
             points_t& vPoints, typename points_t::Entry xPoints, progress_stream_t xProg ) : vSparseCoords()
    {
        // generate the overall sparse coordinates 
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            xProg << Verbosity( 0 );
            if constexpr(DEPENDANT_DIMENSION)
            {
                if(uiI == 1)
                {
                    coordinate_t uiAbsoluteStart = std::numeric_limits<coordinate_t>::max();
                    coordinate_t uiAbsoluteEnd = 0;
                    vPoints.iterate(
                        [ & ]( const point_t& xPoint ) {
                            uiAbsoluteStart = std::min(uiAbsoluteStart, xPoint.vPos[1]);
                            uiAbsoluteEnd = std::max(uiAbsoluteEnd, xPoint.vPos[1]);
                        },
                        xPoints );

                    coordinate_t uiNumBins = rSparseCoords.axisSize(vSparseCoords[0]);
                    typename points_t::Entry xCurrPoints{ };
                    xCurrPoints.uiEndIndex = xPoints.uiStartIndex;
                    for(size_t uiBin = 0; uiBin < uiNumBins; uiBin++)
                    {
                        // collect points for overlay uiI (points still sorted by dim 0)
                        xCurrPoints.uiStartIndex = xCurrPoints.uiEndIndex;
                        while( xCurrPoints.uiEndIndex < xPoints.uiEndIndex &&
                               rSparseCoords.replace( vPoints.vData[ xCurrPoints.uiEndIndex ].vPos[0],
                                                      vSparseCoords[0] ) == uiBin )
                            ++xCurrPoints.uiEndIndex;

                        auto xCurr = makeSparseCoords(rSparseCoords, vPoints, xCurrPoints, 1, xProg,
                                                                uiAbsoluteStart, uiAbsoluteEnd);
                        if(vSparseCoords[ uiI ].uiStartIndex == std::numeric_limits<coordinate_t>::max() || 
                           rSparseCoords.axisSize(vSparseCoords[ uiI ]) < rSparseCoords.axisSize(xCurr) )
                            vSparseCoords[ uiI ] = xCurr;
                        sparse_coord_t::append(xSparseCoordsDependantDimension, xCurr);
                        xProg << Verbosity( 1 );
                    }
                    assert(xCurrPoints.uiEndIndex == xPoints.uiEndIndex);

                    continue;
                }
            }
            vSparseCoords[ uiI ] = makeSparseCoords(rSparseCoords, vPoints, xPoints, uiI, xProg);
        }
        xProg << Verbosity( 2 );
        if(xProg.active())
        {
            vSparseCoords[0].stream(std::cout << "sparse coords 0: ", rSparseCoords) << std::endl;
            if constexpr(DEPENDANT_DIMENSION)
                xSparseCoordsDependantDimension.stream(std::cout << "sparse coords 1: ", rSparseCoords) << std::endl;
            else
                vSparseCoords[1].stream(std::cout << "sparse coords 1: ", rSparseCoords) << std::endl;
            for( size_t uiI = 2; uiI < D; uiI++ )
                vSparseCoords[uiI].stream(std::cout << "sparse coords " << uiI << ": ", rSparseCoords) << std::endl;
        }


        // generate overlay grid
        xOverlays = rOverlays.add( rSparseCoords.axisSizes( vSparseCoords ) );

        // sort points so that they match the overlay grid order
        sort_points_bin( vPoints.vData.begin( ) + xPoints.uiStartIndex, vPoints.vData.begin( ) + xPoints.uiEndIndex,
                         PointsBinComperator( *this, rOverlays, rSparseCoords ) );

        // generate all overlays

        typename points_t::Entry xCurrPoints{ };
        xCurrPoints.uiStartIndex = xPoints.uiStartIndex;
        xCurrPoints.uiEndIndex = xPoints.uiStartIndex;
        for( coordinate_t uiI = 0; uiI < rOverlays.sizeOf( xOverlays ); uiI++ )
        {
            // collect points for overlay uiI
            while( xCurrPoints.uiEndIndex < xPoints.uiEndIndex &&
                   overlayIndex( rOverlays, rSparseCoords, vPoints.vData[ xCurrPoints.uiEndIndex ].vPos ) == uiI )
                ++xCurrPoints.uiEndIndex;
            xProg << Verbosity( 1 ) << "generating overlay for points " << xCurrPoints.uiStartIndex << " - "
                  << xCurrPoints.uiEndIndex << " of " << xPoints.uiStartIndex << " - " << xPoints.uiEndIndex << "\n";

            // get bottom left position (compressed)
            pos_t vGridPos = rOverlays.posOf( uiI + xOverlays.uiStartIndex, xOverlays );
            if constexpr(DEPENDANT_DIMENSION)
                if(!exists(rSparseCoords, vGridPos))
                {
                    xProg << Verbosity( 2 ) << "skipping nonexistant overlap\n";
                    continue;
                }
            pos_t vPosBottomLeftActual = actualFromGridPos(rSparseCoords, vGridPos);
            pos_t vPosTopRightActual = actualTopRightFromGridPos(rSparseCoords, vGridPos);

            xProg << Verbosity( 2 ) << "overlay anchor is " << vGridPos 
                  << " actual pos is " << vPosBottomLeftActual << "\n";

            // collect direct predecessor overlays for each dimension
            std::array<std::vector<coordinate_t>, D> vPredecessors = getPredecessor(
                rOverlays,
                rSparseCoords,
                vGridPos,
                vPosBottomLeftActual,
                vPosTopRightActual,
                xProg);

            for( size_t uiD = 0; uiD < D; uiD++ )
                for( size_t uiJ = 0; uiJ < vPredecessors[ uiD ].size(); uiJ++ )
                {
                    xProg << Verbosity( 2 ) << "predecessor " << uiJ << " dim " << uiD << " is "
                            << vPredecessors[ uiD ][ uiJ ] << "\n";
                    assert( vPredecessors[ uiD ][ uiJ ] < uiI );

                    if( xProg.active( ) )
                        rOverlays.vData[ vPredecessors[ uiD ][ uiJ ] ].stream(
                                std::cout, rSparseCoords, rPrefixSums, vPoints ) << std::endl;
                }

            xProg << Verbosity( 1 ) << "generating overlay ouf of " 
                  << xCurrPoints.uiEndIndex - xCurrPoints.uiStartIndex << " points now...\n";
            // generate the overlay
            rOverlays.vData[ xOverlays.uiStartIndex + uiI ].generate(
                rOverlays, rSparseCoords, rPrefixSums, vPoints, xCurrPoints, vPredecessors, 
                vPosBottomLeftActual, vPosTopRightActual, this, xProg );

            xProg << Verbosity( 2 );
            if( xProg.active( ) )
            {
                std::cout << rOverlays.vData[ xOverlays.uiStartIndex + uiI ] << std::endl;
                rOverlays.vData[ xOverlays.uiStartIndex + uiI ].stream( std::cout, rSparseCoords, rPrefixSums, vPoints )
                    << std::endl;
            }

            // prepare for collecting the next set of points
            xCurrPoints.uiStartIndex = xCurrPoints.uiEndIndex;
            coordinate_t uiNumDone = xCurrPoints.uiStartIndex - xPoints.uiStartIndex;
            coordinate_t uiNumTotal = xPoints.uiEndIndex - xPoints.uiStartIndex;
            if( xProg.printAgain( ) )
                xProg << Verbosity( 0 ) << uiNumDone << " out of " << uiNumTotal << ", thats "
                      << 100 * uiNumDone / (double)uiNumTotal << "%.\n";
        }
    }

    std::array<bool, D> hasPredecessor(const sparse_coord_t& /*rSparseCoords*/, pos_t vPos) const
    {
        std::array<bool, D> vRet;

        for(size_t uiI = 0; uiI < D; uiI++)
            vRet[uiI] = vPos[uiI] > 0;

        return vRet;
    }

    bool exists(const sparse_coord_t& rSparseCoords, pos_t vPos) const
    {
        if constexpr(DEPENDANT_DIMENSION)
            return rSparseCoords.axisSize( sparse_coord_t::at(xSparseCoordsDependantDimension, vPos[0] ) ) > vPos[1];
        return true;
    }

    pos_t actualFromGridPos(const sparse_coord_t& rSparseCoords, pos_t vPos) const
    {
        pos_t vRet = rSparseCoords.invSparse( vPos , vSparseCoords );
    
        // fix dependant dimension
        if constexpr(DEPENDANT_DIMENSION)
            vRet[1] = rSparseCoords.invReplace(vPos[1], sparse_coord_t::at(xSparseCoordsDependantDimension, 
                                    std::min(vPos[0], rSparseCoords.axisSize(vSparseCoords[0])-1)));

        return vRet;
    }

    pos_t actualTopRightFromGridPos(const sparse_coord_t& rSparseCoords, pos_t vPos) const
    {
        for(size_t uiI = 0; uiI < D; uiI++)
            ++vPos[uiI];

        pos_t vRet = actualFromGridPos(rSparseCoords, vPos);

        // fix dependant dimension
        if constexpr(DEPENDANT_DIMENSION)
            vRet[1] = rSparseCoords.invReplace(vPos[1],
                                               sparse_coord_t::at(xSparseCoordsDependantDimension, vPos[0]-1));
        
        return vRet;
    }

    pos_t overlayCoord(const sparse_coord_t& rSparseCoords, const pos_t& vPos) const
    {
        pos_t vRet = rSparseCoords.sparse(vPos, vSparseCoords);

        // fix dependant dimension
        if constexpr(DEPENDANT_DIMENSION)
            if(vRet[0] != std::numeric_limits<coordinate_t>::max())
                vRet[1] = rSparseCoords.replace(vPos[1], sparse_coord_t::at(xSparseCoordsDependantDimension, vRet[0]));

        return vRet;
    }

    struct OverlayInfo
    {
        pos_t vBottomLeft, vTopRight;
        pos_t vGridPos;
        coordinate_t uiIdx;
        std::array<std::vector<coordinate_t>, D> vPredIds;
        std::vector<pos_t> vvPoints{};
    };

    std::vector<OverlayInfo> getOverlayInfo(const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords, 
                                            const points_t& vPoints) const
    {
        std::vector<OverlayInfo> vRet;

        for(coordinate_t uiI = 0; uiI < rOverlays.sizeOf(xOverlays); uiI++)
        {
            vRet.emplace_back();
            vRet.back().uiIdx = uiI;
            vRet.back().vGridPos = rOverlays.posOf(uiI, xOverlays);
            vRet.back().vBottomLeft = actualFromGridPos(rSparseCoords, vRet.back().vGridPos);
            vRet.back().vTopRight = actualTopRightFromGridPos(rSparseCoords, vRet.back().vGridPos);
            vRet.back().vPredIds = getPredecessor(rOverlays, rSparseCoords, 
                                                  vRet.back().vGridPos,
                                                  vRet.back().vBottomLeft,
                                                  vRet.back().vTopRight,
                                                  typename type_defs::progress_stream_t( 0 )
                                                );
            vPoints.iterate([&](const point_t& xP){
                vRet.back().vvPoints.push_back(xP.vPos);
            }, rOverlays.vData[uiI].xPoints);
        }

        return vRet;
    }

    coordinate_t overlayIndex( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords, 
                               const pos_t& vPos ) const
    {
        return rOverlays.indexOf( overlayCoord(rSparseCoords, vPos), xOverlays );
    }

    val_t get( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
               const prefix_sum_grid_t& rPrefixSums, const pos_t& vPos, progress_stream_t& xProg ) const
    {
        auto vSparsePos = overlayCoord(rSparseCoords, vPos);
        xProg << Verbosity( 2 );
        if( xProg.active( ) )
            xProg << "\t" << vPos << " -> " << vSparsePos << "; that's overlay "
                  << rOverlays.indexOf( vSparsePos, xOverlays ) << "\n";
        for(size_t uiI = 0; uiI < D; uiI++)
            if(vSparsePos[uiI] == std::numeric_limits<coordinate_t>::max())
                return 0;
        return rOverlays.get( vSparsePos, xOverlays ).get( rSparseCoords, rPrefixSums, vPos, 
                                                           actualFromGridPos(rSparseCoords, vSparsePos), xProg );
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
        for( size_t uiI = 0; uiI < D; uiI++ )
            vSparseCoords[ uiI ].stream( os, rSparseCoords ) << " ";
        os << std::endl;
        if constexpr(DEPENDANT_DIMENSION)
        {
            os << "\txSparseCoordsDependantDimension: ";
            xSparseCoordsDependantDimension.stream( os, rSparseCoords ) << " ";
            os << std::endl;
        }

        os << "\txOverlayGrid: ";
        xOverlays.stream( os, rOverlays, rSparseCoords, rPrefixSums, vPoints, vDesc ) << std::endl;

        os << ">";

        return os;
    }
};


} // namespace sps