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
            return rDataset.overlayCoord(rSparseCoords, xIt->vPos)[ uiDimension ];
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

#ifndef NDEBUG
            coordinate_t uiTest = rOverlays.indexOf( vPos, rDataset.xOverlays );
            assert( uiTest == uiIndex );
#endif

            pos_t vActPos = rSparseCoords.invSparse( vPos, rDataset.vSparseCoords );

#ifndef NDEBUG
            uiTest = rDataset.overlayIndex( rOverlays, rSparseCoords, vActPos );
            assert( uiTest == uiIndex );
#endif

            return point_t( vActPos, std::numeric_limits<size_t>::max( ) );
        };
    };
    sort_func_t<points_it_t, PointsBinComperator> sort_points_bin = sort_func_t<points_it_t, PointsBinComperator>( );

    struct PointsComperator
    {
        const Dataset& rDataset;
        sparse_coord_t& rSparseCoords;
        const size_t uiDim;

        PointsComperator( const Dataset& rDataset, sparse_coord_t& rSparseCoords, size_t uiDim )
            : rDataset( rDataset ), rSparseCoords( rSparseCoords ), uiDim( uiDim )
        {}

        bool operator( )( const point_t& a, const point_t& b ) const
        {
            return rDataset.overlayCoord(rSparseCoords, a.vPos)[uiDim]
                 < rDataset.overlayCoord(rSparseCoords, b.vPos)[uiDim];
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
    Dataset( )
    {}

    Dataset( overlay_grid_t& rOverlays, sparse_coord_t& rSparseCoords, prefix_sum_grid_t& rPrefixSums,
             points_t& vPoints, typename points_t::Entry xPoints, progress_stream_t xProg )
    {
        // generate the overall sparse coordinates 
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            size_t uiNumCoords = 0;
            coordinate_t uiLast = std::numeric_limits<coordinate_t>::max();
            sort_points(vPoints.vData.begin( ) + xPoints.uiStartIndex, vPoints.vData.begin( ) + xPoints.uiEndIndex,
                        PointsComperator(*this, rSparseCoords, uiI));
            vPoints.iterate(
                [ & ]( const point_t& xPoint ) {
                    coordinate_t uiCurr = overlayCoord(rSparseCoords, xPoint.vPos)[ uiI ];
                    if(uiCurr != uiLast)
                    {
                        uiLast = uiCurr;
                        ++uiNumCoords;
                    }
                },
                xPoints );
            size_t uiNumPerDimension = (size_t)std::pow( uiNumCoords, 1.0 / (float)( D ) ); // 
            xProg << Verbosity( 0 ) << "generating " << uiNumPerDimension << " overlays in dimension " << uiI 
                                    << " for " << uiNumCoords << " different coordinates.\n";
            size_t uiBlockSize = std::max( 1ul, uiNumCoords / uiNumPerDimension );
            vSparseCoords[ uiI ] = rSparseCoords.add( DimIterator( *this, rSparseCoords, vPoints.cbegin( xPoints ), 
                                                                        vPoints.cend( xPoints ), uiI, uiBlockSize ),
                                                      DimIterator( *this, rSparseCoords, vPoints.cend( xPoints ), 
                                                                        vPoints.cend( xPoints ), uiI, uiBlockSize ) );
        }
        // generate overlay grid
        xOverlays = rOverlays.add( rSparseCoords.axisSizes( vSparseCoords ) );

        // sort points so that they match the overlay grid order
        sort_points_bin( vPoints.vData.begin( ) + xPoints.uiStartIndex, vPoints.vData.begin( ) + xPoints.uiEndIndex,
                         PointsBinComperator( *this, rOverlays, rSparseCoords ) );

        // generate all overlays

        coordinate_t uiCenterBin = rSparseCoords.replace( std::numeric_limits<coordinate_t>::max() / 2, 
                                                          vSparseCoords[1] );

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
            pos_t vPosBottomLeft = rOverlays.posOf( uiI + xOverlays.uiStartIndex, xOverlays );
            pos_t vPosBottomLeftActual = actualFromGridPos(rSparseCoords, vPosBottomLeft);
            pos_t vPosTopRightActual = actualTopRightFromGridPos(rSparseCoords, vPosBottomLeft);

            xProg << Verbosity( 2 ) << "overlay anchor is " << vPosBottomLeft 
                  << " actual pos is " << vPosBottomLeftActual << "\n";

            bool bRightOfDiagonal = vPosBottomLeftActual[0] >= vPosBottomLeftActual[1];
            bool bLeftOfDiagonal = vPosBottomLeftActual[1] >= vPosBottomLeftActual[0]; 

            // collect direct predecessor overlays for each dimension
            std::array<std::vector<coordinate_t>, D> vPredecessors;
            for( size_t uiD = 0; uiD < D; uiD++ )
            {
                if( vPosBottomLeftActual[ uiD ] > 0 )
                {
                    // there can be multiple predecessors
                    if( COORD_TRANSFROM && ( (uiD == 0 && bLeftOfDiagonal) || (uiD == 1 && bRightOfDiagonal) ) ) 
                    {
                        pos_t vItrPos = vPosBottomLeft;
                        assert(vItrPos[ 0 ] > 0);
                        --vItrPos[ 0 ];

                        if(vPosBottomLeftActual[1] == vPosBottomLeftActual[0] && uiD == 0)
                        {
                            ++vItrPos[ 1 ];
                            if(vItrPos[1] == uiCenterBin)
                                vItrPos[1] = std::numeric_limits<coordinate_t>::max();
                        }

                        while(true)
                        {
                            pos_t vItrPosActual = actualFromGridPos(rSparseCoords, vItrPos);
                            if(vItrPosActual[uiD] >= vPosTopRightActual[uiD])
                                break;
                            
                            pos_t vPosTopRightActual = actualTopRightFromGridPos(rSparseCoords, vItrPos);
                            if(vPosTopRightActual[uiD] > vPosBottomLeftActual[uiD])
                                vPredecessors[ uiD ].push_back( 
                                    overlayIndex( rOverlays, rSparseCoords, vItrPosActual ) );

                            ++vItrPos[1];
                            if(vItrPos[1] == uiCenterBin)
                                vItrPos[1] = std::numeric_limits<coordinate_t>::max();
                        }
                    }
                    // there can only be one predecessor
                    else 
                    {
                        --vPosBottomLeftActual[ uiD ];

                        vPredecessors[ uiD ].push_back(overlayIndex(rOverlays, rSparseCoords, vPosBottomLeftActual ));

                        ++vPosBottomLeftActual[ uiD ];
                    }
                }
                else
                    xProg << Verbosity( 2 ) << "predecessor dim " << uiD << " nonexistant\n";
                

                for( size_t uiJ = 0; uiJ < vPredecessors[ uiD ].size(); uiJ++ )
                {
                    xProg << Verbosity( 2 ) << "predecessor " << uiJ << " dim " << uiD << " is "
                            << vPredecessors[ uiD ][ uiJ ] << "\n";
                    assert( vPredecessors[ uiD ][ uiJ ] < uiI );

                    if( xProg.active( ) )
                        rOverlays.vData[ vPredecessors[ uiD ][ uiJ ] ].stream(
                                std::cout, rSparseCoords, rPrefixSums, vPoints ) << std::endl;
                }
            }

            xProg << Verbosity( 1 ) << "generating overlay ouf of " 
                  << xCurrPoints.uiEndIndex - xCurrPoints.uiStartIndex << " points now...\n";
            // generate the overlay
            // @todo generate needs the new parameters & each overlay only needs the sparse coords that are between its beginning and end pos
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

    pos_t overlayCoord(const sparse_coord_t& rSparseCoords, pos_t vPos) const
    {
        if constexpr(!COORD_TRANSFROM)
            return vPos;
            
        coordinate_t uiX = vPos[0];
        coordinate_t uiY = vPos[1];
        vPos[0] = std::min(uiX, uiY);
        vPos[1] = std::numeric_limits<coordinate_t>::max() / 2;
        // @note: if vSparseCoords[0] is not already defined uiBinStartPos becomes
        //        std::numeric_limits<coordinate_t>::max(), and then some awkward number due to the if's below
        coordinate_t uiBinStartPos = rSparseCoords.invReplace(
                                                    rSparseCoords.replace( vPos[0], vSparseCoords[0] ),
                                                    vSparseCoords[0]
                                                );
        if(uiX > uiY + uiBinStartPos)
        {
            assert( uiX - (uiY + uiBinStartPos) <= vPos[1] ); // assert that there are no overflows
            vPos[1] += uiX - (uiY + uiBinStartPos);
            assert(vPos[1] > std::numeric_limits<coordinate_t>::max() / 2);
        }
        else if(uiY > uiX + uiBinStartPos)
        {
            assert( vPos[1] >= uiY - (uiX + uiBinStartPos) ); // assert that there are no overflows
            vPos[1] -= uiY - (uiX + uiBinStartPos);
            assert(vPos[1] < std::numeric_limits<coordinate_t>::max() / 2);
        }

        return vPos;
    }

    pos_t overlayCoordInv(const sparse_coord_t& /*rSparseCoords*/, pos_t vPos) const
    {
        if constexpr(!COORD_TRANSFROM)
            return vPos;

        coordinate_t uiX = vPos[0];
        coordinate_t uiY = vPos[1];
        if(uiY >= std::numeric_limits<coordinate_t>::max() / 2)
        {
            vPos[1] = uiX;
            assert( uiX + uiY >= std::numeric_limits<coordinate_t>::max() / 2 );
            if(uiY == std::numeric_limits<coordinate_t>::max())
                vPos[0] = std::numeric_limits<coordinate_t>::max();
            else
                vPos[0] = uiX + uiY - std::numeric_limits<coordinate_t>::max() / 2;
        }
        else
        {
            vPos[0] = uiX;
            assert( uiX + std::numeric_limits<coordinate_t>::max() / 2 >= uiY );
            if(uiY == 0)
                vPos[1] = std::numeric_limits<coordinate_t>::max();
            else
                vPos[1] = uiX + std::numeric_limits<coordinate_t>::max() / 2 - uiY;
        }
        return vPos;
    }

    pos_t fixOverlayOrder(const sparse_coord_t& rSparseCoords, pos_t vPos) const
    {
        if constexpr(!COORD_TRANSFROM)
            return vPos;

        coordinate_t uiCenterBin = rSparseCoords.replace( std::numeric_limits<coordinate_t>::max() / 2, 
                                                          vSparseCoords[1] );

        if(vPos[1] < uiCenterBin)
        {
            coordinate_t uiNumBins = rSparseCoords.replace( vSparseCoords[1].uiEndCord, vSparseCoords[1] ) + 1;
            vPos[1] = uiNumBins - (vPos[1] + 1);
            assert(vPos[1] < uiCenterBin);
        }
        else
            vPos[1] -= uiCenterBin;

        return vPos;
    }

    pos_t fixOverlayOrderInv(const sparse_coord_t& rSparseCoords, pos_t vPos) const
    {
        if constexpr(!COORD_TRANSFROM)
            return vPos;

        coordinate_t uiCenterBin = rSparseCoords.replace( std::numeric_limits<coordinate_t>::max() / 2, 
                                                          vSparseCoords[1] );
        coordinate_t uiNumBins = rSparseCoords.replace( vSparseCoords[1].uiEndCord, vSparseCoords[1] ) + 1;

        if(vPos[1] < uiCenterBin)
            vPos[1] += uiCenterBin;
        else if(vPos[1] == uiNumBins)
            vPos[1] = 0;
        else if(vPos[1] != std::numeric_limits<coordinate_t>::max())
        {
            assert(uiNumBins >= (vPos[1] + 1));
            vPos[1] = uiNumBins - (vPos[1] + 1);
            assert(vPos[1] < uiCenterBin);
        }
        // else if(vPos[1] == std::numeric_limits<coordinate_t>::max())
        //      keep vPos[1] the same

        return vPos;
    }

    pos_t actualFromGridPos(const sparse_coord_t& rSparseCoords, pos_t vPos) const
    {
        return overlayCoordInv(
                rSparseCoords,
                fixOverlayOrderInv( 
                    rSparseCoords,
                    rSparseCoords.invSparse( vPos , vSparseCoords )
                )
            );
    }
    pos_t actualTopRightFromGridPos(const sparse_coord_t& rSparseCoords, pos_t vPos) const
    {
        ++vPos[0];

        coordinate_t uiCenterBin = rSparseCoords.replace( std::numeric_limits<coordinate_t>::max() / 2,
                                                            vSparseCoords[1] );
        if(vPos[1] == 0)
        {
            coordinate_t uiNumBins = rSparseCoords.replace( vSparseCoords[1].uiEndCord, vSparseCoords[1] ) + 1;
            vPos[1] = uiNumBins - uiCenterBin;
        }
        else
        {
            ++vPos[1];
            // if we now have the center position we actually went out on the right side of the original coordinates
            if(vPos[1] == uiCenterBin)
                vPos[1] = std::numeric_limits<coordinate_t>::max();
        }
        
        for(size_t uiI = 2; uiI < D; uiI++)
            ++vPos[uiI];

        return actualFromGridPos(rSparseCoords, vPos);
    }

    coordinate_t overlayIndex( overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords, const pos_t& vPos ) const
    {
        auto vPosSparse = rSparseCoords.sparse(
                    fixOverlayOrder(rSparseCoords, overlayCoord(rSparseCoords, vPos)), 
                    vSparseCoords
                );
        return rOverlays.indexOf( vPosSparse, xOverlays );
    }

    val_t get( const overlay_grid_t& rOverlays, const sparse_coord_t& rSparseCoords,
               const prefix_sum_grid_t& rPrefixSums, const pos_t& vPos, progress_stream_t& xProg ) const
    {
        auto vSparsePos = rSparseCoords.sparse( 
                fixOverlayOrder( rSparseCoords, overlayCoord(rSparseCoords, vPos) ),
                vSparseCoords
            );
        xProg << Verbosity( 2 );
        if( xProg.active( ) )
            xProg << "\t" << vPos << " -> " << vSparsePos << "; that's overlay "
                  << rOverlays.indexOf( vSparsePos, xOverlays ) << "\n";
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

        os << "\txOverlayGrid: ";
        xOverlays.stream( os, rOverlays, rSparseCoords, rPrefixSums, vPoints, vDesc ) << std::endl;

        os << ">";

        return os;
    }
};


} // namespace sps