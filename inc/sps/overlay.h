#pragma once

#include "sps/nd_grid.h"
#include "sps/points.h"
#include "sps/sparse_coordinate.h"
#include <cassert>
#include <functional>
#include <string>



namespace sps
{

template <typename type_defs> class Dataset;

template <typename type_defs> class Overlay
{
    EXTRACT_TYPE_DEFS; // macro call

    using sparse_coord_t = SparseCoord<type_defs>;
    using prefix_sum_grid_t = NDGrid<type_defs, val_t>;
    using overlay_grid_t = NDGrid<type_defs, AlignedPower2<Overlay>>;


    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;
    using desc_t = Desc<type_defs>;
    using dataset_t = Dataset<type_defs>;

    using cord_it_t = typename sparse_coord_t::EntryIterator;
    using point_it_t = typename points_t::EntryIterator;

    using red_pos_t = std::array<coordinate_t, D - 1>;

    using entry_arr_t = std::array<typename sparse_coord_t::Entry, D>;
    using red_entry_arr_t = std::array<typename sparse_coord_t::Entry, D - 1>;

    class MergeIterator
    {
        std::vector<cord_it_t> vBegin;
        std::vector<cord_it_t> vEnd;
        coordinate_t uiEnd;

      public:
        MergeIterator( std::vector<cord_it_t> vBegin, std::vector<cord_it_t> vEnd, coordinate_t uiEnd )
            : vBegin( vBegin ), vEnd( vEnd ), uiEnd(uiEnd)
        {}

        size_t getSmallestValid() const
        {
            size_t uiSmallestValid = 0;
            while(uiSmallestValid < vBegin.size() && !(vBegin[uiSmallestValid] != vEnd[uiSmallestValid]))
                ++uiSmallestValid;
            assert(uiSmallestValid < vBegin.size());
            for(size_t uiI = uiSmallestValid + 1; uiI < vBegin.size(); uiI++)
                if(vBegin[uiI] != vEnd[uiI] && *vBegin[uiI] < *vBegin[uiSmallestValid])
                    uiSmallestValid = uiI;
            return uiSmallestValid;
        }

        void operator++( )
        {
            size_t uiSmallestValid = getSmallestValid();
            
            for(size_t uiI = 0; uiI < vBegin.size(); uiI++)
                if( uiI != uiSmallestValid && vBegin[uiI] != vEnd[uiI] && !(*vBegin[uiSmallestValid] < *vBegin[uiI]) )
                    ++vBegin[uiI];
            ++vBegin[uiSmallestValid];
        }

        const coordinate_t operator*( ) const
        {
            return (*vBegin[getSmallestValid()]).first;
        }

        bool operator!=( const MergeIterator& rOther ) const
        {
            for(size_t uiI = 0; uiI < vBegin.size(); uiI++)
                if(vBegin[uiI] != rOther.vBegin[uiI])
                {
                    if(**this >= uiEnd)
                        return false;
                    return true;
                }
            return false;
        }

        friend std::ostream& operator<<( std::ostream& os, const MergeIterator& rIt )
        {
            return os << rIt.vBegin;
        }
    };
    class CordIterator
    {
        cord_it_t xA;

      public:
        CordIterator( cord_it_t xA )
            : xA( xA )
        {}

        void operator++( )
        {
            ++xA;
        }

        const coordinate_t operator*( ) const
        {
            return ( *xA ).first;
        }

        bool operator!=( const CordIterator& rOther ) const
        {
            return xA != rOther.xA;
        }

        friend std::ostream& operator<<( std::ostream& os, const CordIterator& rIt )
        {
            os << rIt.xA;

            return os;
        }
    };


    class PointIterator
    {
        point_it_t xIt;
        size_t uiDim;

      public:
        PointIterator( point_it_t xIt, size_t uiDim )
            : xIt( xIt ),
              uiDim( uiDim )
        {}

        void operator++( )
        {
            coordinate_t uiLast = **this;
            while( !xIt.eof( ) && uiLast == **this )
                ++xIt;
        }

        const coordinate_t operator*( ) const
        {
            return xIt->vPos[ uiDim ];
        }

        bool operator!=( const PointIterator& rOther ) const
        {
            return xIt != rOther.xIt;
        }

        friend std::ostream& operator<<( std::ostream& os, const PointIterator& rIt )
        {
            return os;
        }
    };


  public:
    std::array<red_entry_arr_t, D> vSparseCoordsOverlay;
    entry_arr_t vSparseCoordsInternal;
    std::array<typename prefix_sum_grid_t::template Entry<D - 1>, D> vOverlayEntries;
    typename prefix_sum_grid_t::template Entry<D> xInternalEntires;
    typename points_t::Entry xPoints;

    Overlay( )
        : vSparseCoordsOverlay{ }, vSparseCoordsInternal{ }, vOverlayEntries{ }, xInternalEntires{ }, xPoints{}
    {}

    void generate( const overlay_grid_t& rOverlays,
                   sparse_coord_t& rSparseCoords, prefix_sum_grid_t& rPrefixSums, points_t& vPoints,
                   typename points_t::Entry xPoints, std::array<std::vector<coordinate_t>, D> vPredecessors,
                   pos_t vMyBottomLeft, pos_t vPosTopRight, dataset_t* pDataset, 
                   size_t uiOverlaysNow, size_t uiOverlaysTotal, progress_stream_t& xProg )
    {
        this->xPoints = xPoints;
        // construct sparse coordinates for each dimension
        xProg << Verbosity(1) << "constructing sparse coordinates for overlay bottom left= " << vMyBottomLeft
              << " top right " << vPosTopRight << "\n";
        for( size_t uiI = 0; uiI < D; uiI++ )
        {
            xProg << Verbosity(2) << "dim " << uiI << "\n";
            for( size_t uiJ = 0; uiJ < D-1; uiJ++ )
            {
                size_t uiJAct = uiJ + (uiJ >= uiI ? 1 : 0);
                xProg << Verbosity(2) << "sub dim " << uiJAct << " (" << uiJ << ")" << "\n";
                std::vector<cord_it_t> vBegin {};
                std::vector<cord_it_t> vEnd {};

                // add coordinates from previous overlays to the overlay entries
                for( size_t uiI2 = 0; uiI2 < D; uiI2++ )
                    if( uiJAct != uiI2 )
                        for(coordinate_t uiPred : vPredecessors[ uiI2 ])
                        {
                            const Overlay* pPred = &rOverlays.vData[uiPred];
                            vBegin.push_back( rSparseCoords.cbegin( 
                                pPred->vSparseCoordsOverlay[ uiI ][ uiJ ]) );
                            vEnd.push_back( rSparseCoords.cend( 
                                pPred->vSparseCoordsOverlay[ uiI ][ uiJ ] ) );
                            if(xProg.active())
                                pPred->vSparseCoordsOverlay[ uiI ][ uiJ ].stream(
                                        std::cout << "from predecessor in dim " << uiI2 
                                                    << " overlay: ", rSparseCoords) << std::endl;
                            if constexpr(DEPENDANT_DIMENSION)
                                if(uiI == 1 && uiI2 != 1 /* <- cause this is done below anyways */)
                                    // predecessors in dim 1 could reach below the start of this overlay ->
                                    // in that case their overlay coords are not sufficient & we also have to 
                                    // use their internal coords (@todo actually we only need some of their 
                                    // internal coords -> could save mem here)
                                {
                                    // @todo pick out the correct points and no more
                                    vBegin.push_back(rSparseCoords.cbegin( pPred->vSparseCoordsInternal[ uiJAct ] ));
                                    vEnd.push_back(rSparseCoords.cend( pPred->vSparseCoordsInternal[ uiJAct ] ));
                                    if(xProg.active())
                                        pPred->vSparseCoordsInternal[ uiJAct ].stream(
                                                std::cout << "from internal dim 1: ", rSparseCoords) << std::endl;
                                }
                        }


                // add coordinates from the points of the previous overlay to the overlay entries
                for(coordinate_t uiPred : vPredecessors[ uiI ])
                {
                    const Overlay* pPred = &rOverlays.vData[uiPred];
                    vBegin.push_back(rSparseCoords.cbegin( pPred->vSparseCoordsInternal[ uiJAct ] ));
                    vEnd.push_back(rSparseCoords.cend( pPred->vSparseCoordsInternal[ uiJAct ] ));
                    if(xProg.active())
                        pPred->vSparseCoordsInternal[ uiJAct ].stream(
                                std::cout << "from internal: ", rSparseCoords) << std::endl;
                }

                xProg << "from bottom left: { " << vMyBottomLeft[ uiJAct ] - 1 << " }\n";

                MergeIterator xBegin( vBegin, vEnd, vPosTopRight[uiJAct] );
                MergeIterator xEnd( vEnd, vEnd, vPosTopRight[uiJAct] );

                // skip over positions that are before and after this overlay
                while(xBegin != xEnd && *xBegin < vMyBottomLeft[uiJAct])
                    ++xBegin;

                if(vPredecessors[ uiJAct ].size() > 0) 
                    vSparseCoordsOverlay[ uiI ][ uiJ ] = rSparseCoords.addStart( xBegin, xEnd, 
                                                                                 vMyBottomLeft[ uiJAct ] - 1 );
                else
                    vSparseCoordsOverlay[ uiI ][ uiJ ] = rSparseCoords.add( xBegin, xEnd );

                if(xProg.active())
                    vSparseCoordsOverlay[ uiI ][ uiJ ].stream(std::cout << "result: ", rSparseCoords) << std::endl;
            }
        }

        if(xPoints.size() > 0)
        {
            xProg << Verbosity(1) << "constructing sparse coordinates for points\n";
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
                xProg << Verbosity(2) << "dim " << uiI << "\n";
                vPoints.sortByDim( uiI, xPoints );

                if(xProg.active())
                    xPoints.stream(std::cout << "from points: ", vPoints) << std::endl;

                vSparseCoordsInternal[ uiI ] = rSparseCoords.add(
                    PointIterator( vPoints.cbegin( xPoints ), uiI ),
                    PointIterator( vPoints.cend( xPoints ), uiI ) );
                if(xProg.active())
                    vSparseCoordsInternal[uiI].stream(std::cout << "result: ", rSparseCoords) << std::endl;
            }

            // construct internal grid
            xProg << Verbosity(1) << "constructing internal grid\n";
            pos_t vInternalAxisSizes = rSparseCoords.axisSizes( vSparseCoordsInternal );
            xProg << Verbosity(2) << "axis sizes: " << vInternalAxisSizes << "\n";
            xInternalEntires = rPrefixSums.add( vInternalAxisSizes );
            vPoints.iterate(
                [ & ]( const point_t& xPoint ) {
                    rPrefixSums.get( rSparseCoords.sparse( xPoint.vPos, vSparseCoordsInternal ),
                                    xInternalEntires ) += 1;
                },
                xPoints );

            xProg << "vSparseCoordsOverlay " << vSparseCoordsOverlay << "\n";
            xProg << "vSparseCoordsInternal " << vSparseCoordsInternal << "\n";
            xProg << "rSparseCoords " << rSparseCoords << "\n";

            // compute internal prefix sum
            coordinate_t uiNumTotal = rPrefixSums.sizeOf( xInternalEntires ) * D;
            coordinate_t uiNumDone = 0;
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
                xProg << Verbosity(3) << "computing prefix sums over dimension " << uiI << "\n";
                red_entry_arr_t vRelevantSparseCoordsInternal = relevant( vSparseCoordsInternal, uiI );
                rSparseCoords.template iterate<D - 1>(
                    [ & ]( const red_pos_t&, const red_pos_t& vTo ) {
                        pos_t vFullTo = expand( vTo, uiI );
                        val_t uiPrefixSum = 0;
                        xProg << Verbosity(3) << "starting...: " << vFullTo << ": " << uiPrefixSum << "\n";
                        rSparseCoords.iterate(
                            [ & ]( coordinate_t, coordinate_t uiTo ) {
                                vFullTo[ uiI ] = uiTo;
                                uiPrefixSum += rPrefixSums.get( vFullTo, xInternalEntires );
                                xProg << Verbosity( 3 ) << vFullTo << ": " << uiPrefixSum << "\n";
                                rPrefixSums.get( vFullTo, xInternalEntires ) = uiPrefixSum;
                                
                                ++uiNumDone;
                                if( xProg.printAgain( ) )
                                    xProg << Verbosity( 0 ) 
                                          << uiOverlaysNow << " out of " << uiOverlaysTotal << " overlays, thats "
                                          << 100.0 * ((double)uiOverlaysNow / (double)uiOverlaysTotal) << "%. " 
                                          << uiNumDone << " out of " << uiNumTotal << " prefix sums, thats "
                                          << 100.0 * ((double)uiNumDone / (double)uiNumTotal) << "%.\n";
                            },
                            vSparseCoordsInternal[ uiI ] );
                    },
                    vRelevantSparseCoordsInternal );
            }
        }

        // construct overlay sum grid
        xProg << Verbosity(1) << "constructing overlay sum grid\n";
        for( size_t uiI = 0; uiI < D; uiI++ )
            if( vPredecessors[ uiI ].size() > 0 )
            {
                xProg << Verbosity(2) << "dim " << uiI << "\n";

                red_pos_t vAxisSizes = rSparseCoords.axisSizes( vSparseCoordsOverlay[uiI] );

                vOverlayEntries[ uiI ] = rPrefixSums.add( vAxisSizes );
                assert(vMyBottomLeft[uiI] > 0);
                if(rPrefixSums.sizeOf(vOverlayEntries[ uiI ]) > 0)
                    rSparseCoords.template iterate<D - 1>(
                        [ & ]( const red_pos_t& vFrom, const red_pos_t& vTo ) {
                            pos_t vFullFrom = expand( vFrom, uiI );
                            vFullFrom[ uiI ] = vMyBottomLeft[ uiI ] - 1;

                            xProg << Verbosity(3) << "query " << vFullFrom << "\n";
                            auto uiRet = pDataset->get( rOverlays, rSparseCoords, rPrefixSums, vFullFrom, xProg );
                            //auto uiRet = vPredecessors[ uiI ]->get( rSparseCoords, rPrefixSums, vFullFrom, xProg );
                            xProg << Verbosity(3) << "query " << vFullFrom << ": " << uiRet << "\n";

                            rPrefixSums.get( vTo, vOverlayEntries[ uiI ] ) = uiRet;
                        },
                        vSparseCoordsOverlay[uiI] );
            }
        xProg << Verbosity(1) << "done\n";
    }


    val_t get( const sparse_coord_t& rSparseCoords, const prefix_sum_grid_t& rPrefixSums, pos_t vCoords,
               pos_t vMyBottomLeft, progress_stream_t& xProg ) const
    {
        val_t uiRet = 0;

        xProg << Verbosity(2) << "\tquerying overlay for " << vCoords << "...\n";

        for( size_t uiI = 0; uiI < D; uiI++ )
            --vMyBottomLeft[uiI]; // will turn zero values into max values
        xProg << Verbosity(3) << "\t\tvMyBottomLeft " << vMyBottomLeft << "\n";

        forAllCombinations<pos_t>(
            [&](pos_t vPos, size_t uiDistToTo)
            {
                if(uiDistToTo == 0)
                    return;
                size_t uiI = 0;
                while(uiI < D && (vPos[uiI] != vMyBottomLeft[uiI] || 
                        vSparseCoordsOverlay[ uiI ][0].uiStartIndex == std::numeric_limits<coordinate_t>::max() ))
                    ++uiI;

                if(uiI == D)
                    return;

                xProg << Verbosity(3) << "\t\tquery: " << vPos << " in overlay " << uiI << "\n";

                red_pos_t vRelevant = relevant( vPos, uiI );
                red_pos_t vSparse = rSparseCoords.sparse(vRelevant, vSparseCoordsOverlay[ uiI ] );

                xProg << "\t\trelevant: " << vRelevant << " sparse: " << vSparse << "\n";

                val_t uiCurr = rPrefixSums.get( vSparse, vOverlayEntries[ uiI ] );

                xProg << "\t\tis " << (uiDistToTo % 2 == 0 ? "-" : "+") << uiCurr << "\n";
                uiRet += uiCurr * (uiDistToTo % 2 == 0 ? -1 : 1);
            },
            vMyBottomLeft,
            vCoords,
            [](coordinate_t uiPos){return uiPos != std::numeric_limits<coordinate_t>::max();}
        );

        if( rPrefixSums.sizeOf( xInternalEntires ) > 0 )
        {
            auto vSparseCoords = rSparseCoords.sparse( vCoords, vSparseCoordsInternal );
            auto uiCurr = rPrefixSums.get( vSparseCoords, xInternalEntires );
            xProg << Verbosity(2) << "\tquerying internal " << vCoords << " -> " << vSparseCoords << ": +" 
                  << uiCurr << "\n";
            uiRet += uiCurr;
        }
        else
            xProg << Verbosity(2) << "\tno internal entries\n";
        

        return uiRet;
    }

    template <typename T> std::array<T, D - 1> relevant( const std::array<T, D>& vAllEntries, size_t uiI ) const
    {
        std::array<T, D - 1> vRelevantEntries;
        for( size_t uiJ = 0; uiJ < uiI; uiJ++ )
            vRelevantEntries[ uiJ ] = vAllEntries[ uiJ ];
        for( size_t uiJ = uiI + 1; uiJ < D; uiJ++ )
            vRelevantEntries[ uiJ - 1 ] = vAllEntries[ uiJ ];
        return vRelevantEntries;
    }

    template <typename T> std::array<T, D> expand( const std::array<T, D - 1>& vCompressed, size_t uiI ) const
    {
        std::array<T, D> vAllEntries;
        for( size_t uiJ = 0; uiJ < uiI; uiJ++ )
            vAllEntries[ uiJ ] = vCompressed[ uiJ ];
        for( size_t uiJ = uiI + 1; uiJ < D; uiJ++ )
            vAllEntries[ uiJ ] = vCompressed[ uiJ - 1 ];
        return vAllEntries;
    }

    friend std::ostream& operator<<( std::ostream& os, const Overlay& rOverlay )
    {
        os << "<" << std::endl;
        os << "\tvSparseCoordsOverlay: ";
        os << rOverlay.vSparseCoordsOverlay << std::endl;

        os << "\tvSparseCoordsInternal: ";
        os << rOverlay.vSparseCoordsInternal << std::endl;

        os << "\tvOverlayEntries: ";
        os << rOverlay.vOverlayEntries << std::endl;

        os << "\txInternalEntires: ";
        os << rOverlay.xInternalEntires << std::endl;

        os << "\txPoints: ";
        os << rOverlay.xPoints << std::endl;

        os << ">";

        return os;
    }

    std::ostream& stream( std::ostream& os, pos_t vGridPos, const sparse_coord_t& rSparseCoords, 
                          const prefix_sum_grid_t& rPrefixSums, const dataset_t& rData ) const
    {
        os << "<";
        if(rData.exists(rSparseCoords, vGridPos))
        {
            os << std::endl;
            auto vPosActual = rData.actualFromGridPos(rSparseCoords, vGridPos);
            os << "\tbottom left: " << vPosActual << "\n";
            auto vPosTR = rData.actualTopRightFromGridPos(rSparseCoords, vGridPos);
            os << "\ttop right: " << vPosTR << "\n";
            os << "\tvSparseCoordsOverlay: ";
            for(size_t uiI = 0; uiI < D; uiI++)
            {
                os << uiI << ":<";
                for(size_t uiJ = 0; uiJ < D-1; uiJ++)
                {
                    os << " " << (uiJ + (uiJ >= uiI ? 1 : 0)) << ":";
                    vSparseCoordsOverlay[uiI][uiJ].stream(os, rSparseCoords) << " ";
                }
                os << "> ";
            }
            os << std::endl;

            os << "\tvSparseCoordsInternal: ";
            for(size_t uiI = 0; uiI < D; uiI++)
                vSparseCoordsInternal[uiI].stream(os, rSparseCoords) << " ";
            os << std::endl;

            os << "\tvOverlayEntries: ";
            for(size_t uiI = 0; uiI < D; uiI++)
                vOverlayEntries[uiI].streamOp(os, rPrefixSums) << " ";
            os << std::endl;

            os << "\txInternalEntires: ";
            xInternalEntires.streamOp(os, rPrefixSums) << std::endl;
        }
        else
            os << "n/a";

        return os;
    }
    std::ostream& stream( std::ostream& os, pos_t vGridPos, const sparse_coord_t& rSparseCoords, 
                          const prefix_sum_grid_t& rPrefixSums, const dataset_t& rData, const points_t& vPoints ) const
    {
        this->stream(os, vGridPos, rSparseCoords, rPrefixSums, rData);

        if(rData.exists(rSparseCoords, vGridPos))
        {
            os << "\txPoints: ";
            xPoints.stream(os, vPoints) << std::endl;
        }

        os << ">";

        return os;
    }

    std::ostream& stream( std::ostream& os, pos_t vGridPos, const sparse_coord_t& rSparseCoords, 
                          const prefix_sum_grid_t& rPrefixSums, const dataset_t& rData, const points_t& vPoints, const desc_t& vDesc ) const
    {
        this->stream(os, vGridPos, rSparseCoords, rPrefixSums, rData);

        if(rData.exists(rSparseCoords, vGridPos))
        {
            os << "\txPoints: ";
            xPoints.stream(os, vPoints, vDesc) << std::endl;
        }

        os << ">";

        return os;
    }

};


} // namespace sps
