#pragma once

#include "sps/nd_grid.h"
#include "sps/points.h"
#include "sps/sparse_coordinate.h"
#include <cassert>
#include <functional>
#include <string>



namespace sps
{

template <typename type_defs> class Overlay
{
    EXTRACT_TYPE_DEFS; // macro call

    using sparse_coord_t = SparseCoord<type_defs>;
    using prefix_sum_grid_t = NDGrid<type_defs, val_t>;


    using point_t = Point<type_defs>;
    using points_t = Points<type_defs>;
    using desc_t = Desc<type_defs>;

    using cord_it_t = typename sparse_coord_t::EntryIterator;
    using point_it_t = typename points_t::EntryIterator;

    using red_pos_t = std::array<coordinate_t, D - 1>;

    using entry_arr_t = std::array<typename sparse_coord_t::Entry, D>;
    using red_entry_arr_t = std::array<typename sparse_coord_t::Entry, D - 1>;

    class MergeIterator
    {
        cord_it_t xA, xB;
        cord_it_t xAE, xBE;

      public:
        MergeIterator( cord_it_t xA, cord_it_t xB, cord_it_t xAE, cord_it_t xBE )
            : xA( xA ), xB( xB ), xAE( xAE ), xBE( xBE )
        {}

        void operator++( )
        {
            bool bAValid = xA != xAE;
            bool bBValid = xB != xBE;
            if( bAValid && bBValid && ( *xA ).first == ( *xB ).first )
            {
                ++xA;
                ++xB;
            }
            else if( bAValid && ( !bBValid || ( *xA ).first < ( *xB ).first ) )
                ++xA;
            else if( bBValid && ( !bAValid || ( *xB ).first < ( *xA ).first ) )
                ++xB;
            else
                assert( false );
        }

        const coordinate_t operator*( ) const
        {
            bool bAValid = xA != xAE;
            bool bBValid = xB != xBE;
            assert( bAValid || bBValid );
            if( bAValid && ( !bBValid || ( *xA ).first < ( *xB ).first ) )
                return ( *xA ).first;
            return ( *xB ).first;
        }

        bool operator!=( const MergeIterator& rOther ) const
        {
            return xA != rOther.xA || xB != rOther.xB;
        }

        friend std::ostream& operator<<( std::ostream& os, const MergeIterator& rIt )
        {
            os << rIt.xA;

            os << " & ";
            os << rIt.xB;

            return os;
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
        PointIterator( const sparse_coord_t& rSparseCoords, const entry_arr_t& rGlobalSparseCords,
                       const Overlay& rOverlay, point_it_t xIt, size_t uiDim )
            : rSparseCoords( rSparseCoords ),
              rGlobalSparseCords( rGlobalSparseCords ),
              rOverlay( rOverlay ),
              xIt( xIt ),
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
    val_t uiBottomLeft;

    Overlay( )
        : vSparseCoordsOverlay{ }, vSparseCoordsInternal{ }, vOverlayEntries{ }, xInternalEntires{ }, uiBottomLeft{ }
    {}

    void generate( sparse_coord_t& rSparseCoords, prefix_sum_grid_t& rPrefixSums, points_t& vPoints,
                   typename points_t::Entry xPoints, std::array<Overlay*, D> vPredecessors,
                   std::array<coordinate_t, D> vPos, val_t uiBottomLeft, const entry_arr_t& rGlobalSparseCords,
                   std::optional<progress_stream_t> xProg )
    {
        this->xPoints = xPoints;
        this->uiBottomLeft = uiBottomLeft;
        // construct sparse coordinates for each dimension
        xProg << "constructing sparse coordinates for overlay\n";
        for( size_t uiI = 0; uiI < D; uiI++ )
            if( vPredecessors[ uiI ] != nullptr )
            {
                xProg << "dim " << uiI << "\n";
                for( size_t uiJ = 0; uiJ < D-1; uiJ++ )
                {
                    // add coordinates from previous overlay to the overlay entries
                    auto xA = rSparseCoords.cbegin( vPredecessors[ uiI ]->vSparseCoordsOverlay[ uiI ][ uiJ ] );
                    auto xAE = rSparseCoords.cend( vPredecessors[ uiI ]->vSparseCoordsOverlay[ uiI ][ uiJ ] );

                    auto xB = rSparseCoords.cbegin( 
                        vPredecessors[ uiI ]->vSparseCoordsInternal[ uiJ + (uiJ >= uiI ? 1 : 0) ] );
                    auto xBE = rSparseCoords.cend( 
                        vPredecessors[ uiI ]->vSparseCoordsInternal[ uiJ + (uiJ >= uiI ? 1 : 0) ] );

                    if(xProg)
                    {
                        vPredecessors[ uiI ]->vSparseCoordsOverlay[ uiI ][ uiJ ].stream(
                                std::cout << "from overlay: ", rSparseCoords) << std::endl;
                        vPredecessors[ uiI ]->vSparseCoordsInternal[ uiJ + (uiJ >= uiI ? 1 : 0) ].stream(
                                std::cout << "from internal: ", rSparseCoords) << std::endl;
                    }

                    MergeIterator xBegin( xA, xB, xAE, xBE );
                    MergeIterator xEnd( xAE, xBE, xAE, xBE );

                    vSparseCoordsOverlay[ uiI ][ uiJ ] = rSparseCoords.add( xBegin, xEnd );

                    if(xProg)
                        vSparseCoordsOverlay[ uiI ][ uiJ ].stream(std::cout << "result: ", rSparseCoords) << std::endl;
                }
            }

        if(xPoints.size() > 0)
        {
            xProg << "constructing sparse coordinates for points\n";
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
                xProg << "dim " << uiI << "\n";
                // add coordinates from previous overlay to the overlay entries
                vPoints.sortByDim( uiI, xPoints );

                if(xProg)
                    xPoints.stream(std::cout << "from points: ", vPoints) << std::endl;

                vSparseCoordsInternal[ uiI ] = rSparseCoords.add(
                    PointIterator( rSparseCoords, rGlobalSparseCords, *this, vPoints.cbegin( xPoints ), uiI ),
                    PointIterator( rSparseCoords, rGlobalSparseCords, *this, vPoints.cend( xPoints ), uiI ) );
                if(xProg)
                    vSparseCoordsInternal[uiI].stream(std::cout << "result: ", rSparseCoords) << std::endl;
            }

            // construct internal grid
            xProg << "constructing internal grid\n";
            pos_t vInternalAxisSizes = rSparseCoords.axisSizes( vSparseCoordsInternal );
            xProg << "axis sizes: " << vInternalAxisSizes << "\n";
            xInternalEntires = rPrefixSums.add( vInternalAxisSizes );
            //typename prefix_sum_grid_t::template Entry<D> xInternalEntiresHelper = 
            //    rPrefixSums.add( vInternalAxisSizes );
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
            for( size_t uiI = 0; uiI < D; uiI++ )
            {
                xProg << "integrating over dimension " << uiI << "\n";
                red_entry_arr_t vRelevantSparseCoordsInternal = relevant( vSparseCoordsInternal, uiI );
                rSparseCoords.template iterate<D - 1>(
                    [ & ]( const red_pos_t&, const red_pos_t& vTo ) {
                        pos_t vFullTo = expand( vTo, uiI );
                        val_t uiPrefixSum = 0;
                        xProg << "starting...: " << vFullTo << ": " << uiPrefixSum << "\n";
                        rSparseCoords.iterate(
                            [ & ]( coordinate_t, coordinate_t uiTo ) {
                                vFullTo[ uiI ] = uiTo;
                                uiPrefixSum += rPrefixSums.get( vFullTo, xInternalEntires );
                                xProg << vFullTo << ": " << uiPrefixSum << "\n";
                                rPrefixSums.get( vFullTo, xInternalEntires ) = uiPrefixSum;
                            },
                            vSparseCoordsInternal[ uiI ] );
                    },
                    vRelevantSparseCoordsInternal );
            }
            //xProg << "clearing tmp storage\n";
            //rPrefixSums.remove( xInternalEntiresHelper );
        }

        // construct overlay sum grid
        xProg << "constructing overlay sum grid\n";
        for( size_t uiI = 0; uiI < D; uiI++ )
            if( vPredecessors[ uiI ] != nullptr )
            {
                xProg << "dim " << uiI << "\n";

                red_pos_t vAxisSizes = rSparseCoords.axisSizes( vSparseCoordsOverlay[uiI] );
                vOverlayEntries[ uiI ] = rPrefixSums.add( vAxisSizes );

                rSparseCoords.template iterate<D - 1>(
                    [ & ]( const red_pos_t& vFrom, const red_pos_t& vTo ) {
                        pos_t vFullFrom = expand( vFrom, uiI );
                        for( size_t uiI = 0; uiI < D; uiI++ )
                            vFullFrom[ uiI ]++;
                        vFullFrom[ uiI ] = vPos[ uiI ];

                        rPrefixSums.get( vTo, vOverlayEntries[ uiI ] ) =
                            vPredecessors[ uiI ]->get( rSparseCoords, rPrefixSums, vFullFrom, xProg );
                    },
                    vSparseCoordsOverlay[uiI] );
            }
        xProg << "done\n";
    }

    val_t get( const sparse_coord_t& rSparseCoords, const prefix_sum_grid_t& rPrefixSums, pos_t vCoords,
               std::optional<progress_stream_t> xProg = { } ) const
    {
        val_t uiRet = 0;

        xProg << "querying overlay...\n";
        for( size_t uiI = 0; uiI < D; uiI++ )
            if( rPrefixSums.sizeOf( vOverlayEntries[ uiI ] ) > 0 )
            {
                auto vRelevant = relevant( vCoords, uiI );
                auto vSparse = rSparseCoords.sparse(vRelevant , vSparseCoordsOverlay[ uiI ] );
                uiRet += rPrefixSums.get( vSparse, vOverlayEntries[ uiI ] );
                xProg << "dim " << uiI << " " << vRelevant << " -> " << vSparse << ": " << uiRet << "\n";
            }

        assert( uiRet >= ( D - 1 ) * uiBottomLeft );
        xProg << "bottom left " << uiBottomLeft << "\n";
        uiRet -= ( D - 1 ) * uiBottomLeft;

        if( rPrefixSums.sizeOf( xInternalEntires ) > 0 )
        {
            auto vSparseCoords = rSparseCoords.sparse( vCoords, vSparseCoordsInternal );
            uiRet += rPrefixSums.get( vSparseCoords, xInternalEntires );
            xProg << "querying internal " << vCoords << " -> " << vSparseCoords << ": " << uiRet << "\n";
        }

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

        os << "\tuiBottomLeft: ";
        os << rOverlay.uiBottomLeft << std::endl;
        os << ">";

        return os;
    }

    std::ostream& stream( std::ostream& os, const sparse_coord_t& rSparseCoords, 
                          const prefix_sum_grid_t& rPrefixSums, const points_t& vPoints, const desc_t& vDesc ) const
    {
        os << "<" << std::endl;
        os << "\tvSparseCoordsOverlay: ";
        for(size_t uiI = 0; uiI < D; uiI++)
            for(size_t uiJ = 0; uiJ < D-1; uiJ++)
                vSparseCoordsOverlay[uiI][uiJ].stream(os, rSparseCoords) << " ";
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

        os << "\txPoints: ";
        xPoints.stream(os, vPoints, vDesc) << std::endl;

        os << "\tuiBottomLeft: ";
        os << uiBottomLeft << std::endl;
        os << ">";

        return os;
    }

    std::ostream& stream( std::ostream& os, const sparse_coord_t& rSparseCoords, 
                          const prefix_sum_grid_t& rPrefixSums, const points_t& vPoints ) const
    {
        os << "<" << std::endl;
        os << "\tvSparseCoordsOverlay: ";
        for(size_t uiI = 0; uiI < D; uiI++)
            for(size_t uiJ = 0; uiJ < D-1; uiJ++)
                vSparseCoordsOverlay[uiI][uiJ].stream(os, rSparseCoords) << " ";
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

        os << "\txPoints: ";
        xPoints.stream(os, vPoints) << std::endl;

        os << "\tuiBottomLeft: ";
        os << uiBottomLeft << std::endl;
        os << ">";

        return os;
    }
};


} // namespace sps
