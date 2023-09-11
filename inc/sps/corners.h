#pragma once

#include "sps/desc.h"
#include "sps/point.h"
#include "sps/type_defs.h"
#include "sps/util.h"

#include <cassert>
#include <functional>
#include <string>


namespace sps
{

template <typename type_defs> class Corners
{
    EXTRACT_TYPE_DEFS; // macro call

    using corner_t = AlignedPower2<Corner<type_defs>>;

#if 0
    template<int s> struct CheckSizeOfPoint;
    CheckSizeOfPoint<sizeof(Point<type_defs>)> xCheckSizeOfPoint;
    CheckSizeOfPoint<sizeof(corner_t)> xCheckSizeOfAlignedPoint;
#endif

  public:
    EXTRACT_VEC_GENERATOR( points, corner_t ); // macro call
    points_file_t xFile;
    points_vec_t vData;

    static constexpr bool THREADSAVE = points_THREADSAVE;

  private:
    struct PointsComperator
    {
        const size_t uiDim;

        PointsComperator( size_t uiDim ) : uiDim( uiDim )
        {}

        bool operator( )( const corner_t& a, const corner_t& b ) const
        {
            return a.vPos[ uiDim ] < b.vPos[ uiDim ];
        }

        corner_t min_value( ) const
        {
            corner_t xRet{ };
            return xRet;
        };

        corner_t max_value( ) const
        {
            corner_t xRet{ };
            for( size_t uiI = 0; uiI < D; uiI++ )
                xRet.vPos[ uiI ] = std::numeric_limits<coordinate_t>::max( );
            return xRet;
        };
    };
    struct PointsComperator2
    {
        const size_t uiDim1, uiDim2;

        PointsComperator2( size_t uiDim1, size_t uiDim2 ) : uiDim1( uiDim1 ), uiDim2( uiDim2 )
        {}

        bool operator( )( const corner_t& a, const corner_t& b ) const
        {
            if( a.vPos[ uiDim1 ] < b.vPos[ uiDim1 ] )
                return true;
            if( a.vPos[ uiDim1 ] > b.vPos[ uiDim1 ] )
                return false;
            return a.vPos[ uiDim2 ] < b.vPos[ uiDim2 ];
        }

        corner_t min_value( ) const
        {
            corner_t xRet{ };
            return xRet;
        };

        corner_t max_value( ) const
        {
            corner_t xRet{ };
            for( size_t uiI = 0; uiI < D; uiI++ )
                xRet.vPos[ uiI ] = std::numeric_limits<coordinate_t>::max( );
            return xRet;
        };
    };

  public:
    struct Entry
    {
        coordinate_t uiStartIndex;
        coordinate_t uiEndIndex;

        friend std::ostream& operator<<( std::ostream& os, const Entry& rEntry )
        {
            std::cout << "s";
            std::cout << rEntry.uiStartIndex;

            std::cout << " e";
            std::cout << rEntry.uiEndIndex;

            return os;
        }

        coordinate_t size( ) const
        {
            return uiEndIndex - uiStartIndex;
        }

        std::ostream& stream( std::ostream& os, const Corners& rCorners ) const
        {
            os << "{ ";
            for( size_t uiI = uiStartIndex; uiI < uiEndIndex; uiI++ )
            {
                if( uiI > uiStartIndex )
                    os << ", ";
                rCorners.vData[ uiI ].stream( os );
            }
            os << " }";

            return os;
        }
    };

    class EntryIterator
    {
        const Corners& rCorners;
        const Entry& rInfo;
        size_t uiI;

      public:
        EntryIterator( const Corners& rCorners, const Entry& rInfo )
            : rCorners( rCorners ), rInfo( rInfo ), uiI( rInfo.uiStartIndex )
        {}

        void operator++( )
        {
            assert( !eof( ) );
            uiI++;
        }

        const corner_t& operator*( ) const
        {
            return rCorners.vData[ uiI ];
        }

        const corner_t* operator->( ) const
        {
            return &rCorners.vData[ uiI ];
        }

        bool operator!=( const EntryIterator& rOther ) const
        {
            return uiI != rOther.uiI;
        }

        bool eof( ) const
        {
            return uiI == rInfo.uiEndIndex;
        }

        friend std::ostream& operator<<( std::ostream& os, const EntryIterator& rIt )
        {
            os << rIt.rInfo;

            os << " ";
            os << rIt.uiI;

            return os;
        }

        friend class Corners;
    };

    Corners( std::string sPrefix = "", bool bWrite = false )
        : xFile( points_vec_generator.file( sPrefix + ".corners", bWrite ) ), vData( points_vec_generator.vec( xFile ) )
    {}

    template <bool trigger = !IS_ORTHOTOPE> typename std::enable_if_t<trigger> add( pos_t vPos, val_t uiVal )
    {
        vData.push_back( corner_t( vPos, uiVal ) );
    }

    template <bool trigger = IS_ORTHOTOPE>
    typename std::enable_if_t<trigger> add( pos_t vStart, pos_t vEnd, val_t uiVal )
    {
        forAllCombinationsN<pos_t, ORTHOTOPE_DIMS>(
            [ & ]( size_t uiI, size_t, size_t, pos_t vPos ) {
                for( size_t uiD = ORTHOTOPE_DIMS; uiD < D; uiD++ )
                {
                    assert( vStart[ uiD ] == vEnd[ uiD ] );
                    vPos[ uiD ] = vStart[ uiD ];
                }
                vData.push_back( corner_t( vPos, uiVal, uiI ) );
            },
            vStart, vEnd );
    }

    EntryIterator cbegin( const Entry& rInfo ) const
    {
        return EntryIterator( *this, rInfo );
    }

    EntryIterator cend( const Entry& rInfo ) const
    {
        EntryIterator xRet( *this, rInfo );
        xRet.uiI = rInfo.uiEndIndex;
        return xRet;
    }

    void iterate( std::function<void( const corner_t& )> fDo, const Entry& rEntry ) const
    {
        auto itEnd = vData.cbegin( ) + rEntry.uiEndIndex;
        for( auto cIter = vData.cbegin( ) + rEntry.uiStartIndex; cIter != itEnd; cIter++ )
            fDo( *cIter );
    }

    void iterate( std::function<void( corner_t& )> fDo, const Entry& rEntry )
    {
        auto itEnd = vData.begin( ) + rEntry.uiEndIndex;
        for( auto cIter = vData.begin( ) + rEntry.uiStartIndex; cIter != itEnd; cIter++ )
            fDo( *cIter );
    }

    void forEqualRange( std::function<bool( const corner_t& )> fBefore, std::function<bool( const corner_t& )> fAfter,
                        std::function<void( const corner_t& )> fDo, const Entry& rEntry ) const
    {
        size_t uiStartIdx = rEntry.uiStartIndex;
        size_t uiEndIdx = rEntry.uiEndIndex;
        while( uiStartIdx < uiEndIdx )
        {
            size_t uiMid = ( uiStartIdx + uiEndIdx ) / 2;
            if( fBefore( vData[ uiMid ] ) )
                uiStartIdx = uiMid + 1;
            else
                uiEndIdx = uiMid;
        }

        while( uiStartIdx < rEntry.uiEndIndex && fBefore( vData[ uiStartIdx ] ) )
            ++uiStartIdx;

        while( uiStartIdx < rEntry.uiEndIndex && !fAfter( vData[ uiStartIdx ] ) )
        {
            fDo( vData[ uiStartIdx ] );
            ++uiStartIdx;
        }
    }

    bool popEntry( const Entry& rEntry )
    {
        if( rEntry.uiEndIndex != vData.size( ) )
            return false;

        vData.resize( rEntry.uiStartIndex );
        return true;
    }

    Entry copyEntry( const Entry& rEntry )
    {
        Entry xRet;
        xRet.uiStartIndex = vData.size( );
        // make sure no new space needs to be allocated during copy operation (this allocation could lead to segfaults)
        vData.reserve( vData.size( ) + rEntry.size( ) );
        iterate( [ & ]( const corner_t& rP ) { vData.push_back( rP ); }, rEntry );
        xRet.uiEndIndex = vData.size( );
        return xRet;
    }

    void sortByDim( size_t uiDim, const Entry& rEntry )
    {
        points_sort_func_t<typename points_vec_t::iterator, PointsComperator>( )(
            vData.begin( ) + rEntry.uiStartIndex, vData.begin( ) + rEntry.uiEndIndex, PointsComperator( uiDim ) );
    }

    void sortByDim( size_t uiDim1, size_t uiDim2, const Entry& rEntry )
    {
        points_sort_func_t<typename points_vec_t::iterator, PointsComperator2>( )(
            vData.begin( ) + rEntry.uiStartIndex,
            vData.begin( ) + rEntry.uiEndIndex,
            PointsComperator2( uiDim1, uiDim2 ) );
    }

    size_t size( ) const
    {
        return vData.size( );
    }

    void clear( )
    {
        vData.clear( );
    }

    friend std::ostream& operator<<( std::ostream& os, const Corners& vCorners )
    {
        size_t uiX = 0;
        for( const auto& rP : vCorners.vData )
            os << uiX++ << ": " << rP << std::endl;
        return os;
    }
};


} // namespace sps
