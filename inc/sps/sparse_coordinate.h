#pragma once

#include "sps/type_defs.h"
#include <cassert>
#include <functional>
#include <string>


namespace sps
{

template <typename type_defs> class SparseCoord
{
    EXTRACT_TYPE_DEFS; // macro call

    EXTRACT_VEC_GENERATOR( coord, coordinate_t ); // macro call

  public:
    static constexpr bool THREADSAVE = coord_THREADSAVE;
    std::mutex xResizeLock;
    coord_file_t xFile;
    coord_vec_t vData;
    size_t uiCapChangeLocks = 0;
    std::condition_variable xCapacityChangeVar;

    class CapacityChangeLock{
            SparseCoord& rCoords;
        public:
            CapacityChangeLock(SparseCoord& rCoords) : rCoords(rCoords)
            {
                std::lock_guard<std::mutex> xGuard( rCoords.xResizeLock );
                ++rCoords.uiCapChangeLocks;
            }

            ~CapacityChangeLock()
            {
                std::lock_guard<std::mutex> xGuard( rCoords.xResizeLock );
                assert(rCoords.uiCapChangeLocks > 0);
                --rCoords.uiCapChangeLocks;
                rCoords.xCapacityChangeVar.notify_one();
            }
    };

    std::shared_ptr<CapacityChangeLock> getCapacityGuard()
    {
        return std::make_shared<CapacityChangeLock>(*this);
    }

    size_t add_size(size_t uiAddSize)
    {
        std::unique_lock<std::mutex> xGuard( xResizeLock );
        if(vData.capacity( ) <= uiAddSize + vData.size( ))
        {
            assert(uiCapChangeLocks > 0);
            --uiCapChangeLocks;
            while(uiCapChangeLocks > 0 && vData.capacity( ) <= uiAddSize + vData.size( ))
                xCapacityChangeVar.wait(xGuard);
            ++uiCapChangeLocks;
            if(vData.capacity( ) <= uiAddSize + vData.size( ))
                vData.reserve( ( vData.size( ) + uiAddSize )*2 );
            xCapacityChangeVar.notify_all();
        }
        size_t uiRet = vData.size();
        vData.resize(uiRet + uiAddSize);
        return uiRet;
    }

    struct Entry
    {
        coordinate_t uiStartIndex = std::numeric_limits<coordinate_t>::max( );
        coordinate_t uiStartCord;
        coordinate_t uiEndCord;

        friend std::ostream& operator<<( std::ostream& os, const Entry& rEntry )
        {
            os << "i";
            os << rEntry.uiStartIndex;

            os << " s";
            os << rEntry.uiStartCord;

            os << " e";
            os << rEntry.uiEndCord;

            return os;
        }

        std::ostream& stream( std::ostream& os, const SparseCoord& rSparseCoords ) const
        {
            os << "{ ";
            bool bFrist = true;
            rSparseCoords.iterate(
                [ & ]( coordinate_t uiFrom, coordinate_t uiTo ) {
                    if( !bFrist )
                        os << ", ";
                    bFrist = false;
                    os << uiFrom << " -> " << uiTo;
                },
                *this );

            os << " }";

            return os;
        }

        std::string str( ) const
        {
            std::stringstream ss;
            ss << *this;
            return ss.str( );
        }
    };
    struct EntryArray : public Entry
    {
        coordinate_t uiNum = 0;

        friend std::ostream& operator<<( std::ostream& os, const EntryArray& rEntry )
        {
            Entry::operator<<( os, rEntry );
            os << " n";
            os << rEntry.uiNum;

            return os;
        }

        std::ostream& stream( std::ostream& os, const SparseCoord& rSparseCoords ) const
        {
            for( size_t uiI = 0; uiI < uiNum; uiI++ )
                SparseCoord::at( *this, uiI ).stream( os, rSparseCoords );
            return os;
        }
    };

    SparseCoord( std::string sPrefix, bool bWrite )
        : xFile( coord_vec_generator.file( sPrefix + ".coords", bWrite ) ), vData( coord_vec_generator.vec( xFile ) )
    {}

    static void append( EntryArray& rArr, const Entry& rE )
    {
        if( rArr.uiStartIndex == std::numeric_limits<coordinate_t>::max( ) )
        {
            rArr.uiStartIndex = rE.uiStartIndex;
            rArr.uiStartCord = rE.uiStartCord;
            rArr.uiEndCord = rE.uiEndCord;
        }
        assert( rArr.uiStartCord == rE.uiStartCord );
        assert( rArr.uiEndCord == rE.uiEndCord );
#ifndef NDEBUG
        coordinate_t uirArrEndIds = rArr.uiStartIndex + rArr.uiNum * ( 1 + rArr.uiEndCord - rArr.uiStartCord );
        assert( uirArrEndIds == rE.uiStartIndex );
#endif
        rArr.uiNum += 1;
    }
    static Entry at( const EntryArray& rArr, size_t uiI )
    {
        assert( uiI < rArr.uiNum );
        Entry rE;
        rE.uiStartIndex = rArr.uiStartIndex + uiI * ( 1 + rArr.uiEndCord - rArr.uiStartCord );
        rE.uiStartCord = rArr.uiStartCord;
        rE.uiEndCord = rArr.uiEndCord;
        return rE;
    }

    template <bool SANITY = true> inline coordinate_t replace( coordinate_t uiX, const Entry& rInfo ) const
    {
        if constexpr( SANITY )
            if( rInfo.uiStartIndex == std::numeric_limits<coordinate_t>::max( ) )
                return std::numeric_limits<coordinate_t>::max( );
        assert( vData.size( ) > rInfo.uiStartIndex + rInfo.uiEndCord - rInfo.uiStartCord );
        if( uiX < rInfo.uiStartCord )
            return std::numeric_limits<coordinate_t>::max( );
        if( uiX >= rInfo.uiEndCord )
            return vData[ rInfo.uiStartIndex + rInfo.uiEndCord - rInfo.uiStartCord ];
        return vData[ rInfo.uiStartIndex + uiX - rInfo.uiStartCord ];
    }

    coordinate_t invReplace( coordinate_t uiX, const Entry& rInfo ) const
    {
        if( rInfo.uiStartIndex == std::numeric_limits<coordinate_t>::max( ) )
            return std::numeric_limits<coordinate_t>::max( );
        assert( vData.size( ) > rInfo.uiStartIndex + rInfo.uiEndCord - rInfo.uiStartCord );
        if( uiX > vData[ rInfo.uiStartIndex + rInfo.uiEndCord - rInfo.uiStartCord ] )
            return std::numeric_limits<coordinate_t>::max( );
        if( uiX == std::numeric_limits<coordinate_t>::max( ) )
            return std::numeric_limits<coordinate_t>::max( );
        assert( uiX <= vData[ rInfo.uiStartIndex + rInfo.uiEndCord - rInfo.uiStartCord ] );
        auto xItBegin = vData.begin( ) + rInfo.uiStartIndex;
        auto xItEnd = vData.begin( ) + rInfo.uiStartIndex + 1 + rInfo.uiEndCord - rInfo.uiStartCord;
        // lowerbound can be used as search because indices must be continuous
        auto xIt = std::lower_bound( xItBegin, xItEnd, uiX );
        if( xIt == xItEnd )
            return std::numeric_limits<coordinate_t>::max( );
        assert( *xIt == uiX );
        return ( xIt - vData.begin( ) ) - rInfo.uiStartIndex + rInfo.uiStartCord;
    }

    coordinate_t axisSize( const Entry& rE ) const
    {
        return replace( rE.uiEndCord, rE ) + 1;
    }

    template <size_t N> std::array<coordinate_t, N> axisSizes( const std::array<Entry, N>& vAxes ) const
    {
        std::array<coordinate_t, N> vAxisSizes;
        for( size_t uiI = 0; uiI < N; uiI++ )
            vAxisSizes[ uiI ] = axisSize( vAxes[ uiI ] );
        return vAxisSizes;
    }


    template <size_t N, bool SANITY = true>
    inline __attribute__( ( always_inline ) ) std::array<coordinate_t, N>
    sparse( const std::array<coordinate_t, N>& vCoords, const std::array<Entry, N>& vAxes ) const
    {
        std::array<coordinate_t, N> vRet;
        for( size_t uiI = 0; uiI < N; uiI++ )
            vRet[ uiI ] = replace<SANITY>( vCoords[ uiI ], vAxes[ uiI ] );
        return vRet;
    }

    template <size_t N>
    std::array<coordinate_t, N> invSparse( const std::array<coordinate_t, N>& vCoords,
                                           const std::array<Entry, N>& vAxes ) const
    {
        std::array<coordinate_t, N> vRet;
        for( size_t uiI = 0; uiI < N; uiI++ )
            vRet[ uiI ] = invReplace( vCoords[ uiI ], vAxes[ uiI ] );
        return vRet;
    }

    template <bool CAPACITY_INC_ALLOWED = false, typename Iterator_t>
    Entry addStartEnd( Iterator_t xBegin, const Iterator_t& xEnd, coordinate_t uiStartWith,
                       coordinate_t uiEndWith = std::numeric_limits<coordinate_t>::max( ) )
    {
        Entry xRet{ };
        std::vector<coordinate_t> vTmp;
        assert( !( xBegin != xEnd ) || uiStartWith <= *xBegin );
        auto uiLast = uiStartWith;
        size_t uiI = 0;
        while( xBegin != xEnd )
        {
            for( coordinate_t uiX = uiLast; uiX < *xBegin; uiX++ )
                vTmp.push_back( uiI );
            if( uiLast < *xBegin )
                uiI++;
            uiLast = *xBegin;
            ++xBegin;
        }
        vTmp.push_back( uiI );
        if( uiEndWith != std::numeric_limits<coordinate_t>::max( ) )
            while( uiLast < uiEndWith )
            {
                vTmp.push_back( uiI );
                uiLast++;
            }


        assert( vTmp.size( ) == 1 + uiLast - uiStartWith );

        if constexpr(CAPACITY_INC_ALLOWED)
        {
            xRet.uiStartIndex = vData.size( );
            vData.resize( vTmp.size( ) + vData.size( ) );
        }
        else
            xRet.uiStartIndex = add_size(vTmp.size( ));


        xRet.uiStartCord = uiStartWith;
        xRet.uiEndCord = uiLast;

        // copy over the tmp vector
        for( size_t uiI = 0; uiI < vTmp.size( ); uiI++ )
            vData[ uiI + xRet.uiStartIndex ] = vTmp[ uiI ];

        return xRet;
    }


    template <bool CAPACITY_INC_ALLOWED = false> Entry addStart( coordinate_t uiStartWith )
    {
        Entry xRet{ };
        
        if constexpr(CAPACITY_INC_ALLOWED)
        {
            xRet.uiStartIndex = vData.size( );
            vData.push_back( 0 );
        }
        else
        {
            xRet.uiStartIndex = add_size(1);
            vData[xRet.uiStartIndex] = 0;
        }
        xRet.uiStartCord = uiStartWith;
        xRet.uiEndCord = uiStartWith;

        return xRet;
    }

    template <bool CAPACITY_INC_ALLOWED = false, typename Iterator_t>
    Entry add( Iterator_t xBegin, const Iterator_t& xEnd )
    {
        if( !( xBegin != xEnd ) )
            return Entry{ };
        return addStartEnd<CAPACITY_INC_ALLOWED>( xBegin, xEnd, *xBegin );
    }

    template <bool CAPACITY_INC_ALLOWED = false> Entry add_vec( std::vector<size_t> vVec )
    {
        return add<CAPACITY_INC_ALLOWED>( vVec.begin( ), vVec.end( ) );
    }

    class EntryIterator
    {
        const SparseCoord& rCord;
        const Entry& rInfo;
        size_t uiI;


      public:
        EntryIterator( const SparseCoord& rCord, const Entry& rInfo ) : rCord( rCord ), rInfo( rInfo ), uiI( 0 )
        {}

        void operator++( )
        {
            if( uiI <= rInfo.uiEndCord - rInfo.uiStartCord )
            {
                coordinate_t uiLast = ( **this ).second;
                while( uiI <= rInfo.uiEndCord - rInfo.uiStartCord && ( **this ).second == uiLast )
                    uiI++;
            }
            else
                throw std::runtime_error( "incrementing eof iterator" );
        }

        const std::pair<coordinate_t, coordinate_t> operator*( ) const
        {
            assert( uiI <= rInfo.uiEndCord - rInfo.uiStartCord );
            return std::make_pair( uiI + rInfo.uiStartCord, rCord.vData[ uiI + rInfo.uiStartIndex ] );
        }


        bool operator!=( const EntryIterator& rOther ) const
        {
            return uiI != rOther.uiI;
        }

        friend std::ostream& operator<<( std::ostream& os, const EntryIterator& rIt )
        {
            os << rIt.uiI;

            os << " ";
            os << rIt.rInfo;

            return os;
        }

        friend class SparseCoord;
    };

    EntryIterator cbegin( const Entry& rInfo ) const
    {
        return EntryIterator( *this, rInfo );
    }

    EntryIterator cend( const Entry& rInfo ) const
    {
        EntryIterator xRet( *this, rInfo );
        if( rInfo.uiStartIndex != std::numeric_limits<coordinate_t>::max( ) )
            xRet.uiI += 1 + rInfo.uiEndCord - rInfo.uiStartCord;
        return xRet;
    }

    void iterate( std::function<void( coordinate_t, coordinate_t )> fDo, const Entry& rInfo ) const
    {
        auto xIt = this->cbegin( rInfo );
        auto xItEnd = this->cend( rInfo );
        // if(xIt != xItEnd)
        while( xIt != xItEnd )
        {
            fDo( ( *xIt ).first, ( *xIt ).second );
            ++xIt;
        }
        // else
        //     fDo( rInfo.uiStartCord, 0 );
    }

    template <size_t I, size_t N>
    inline void
    iterateHelper( std::function<void( const std::array<coordinate_t, N>&, const std::array<coordinate_t, N>& )> fDo,
                   const std::array<Entry, N>& rInfos, std::array<coordinate_t, N>& rFrom,
                   std::array<coordinate_t, N>& rTo ) const
    {
        if constexpr( I == N )
            fDo( rFrom, rTo );
        else
            iterate(
                [ & ]( coordinate_t uiFrom, coordinate_t uiTo ) {
                    rFrom[ I ] = uiFrom;
                    rTo[ I ] = uiTo;
                    iterateHelper<I + 1, N>( fDo, rInfos, rFrom, rTo );
                },
                rInfos[ I ] );
    }


    template <size_t N>
    void iterate( std::function<void( const std::array<coordinate_t, N>&, const std::array<coordinate_t, N>& )> fDo,
                  const std::array<Entry, N>& rInfos ) const
    {
        std::array<coordinate_t, N> rFrom;
        std::array<coordinate_t, N> rTo;
        iterateHelper<0, N>( fDo, rInfos, rFrom, rTo );
    }

    void clear( )
    {
        vData.clear( );
    }

    friend std::ostream& operator<<( std::ostream& os, const SparseCoord& rCoords )
    {
        os << rCoords.vData;

        return os;
    }

    std::string str( ) const
    {
        std::stringstream ss;
        ss << *this;
        return ss.str( );
    }
};


} // namespace sps
