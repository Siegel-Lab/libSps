#pragma once

#include "sps/type_defs.h"
#include <cassert>
#include <functional>
#include <string>

#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif

namespace sps
{

template <typename type_defs> class SparseCoord
{
    EXTRACT_TYPE_DEFS; // macro call


    EXTRACT_VEC_GENERATOR( coord, coordinate_t ); // macro call

  public:
    coord_file_t xFile;
    coord_vec_t vData;

    struct Entry
    {
        coordinate_t uiStartIndex = std::numeric_limits<coordinate_t>::max();
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
            rSparseCoords.iterate([&](coordinate_t uiFrom, coordinate_t uiTo){
                if(!bFrist)
                    os << ", ";
                bFrist = false;
                os << uiFrom << " -> " << uiTo;
            }, *this);

            os << " }";

            return os;
        }

        std::string str() const
        {
            std::stringstream ss;
            ss << *this;
            return ss.str();
        }
    };

    SparseCoord( std::string sPrefix, bool bWrite )
        : xFile( coord_vec_generator.file( sPrefix + ".coords", bWrite ) ), vData( coord_vec_generator.vec( xFile ) )
    {}

    coordinate_t replace( coordinate_t uiX, const Entry& rInfo ) const
    {
        if( rInfo.uiStartIndex == std::numeric_limits<coordinate_t>::max( ) )
            return std::numeric_limits<coordinate_t>::max( );
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
        assert( uiX <= vData[ rInfo.uiStartIndex + rInfo.uiEndCord - rInfo.uiStartCord ] );
        auto xItBegin = vData.begin( ) + rInfo.uiStartIndex;
        auto xItEnd = vData.begin( ) + rInfo.uiStartIndex + 1 + rInfo.uiEndCord - rInfo.uiStartCord;
        // lowerbound can be used as search because indices must be continuous
        auto xIt = std::lower_bound(xItBegin, xItEnd, uiX );
        if(xIt == xItEnd)
            return std::numeric_limits<coordinate_t>::max();
        assert(*xIt == uiX);
        return (xIt - vData.begin( )) - rInfo.uiStartIndex + rInfo.uiStartCord;
    }

    template <size_t N> std::array<coordinate_t, N> axisSizes( const std::array<Entry, N>& vAxes ) const
    {
        std::array<coordinate_t, N> vAxisSizes;
        for( size_t uiI = 0; uiI < N; uiI++ )
            vAxisSizes[ uiI ] = replace( vAxes[ uiI ].uiEndCord, vAxes[ uiI ] ) + 1;
        return vAxisSizes;
    }

    template <size_t N>
    std::array<coordinate_t, N> sparse( const std::array<coordinate_t, N>& vCoords,
                                        const std::array<Entry, N>& vAxes ) const
    {
        std::array<coordinate_t, N> vRet;
        for( size_t uiI = 0; uiI < N; uiI++ )
            vRet[ uiI ] = replace( vCoords[ uiI ], vAxes[ uiI ] );
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

    template <typename Iterator_t> Entry addStart( Iterator_t xBegin, const Iterator_t& xEnd, 
                                                   coordinate_t uiStartWith, coordinate_t uiDiv = 1 )
    {
        Entry xRet{ };
        xRet.uiStartIndex = vData.size( );
        xRet.uiStartCord = uiStartWith;
        assert(!(xBegin != xEnd) || uiStartWith <= *xBegin);
        auto uiLast = uiStartWith;
        size_t uiI = 0;
        while( xBegin != xEnd )
        {
            for( coordinate_t uiX = uiLast; uiX < *xBegin; uiX++ )
                vData.push_back( uiI / uiDiv );
            if(uiLast < *xBegin)
                uiI++;
            uiLast = *xBegin;
            ++xBegin;
        }
        xRet.uiEndCord = uiLast;
        vData.push_back( uiI / uiDiv );
        assert(vData.size() - xRet.uiStartIndex == 1 + xRet.uiEndCord - xRet.uiStartCord);

        return xRet;
    }

    Entry addStart( coordinate_t uiStartWith )
    {
        Entry xRet{ };
        xRet.uiStartIndex = vData.size( );
        xRet.uiStartCord = uiStartWith;
        xRet.uiEndCord = uiStartWith;
        vData.push_back( 0 );

        return xRet;
    }

    template <typename Iterator_t> Entry add( Iterator_t xBegin, const Iterator_t& xEnd, coordinate_t uiDiv = 1 )
    {
        if( ! (xBegin != xEnd) )
            return Entry{};
        return addStart(xBegin, xEnd, *xBegin, uiDiv);
    }

    Entry add_vec(std::vector<size_t> vVec)
    {
        return add(vVec.begin(), vVec.end());
    }

    class EntryIterator
    {
        const SparseCoord& rCord;
        const Entry& rInfo;
        size_t uiI;


      public:
        EntryIterator( const SparseCoord& rCord, const Entry& rInfo ) : rCord( rCord ), rInfo( rInfo ), uiI( 0 )
        {
        }

        void operator++( )
        {
            if(uiI <= rInfo.uiEndCord - rInfo.uiStartCord)
            {
                coordinate_t uiLast = ( **this ).second;
                while( uiI <= rInfo.uiEndCord - rInfo.uiStartCord && ( **this ).second == uiLast )
                    uiI++;
            }
            else
                throw std::runtime_error("incrementing eof iterator");
        }

        const std::pair<coordinate_t, coordinate_t> operator*( ) const
        {
            assert(uiI <= rInfo.uiEndCord - rInfo.uiStartCord);
            return std::make_pair(uiI + rInfo.uiStartCord, rCord.vData[ uiI + rInfo.uiStartIndex ]);
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
        if( rInfo.uiStartIndex != std::numeric_limits<coordinate_t>::max() )
            xRet.uiI += 1 + rInfo.uiEndCord - rInfo.uiStartCord;
        return xRet;
    }

    void iterate( std::function<void( coordinate_t, coordinate_t )> fDo, const Entry& rInfo ) const
    {
        auto xIt = this->cbegin( rInfo );
        auto xItEnd = this->cend( rInfo );
        if(xIt != xItEnd)
            while( xIt != xItEnd )
            {
                fDo( (*xIt).first, (*xIt).second );
                ++xIt;
            }
        else
            fDo( rInfo.uiStartCord, 0 );
    }

    template <size_t I, size_t N>
    inline void iterateHelper( std::function<void( const std::array<coordinate_t, N>&, 
                                                   const std::array<coordinate_t, N>& )> fDo,
                   const std::array<Entry, N>& rInfos, std::array<coordinate_t, N>& rFrom,
                   std::array<coordinate_t, N>& rTo ) const
    {
        if constexpr(I == N)
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

    std::string str() const
    {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }

};


} // namespace sps


#if WITH_PYTHON
template <typename type_defs> void exportSparseCoord( pybind11::module& m, std::string sName )
{
    pybind11::class_<typename sps::SparseCoord<type_defs>::Entry>( m, (sName + "__Entry").c_str( ) )
        .def( "__str__", &sps::SparseCoord<type_defs>::Entry::str );

    pybind11::class_<sps::SparseCoord<type_defs>>( m, sName.c_str( ) )
        .def( pybind11::init<std::string, bool>( ), 
              pybind11::arg( "path" ),
              pybind11::arg( "write_mode" ) = false ) // constructor
        .def( "add", &sps::SparseCoord<type_defs>::add_vec )
        .def( "replace", &sps::SparseCoord<type_defs>::replace )
        .def( "inv_replace", &sps::SparseCoord<type_defs>::invReplace )
        .def( "__str__", &sps::SparseCoord<type_defs>::str )
        .def( "clear", &sps::SparseCoord<type_defs>::clear )
        ;
}
#endif
