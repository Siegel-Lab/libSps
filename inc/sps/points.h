#pragma once

#include "sps/desc.h"
#include "sps/point.h"
#include "sps/type_defs.h"
#include <cassert>
#include <functional>
#include <string>


namespace sps
{

template <typename type_defs> class Points
{
    EXTRACT_TYPE_DEFS; // macro call

    using point_t = Point<type_defs>;
    using desc_t = Desc<type_defs>;

  public:
    EXTRACT_VEC_GENERATOR_ELE( points, point_t, 4 * 1024 ); // macro call
  private:

    using points_it_t = typename points_vec_t::iterator;
    using const_points_it_t = typename points_vec_t::const_iterator;

    struct PointsComperator
    {
        const size_t uiDim;

        PointsComperator( size_t uiDim ) : uiDim( uiDim )
        {}

        bool operator( )( const point_t& a, const point_t& b ) const
        {
            return a.vPos[ uiDim ] < b.vPos[ uiDim ];
        }

        point_t min_value( ) const
        {
            return point_t( );
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
    points_file_t xFile;
    points_vec_t vData;
    
    struct Entry{
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

        coordinate_t size() const{
            return uiEndIndex - uiStartIndex;
        }
        
        std::ostream& stream( std::ostream& os, const Points& rPoints, const desc_t& vDesc ) const
        {
            os << "{ "; 
            for(size_t uiI = uiStartIndex; uiI < uiEndIndex; uiI++)
            {
                if(uiI > uiStartIndex)
                    os << ", ";
                rPoints.vData[uiI].stream(os, vDesc);
            }
            os << " }";

            return os;
        }

        std::ostream& stream( std::ostream& os, const Points& rPoints ) const
        {
            os << "{ "; 
            for(size_t uiI = uiStartIndex; uiI < uiEndIndex; uiI++)
            {
                if(uiI > uiStartIndex)
                    os << ", ";
                os << rPoints.vData[uiI];
            }
            os << " }";

            return os;
        }
    };
    
    class EntryIterator
    {
        const Points& rPoints;
        const Entry& rInfo;
        size_t uiI;
    public:
        EntryIterator(const Points& rPoints, const Entry& rInfo) :
            rPoints(rPoints),
            rInfo(rInfo),
            uiI(rInfo.uiStartIndex)
        {}

        void operator++()
        {
            assert(!eof());
            uiI++;
        }

        const point_t& operator*() const
        {
            return rPoints.vData[uiI];
        }

        const point_t* operator->() const
        {
            return &rPoints.vData[uiI];
        }

        bool operator!=(const EntryIterator& rOther) const
        {
            return uiI != rOther.uiI;
        }

        bool eof() const
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

        friend class Points;
    };

    Points( std::string sPrefix, bool bWrite )
        : xFile( points_vec_generator.file( sPrefix + ".points", bWrite ) ), vData( points_vec_generator.vec( xFile ) )
    {}

    void add( pos_t vPos, size_t uiDescOffset )
    {
        vData.push_back( point_t( vPos, uiDescOffset ) );
    }

    Entry getEntry() const
    {
        Entry xRet{};
        xRet.uiStartIndex = 0;
        xRet.uiEndIndex = vData.size( );
        return xRet;
    }

    EntryIterator cbegin(const Entry& rInfo) const
    {
        return EntryIterator(*this, rInfo);
    }

    EntryIterator cend(const Entry& rInfo) const
    {
        EntryIterator xRet(*this, rInfo);
        xRet.uiI = rInfo.uiEndIndex;
        return xRet;
    }

    void iterate( std::function<void( const point_t& )> fDo, const Entry& rEntry ) const
    {
        auto itEnd = vData.cbegin( ) + rEntry.uiEndIndex;
        for( auto cIter = vData.cbegin( ) + rEntry.uiStartIndex; cIter != itEnd; cIter++ )
            fDo( *cIter );
    }

    void sortByDim( size_t uiDim, const Entry& rEntry )
    {
        sort_points( vData.begin( ) + rEntry.uiStartIndex, vData.begin( ) + rEntry.uiEndIndex, 
                     PointsComperator( uiDim ) );
    }

    size_t size( ) const
    {
        return vData.size( );
    }

    void clear( )
    {
        vData.clear( );
    }

    friend std::ostream& operator<<( std::ostream& os, const Points& vPoints )
    {
        size_t uiX = 0;
        for( const auto& rP : vPoints.vData )
            os << uiX++ << ": " << rP << std::endl;
        return os;
    }
};


} // namespace sps
