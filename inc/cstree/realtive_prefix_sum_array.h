#pragma once

#if WITH_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif
#include <string>


namespace rpsa
{
template <typename vec_t> class DataRpsa
{
  protected:
    vec_t vData;
    using val_t = typename vec_t::value_type;

    virtual size_t lastIndex( )
    {
        return 0;
    }
};

template <typename vec_t> class DimRpsa : public DataRpsa<vec_t>
{
  protected:
    val_t& numDimensions( )
    {
        return vData[ 0 ];
    }

    val_t& getSizeOfDimension( size_t uiD )
    {
        assert( uiD < numDimensions( ) );
        return vData[ 1 + uiD ];
    }

    virtual size_t lastIndex( )
    {
        return DataRpsa::lastIndex( ) + numDimensions( ) + 1;
    }
};

template <typename vec_t> class BinRpsa : public DimRpsa<vec_t>
{
  protected:
    size_t getSizeOfBinsInDimension( size_t uiD )
    {
        return (val_t)std::sqrt( getSizeOfDimension( uiD ) );
    }

    size_t getDataSizeOfBin( )
    {
        size_t uiRet = 1;
        for( size_t uiI = 0; uiI < numDimensions( ); uiI++ )
            uiRet += getSizeOfBinsInDimension( uiI );
        return uiRet;
    }

    size_t getNumberOfBinsAlongDimension( size_t uiD )
    {
        return (val_t)( getSizeOfBinsInDimension( uiD ) / getSizeOfDimension( uiD ) ) + 1;
    }

    size_t getNumberOfBins( )
    {
        size_t uiRet = 1;
        for( size_t uiI = 0; uiI < numDimensions( ); uiI++ )
            uiRet *= getNumberOfBinsAlongDimension( uiI );
        return uiRet;
    }

    size_t getNumberOfBin( std::vector<size_t> vIndexOfBinPerDimension )
    {
        assert( vIndexOfBinPerDimension.size( ) == numDimensions( ) );
        size_t uiRet = 0;
        for( size_t uiI = 0; uiI < numDimensions( ); uiI++ )
            uiRet = uiRet * getNumberOfBinsAlongDimension( uiI ) + vIndexOfBinPerDimension[ uiI ];
        return uiRet;
    }

    size_t getIndexOfBin( std::vector<size_t> vIndexOfBinPerDimension )
    {
        return getNumberOfBin( vIdx ) * getDataSizeOfBin( ) + DimRpsa::lastIndex( );
    }

    val_t& getContSumPointerOfBin( size_t uiBinIdx )
    {
        return vData[ uiBinIdx ];
    }

    val_t& getOverlayPrefixSumFromBin( size_t uiBinIdx, size_t uiDim, size_t uiOverlayPrefixSumIdx )
    {
        size_t uiIdx = uiBinIdx + uiOverlayPrefixSumIdx + 1;
        for( size_t uiI = 0; uiI < uiDim; uiI++ )
            uiIdx += getSizeOfBinsInDimension( uiI );
        return vData[ uiIdx ];
    }

    std::vector<size_t> getBinIndexFromPos( std::vector<val_t> vPos )
    {
        assert( vPos.size( ) == numDimensions( ) );
        std::vector<size_t> vRet;
        for( size_t uiI = 0; uiI < numDimensions( ); uiI++ )
            vRet.push_back( vPos[ uiI ] / getSizeOfBinsInDimension( uiI ) );
        return vRet;
    }

    std::vector<size_t> getOverlayPrefixSumIndexFromPosition( std::vector<val_t> vPos )
    {
        assert( vPos.size( ) == numDimensions( ) );
        std::vector<size_t> vRet;
        for( size_t uiI = 0; uiI < numDimensions( ); uiI++ )
            vRet.push_back( vPos[ uiI ] % getSizeOfBinsInDimension( uiI ) );
        return vRet;
    }

    virtual size_t lastIndex( )
    {
        return DimRpsa::lastIndex( ) + getDataSizeOfBin( ) * getNumberOfBins( );
    }
};

template <typename vec_t> class ContSumRpsa : public BinRpsa<vec_t>
{
  protected:
    val_t& getContSumCoordinate( size_t uiContSumIdx, size_t uiD )
    {
        return vData[uiContSumIdx + 2 + uiD];
    }

    size_t contSumBinarySearch( size_t uiContSumStart, size_t uiContSumEnd )
    {}
};
} // namespace rpsa