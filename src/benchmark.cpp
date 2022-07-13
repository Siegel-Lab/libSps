#include "sps/default.h"
#include "sps/index.h"
#include "sps/simple_vector.h"

int randNum(int max=100)
{
    return std::rand()/((RAND_MAX + 1u)/max);
}

int main( )
{
    using idx_t = typename sps::Index<InMemTypeDef<2, false, true, 0>> ;
    idx_t xIndex( "benchmark1", true );
    //typename sps::Index<CachedTypeDef<2, true, false, 1>> xTree2( "test2" );

    idx_t::coordinate_t uiZero = xIndex.numPoints();

    for(size_t uiI = 0; uiI < 1000; uiI++)
    {
        idx_t::ret_pos_t xAdd;
        for(size_t uiD = 0; uiD < std::tuple_size<idx_t::ret_pos_t>::value; uiD++)
            xAdd[uiD] = randNum();
        xIndex.addPoint(xAdd);
    }

    idx_t::class_key_t uiClass = xIndex.generate(uiZero, xIndex.numPoints());

    std::vector<std::pair<idx_t::ret_pos_t, idx_t::ret_pos_t>> vRegions;
    for(size_t uiI = 0; uiI < 1000; uiI++)
    {
        idx_t::ret_pos_t xA, xB;
        for(size_t uiD = 0; uiD < std::tuple_size<idx_t::ret_pos_t>::value; uiD++)
        {
            xA[uiD] = randNum();
            xB[uiD] = xA[uiD] + randNum();
        }
        vRegions.push_back(std::make_pair(xA, xB));
    }

    auto vRet = xIndex.countMultiple(uiClass, vRegions);

    return 0;
}