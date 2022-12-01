#include "sps/default.h"
#include "sps/index.h"
#include "sps/simple_vector.h"

int main( )
{
    typename sps::Index<InMemTypeDef<2, 1, true>> xTree1( "test1" );
    typename sps::Index<CachedTypeDef<2, 1, false>> xTree2( "test2" );

    typename sps::SimpleVector<InMemTypeDef<2, 2, false>> xTree4( "test4" );

    return 0;
}