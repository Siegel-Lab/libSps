#include "sps/index.h"
#include "sps/default.h"
#include "sps/simple_vector.h"

int main( )
{
    typename sps::Index<InMemTypeDef<2, true, 1>> xTree1( "test1" );
    typename sps::Index<CachedTypeDef<2, true, 1>> xTree2( "test2" );
    typename sps::Index<DiskTypeDef<2, true, 1>> xTree3( "test3" );

    typename sps::SimpleVector<InMemTypeDef<2, false, 2>> xTree4( "test4" );

    return 0;
}