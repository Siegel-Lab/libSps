#include "sps/index.h"
#include "sps/default.h"

int main( )
{
    typename sps::Index<InMemTypeDef<2, true, true>> xTree1( "test1" );
    typename sps::Index<CachedTypeDef<2, true, true>> xTree2( "test2" );
    typename sps::Index<DiskTypeDef<2, true, true>> xTree3( "test3" );

    return 0;
}