#include "sps/main.h"
#include "sps/default.h"

int main( )
{
    typename sps::Main<InMemTypeDef<2, true, true>> xTree1( "test1" );
    typename sps::Main<CachedTypeDef<2, true, true>> xTree2( "test2" );
    typename sps::Main<DiskTypeDef<2, true, true>> xTree3( "test3" );

    return 0;
}