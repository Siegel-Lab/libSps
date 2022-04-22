#include "sps/default.h"
#include "sps/main.h"

int main()
{
    typename sps::Main<InMemTypeDef<2, true>> xTree1("test1");
    typename sps::Main<CachedTypeDef<2, true>> xTree2("test2");
    typename sps::Main<DiskTypeDef<2, true>> xTree3("test3");

    return 0;
}