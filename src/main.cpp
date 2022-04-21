#include "sps/default.h"
#include "sps/main.h"

int main()
{
    typename sps::Main<InMemTypeDef<2, true>> xTree1("test1");
    typename sps::Main<OnDiskTypeDef<2, true>> xTree2("test2");

    return 0;
}