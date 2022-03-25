#include "kdpstree/default.h"
#include "kdpstree/main.h"

int main()
{
    typename kdpstree::Main<InMemTypeDef<2>> xTree1("test1");
    typename kdpstree::Main<OnDiskTypeDef<2>> xTree2("test2");

    return 0;
}