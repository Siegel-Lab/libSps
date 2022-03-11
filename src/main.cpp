#include "kdpstree/default.h"
#include "kdpstree/tree.h"
#include "kdpstree/psarray.h"

int main(){
    typename kdpstree::Tree<InMemTypeDef<2>> xTree1("test1");
    typename kdpstree::Tree<OnDiskTypeDef<2>> xTree2("test2");

    typename kdpstree::Tree<InMemTypeDef<2>> xArray1("array1");
    typename kdpstree::Tree<OnDiskTypeDef<2>> xArray2("array2");

    return 0;
}