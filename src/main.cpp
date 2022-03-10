#include "kdpstree/default.h"
#include "kdpstree/tree.h"

int main(){
    typename kdpstree::Tree<InMemTypeDef<1>> xTree1("test");
    typename kdpstree::Tree<InMemTypeDef<2>> xTree2("test");
    typename kdpstree::Tree<InMemTypeDef<3>> xTree3("test");
    typename kdpstree::Tree<InMemTypeDef<4>> xTree4("test");
    typename kdpstree::Tree<InMemTypeDef<5>> xTree5("test");

    return 0;
}