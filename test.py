from build.libcstree import *

print("a")
tree = CsTree_2D("test/blub")
print("b")

print(tree.str_raw())
print("c")
tree.add_point([0,0], "a")
tree.add_point([1,2], "b")
print("d")
print(tree.str_raw())
print(tree)
print("e")
