from build.libcstree import *

print("a")
tree = CsTree_2D("test/blub", BinCordsGen())
print("b")

print(tree.str_raw())
print("c")
tree.add_point([0,0], "a")
tree.add_point([1,2], "b")
print("d")
tree.sort_points()
print("e")
tree.subdivide(1)
print("f")
print(tree.str_raw())
print(tree)
print("g")
