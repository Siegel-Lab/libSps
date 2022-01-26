from build.libkdtree import *

threads = 32
cache = 10

bins = [range(4), range(4)]

d = [([1,1],"a"), ([2,2],"c"), ([2,2],"a")]

tree = kdtree_2("test/blub", d, bins, cache, threads, 0,0,0)

print(tree)

def assert_print(a, b):
    if a != b:
        print(a)
        assert False
    else:
        print("success")

assert_print(tree.count([0,0], [3,3]), 3)
assert_print(tree.count([2,2], [3,3]), 2)
assert_print(tree.count([0,0], [2,2]), 1)
assert_print(tree.count([0,0], [1,1]), 0)


bins = [range(3), range(3), range(3)]
d = [([1, 1, 1],"a")]

tree = kdtree_3("test/blub", d, bins, cache, threads, 0,0,0)

print(tree)


assert_print(tree.count([0,0,0], [2,2,2]), 1)


x = [x*2 for x in range(4)]
bins = [x, x]

d = [([1,1],"c"), ([4,4],"s"), ([5,5],"a")]

tree = kdtree_2("test/blub", d, bins, cache, threads, 0,0,0)

print(tree)
tree.save()


assert_print(tree.count([0,0], [6,6]), 3)
assert_print(tree.count([2,2], [4,4]), 0)
assert_print(tree.count([0,0], [4,4]), 1)
assert_print(tree.count([0,0], [1,1]), 0)
assert_print(tree.count([2,2], [4,4]), 0)
assert_print(tree.count([2,2], [3,3]), 0)
assert_print(tree.count([4,4], [5,5]), 1)

print(tree)

print("saving")

print("loading")
tree = kdtree_2("test/blub", cache, threads, 0,0,0)

print(tree)


assert_print(tree.count([0,0], [6,6]), 3)
assert_print(tree.count([2,2], [4,4]), 0)
assert_print(tree.count([0,0], [4,4]), 1)
assert_print(tree.count([0,0], [1,1]), 0)
assert_print(tree.count([2,2], [4,4]), 0)
assert_print(tree.count([2,2], [3,3]), 0)
assert_print(tree.count([4,4], [5,5]), 1) # @todo bug
exit()
