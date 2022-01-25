from build.libkdtree import kdtree

threads = 32

bins = [range(4), range(4)]

d = [([1,1],""), ([2,2],""), ([2,2],"")]

tree = kdtree(d, bins, 10, threads)

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
d = [([1, 1, 1],"")]

tree = kdtree(d, bins, 10, threads)

print(tree)


assert_print(tree.count([0,0,0], [2,2,2]), 1)


x = [x*2 for x in range(4)]
bins = [x, x]

d = [([1,1],""), ([4,4],""), ([5,5],"")]

tree = kdtree(d, bins, 10, threads)

print(tree)


assert_print(tree.count([0,0], [6,6]), 3)
assert_print(tree.count([2,2], [4,4]), 0)
assert_print(tree.count([0,0], [4,4]), 1)
assert_print(tree.count([0,0], [1,1]), 0)
assert_print(tree.count([2,2], [4,4]), 0)
assert_print(tree.count([2,2], [3,3]), 0)
assert_print(tree.count([4,4], [5,5]), 1)