from build.libkdtree import kdtree


bins = [range(4), range(4)]

d = [([1,1],""), ([2,2],""), ([2,2],"")]

tree = kdtree(d, bins)

print(tree)

def assert_print(a, b):
    if a != b:
        print(a)
        assert False

assert_print(tree.count([0,0], [3,3]), 3)
assert_print(tree.count([0,0], [2,2]), 3)
assert_print(tree.count([0,0], [1,1]), 1)
assert_print(tree.count([0,0], [0,0]), 0)


bins = [range(3), range(3), range(3)]
d = [([1,1, 1],"")]

tree = kdtree(d, bins)

print(tree)

def assert_print(a, b):
    if a != b:
        print(a)
        assert False

assert_print(tree.count([0,0,0], [3,3,3]), 1)