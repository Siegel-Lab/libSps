import os
os.environ["STXXLLOGFILE"] = "/dev/null"  # prevent creation of stxxl log files
os.environ["STXXLERRLOGFILE"] = "/dev/null" # prevent creation of stxxl log files

from build_dbg.libSps import *
import random

print_all = False


def fixed(arr, points):
    e = arr.add(sorted(points))
    d = {}
    for idx, p in enumerate(set(points)):
        d[p] = idx
    d[max(points) + 1] = idx
    for p in [*points, max(points)+1]:
        truth = d[p]
        x = arr.replace(p, e)
        y = p#arr.inv_replace(x, e)
        if truth != x or y != p:
            print(p, truth, x, y)
            print(str(e), str(arr))
            assert False 
    print("success")


def test(arr, n=30):
    for x in range(1, n):
        for _ in range(min(x*2, 100)):
            arr.clear()
            points = []
            for _ in range(x):
                points.append(random.choice(range(x)))
                if print_all:
                    print("adding", points[-1])
            fixed(arr, points)



random.seed(6846854546132)
#fixed(KdpsTree_2D("test/blub2"), 2, [[0,1], [1,0], [1,2], [0,3], [1,4]])

test(SparseCoords("test/sparse", True))
