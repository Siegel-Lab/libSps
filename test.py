import os
os.environ["STXXLLOGFILE"] = "/dev/null"
os.environ["STXXLERRLOGFILE"] = "/dev/null"

from build.libKdpsTree import *
import random

print_all = False


def fixed(tree, l, points, x=0, d=2):
    for idx, (pos, layer) in enumerate(points):
        tree.add_point(x, pos, layer, "p" + str(idx))
    tree.generate(x, 2, 0, len(points))
    if print_all:
        print(tree)
        print("generated")
    for (p1, layer1) in points:
        for (p2, layer2) in points + [([max(p[i] for (p, _) in points)+1 for i in range(d)], 0)]:
            if all(a < b for a, b in zip(p1, p2)):
                cnt = tree.count(x, p1, p2)
                itr = tree.get(x, p1, p2)
                truth = [sum(1 if all(i >= j and i < k for i, j,
                            k in zip(p, p1, p2)) and lx == ly else 0 for (p, ly) in points) for lx in range(l)]
                if not cnt == [len(x) for x in itr] == truth:
                    print("counts:", cnt, [len(x) for x in itr], truth)
                    print("query", p1, p2)
                    print("points", points)
                    print("iteration result", itr)
                    print(tree)
                    print("failure", l, x)
                    exit()
    print("success", l, x)


def test(tree, l, n=30):
    for x in range(1, n):
        for _ in range(min(x*2, 100)):
            tree.clear()
            points = []
            for _ in range(x):
                points.append(([random.choice(range(x)), random.choice(range(x))], random.choice(range(l))))
                if print_all:
                    print("adding", points[-1])
            fixed(tree, l, points, x)

def test_array(tree, l, n=30):
    for x in range(1, n):
        for _ in range(min(x*2, 100)):
            tree.clear()
            points = []
            for _ in range(x):
                points.append(([random.choice(range(x))], random.choice(range(l))))
                if print_all:
                    print("adding", points[-1])
            fixed(tree, l, points, x, d=1)


random.seed(6846854546132)
#fixed(KdpsTree_2D("test/blub2"), 2, [[0,1], [1,0], [1,2], [0,3], [1,4]])

test(KdpsTreeTest("test/blub1"), SETTINGS.NUM_LAYERS)
test(KdpsTree("test/blub2"), SETTINGS.NUM_LAYERS)
test_array(PsArray("test/blub3"), SETTINGS.NUM_LAYERS)
