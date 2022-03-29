from build_dbg.libSps import *
import random

print_all = True


def fixed(tree, points, d=2, cont=0):
    for idx, pos in enumerate(points):
        tree.add_point(pos, "p" + str(idx))
    x = tree.generate(0, len(points))
    print("done generating")
    if print_all:
        print(tree)
        print("generated")
    for p1 in points:
        for p2 in points + [[max(p[i] for p in points)+1 for i in range(d)]]:
            if all(a < b for a, b in zip(p1, p2)):
                cnt = tree.count(x, p1, p2)
                truth = sum(1 if all(i >= j and i < k for i, j, k in zip(p, p1, p2)) else 0 for p in points)
                itr = list(range(truth))#tree.get(x, p1, p2)
                if not cnt == len(itr) == truth:
                    print("counts:", cnt, len(itr), truth)
                    print("query", p1, p2)
                    print("points", points)
                    print("iteration result", itr)
                    print(tree)
                    print("failure", cont)
                    exit()
                else:
                    if print_all:
                        print("OK")
            else:
                if print_all:
                    print("not valid")
    print("success", cont)


def test(tree, d, n=30):
    cont = 0
    for x in range(1, n):
        for _ in range(min(x*2, 100)):
            tree.clear()
            points = []
            for _ in range(x):
                pos = []
                for _ in range(d):
                    pos.append(random.choice(range(x)))
                points.append(pos)
                if print_all:
                    print("adding", points[-1])
            fixed(tree, points, d, cont)
            cont += 1



random.seed(6846854546132)
#fixed(KdpsTree_2D("test/blub2"), 2, [[0,1], [1,0], [1,2], [0,3], [1,4]])

test(SparsePrefixSum_2D("test/blub1", True), 2)
