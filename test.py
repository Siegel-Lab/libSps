from build.libKdpsTree import *
import random

print_all = True


def fixed(tree, d, points, x=0):
    for idx, p in enumerate(points):
        tree.add_point(p, "p" + str(idx))
    tree.generate_for_points(x, 2, 0, len(points))
    if print_all:
        print(tree)
        print("generated")
    for p1 in points:
        for p2 in points + [[max(p[i] for p in points)+1 for i in range(len(points[0]))]]:
            if all(a < b for a, b in zip(p1, p2)):
                cnt = tree.count(x, p1, p2)
                itr = tree.get(x, p1, p2)
                truth = sum(1 if all(i >= j and i < k for i, j,
                            k in zip(p, p1, p2)) else 0 for p in points)
                if not cnt == len(itr) == truth:
                    print("counts:", cnt, len(itr), truth)
                    print("query", p1, p2)
                    print("points", points)
                    print("iteration result", itr)
                    print(tree)
                    print("failure", d, x)
                    exit()
    print("success", d, x)


def test(tree, d, n=30):
    for x in range(1, n):
        for _ in range(min(x*2, 100)):
            tree.clear()
            points = []
            for _ in range(x):
                points.append([])
                for _ in range(d):
                    points[-1].append(random.choice(range(x)))
                if print_all:
                    print("adding", points[-1])
            fixed(tree, d, points, x)


random.seed(6846854546132)
#fixed(KdpsTree_2D("test/blub2"), 2, [[0,1], [1,0], [1,2], [0,3], [1,4]])
test(KdpsTree("test/blub2"), 2)
