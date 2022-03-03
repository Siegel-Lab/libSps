from build.libKdpsTree import *
import random


def test(tree, d, n=100):
    for x in range(1,n):
        for _ in range(x*2):
            tree.clear()
            points = []
            for _ in range(x):
                points.append([])
                for _ in range(d):
                    points[-1].append(random.choice(range(x)))
                print("adding", points[-1])
                tree.add_point(points[-1], "p" + str(len(points)))
            tree.generate_for_points(x, 2, 0, len(points))
            print(tree)
            print("generated")
            for p1 in points:
                for p2 in points + [[max(p[i] for p in points)+1 for i in range(len(points[0]))]]:
                    if all(a < b for a, b in zip(p1, p2)):
                        cnt = tree.count(x, p1, p2)
                        itr = tree.get(x, p1, p2)
                        truth = sum(1 if all(i >= j and i < k for i, j, k in zip(p, p1, p2)) else 0 for p in points)
                        if not cnt == len(itr) == truth:
                            print("counts:", cnt, len(itr), truth)
                            print("query", p1, p2)
                            print("points", points)
                            print("iteration result", itr)
                            print(tree)
                            print("failure", d, x)
                            exit()
            print("success", d, x)

random.seed(6846854546132)
test(KdpsTree_2D("test/blub2"), 2)
test(KdpsTree_3D("test/blub3"), 3)
test(KdpsTree_4D("test/blub4"), 4)
test(KdpsTree_5D("test/blub5"), 5)
