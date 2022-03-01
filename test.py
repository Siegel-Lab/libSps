from build.libKdpsTree import *
import random


def test(tree, d, n=100):
    for x in range(1,n):
        points = []
        for _ in range(x):
            points.append([])
            for _ in range(d):
                points[-1].append(random.choice(range(n)))
            print("adding", points[-1])
            tree.add_point(points[-1], "p" + str(len(points)))
        tree.generate_for_points(0, 5, 0, len(points))
        for p1 in points:
            for p2 in points:
                if all(a < b for a, b in zip(p1, p2)):
                    cnt = tree.count(0, p1, p2)
                    itr = tree.get(0, p1, p2)
                    truth = sum(1 if all(x >= p1 and x < p2 for x in p) else 0 for p in points)
                    if not cnt == len(itr) == truth:
                        print("!", cnt, len(itr), truth)
                        print(points)
                        print(itr)
                        print(tree)
                        exit()
        print("sucess", x)

test(KdpsTree_2D("test/blub"), 2)
