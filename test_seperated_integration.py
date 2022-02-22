import random
import bisect
from re import L

def rnd_point(d, max):
    p = []
    for _ in range(d):
        p.append(random.randrange(max))
    return p

def rnd_points(d, n):
    ret = []
    for _ in range(n):
        ret.append(rnd_point(d, n*2))
    return ret

def point_query(points, start, end):
    num = 0
    for p in points:
        inside = 1
        for a, b in zip(start, p):
            if a > b:
                inside = 0
        for a, b in zip(end, p):
            if a <= b:
                inside = 0
        num += inside
    return num

def make_index(points):
    ret = []
    for d in range(len(points[0])):
        ret.append([-10])
        for p in points:
            ret[-1].append(p[d])
        ret[-1].sort()
    return ret

def search_dim(l, p, inc):
    if inc:
        return bisect.bisect_left(l, p)-1
    else:
        return bisect.bisect_right(l, p-1)-1


def search_dims(idx, point, start):
    ret = float('inf')
    for l, p, s in zip(idx, point, start):
        ret = min(ret, search_dim(l, p, s == p))
    return ret

def search_area(idx, start, end, p=None, d=0, num_start=0):
    if p == None:
        p = [0]*len(start)
    if d == len(p):
        ret = search_dims(idx, p, start) * (1 if num_start % 2 == 0 else -1)
        print("search_area", ret, p, num_start)
        return ret
    else:
        ret = 0
        p[d] = start[d]
        ret += search_area(idx, start, end, p, d+1, num_start+1)
        p[d] = end[d]
        ret += search_area(idx, start, end, p, d+1, num_start)
        return ret

def test(idx, points):
    zero_pt = [0]*len(points[0])
    max_point = [max(p[i] for p in points)+1 for i in range(len(points[0]))]
    for start in points + [max_point, zero_pt]:
        for end in points + [max_point, zero_pt]:
            valid = True
            for a, b in zip(start, end):
                if a > b:
                    valid = False
            if valid:
                a = point_query(points, start, end)
                print("xxxx")
                b = search_area(idx, start, end)
                if a != b:
                    print("a, b", a, b)
                    print("points", points)
                    print("idx", idx)
                    print("start", start)
                    print("end", end)
                    assert False

def test_multiple():
    for n in range(1, 1000):
        for d in range(1, 5):
            points = rnd_points(d, n)
            idx = make_index(points)
            test(idx, points)


if __name__ == "__main__":
    test_multiple()