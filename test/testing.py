import random
import sps
import sys

class Point:
    def __init__(self, pos):
        self._pos = pos

    def __getitem__(self, key):
        return self._pos[key]

    def __setitem__(self, key, val):
        self._pos[key] = val

    def __iter__(self):
        yield from self._pos
    
    @staticmethod
    def random(data_space):
        return Point([random.randrange(f, t) for f, t in zip(data_space.from_pos(), data_space.to_pos())])

    def max(self, other):
        return Point([max(self_x, other_x) for self_x, other_x in zip(self, other)])

    def min(self, other):
        return Point([min(self_x, other_x) for self_x, other_x in zip(self, other)])

    def __str__(self):
        return str(self._pos)

    def __len__(self):
        return len(self._pos)

    def d(self):
        return len(self)

    def orto(self):
        return 0

    def to_list(self):
        return self._pos

    def add(self, v):
        return Point([x + v for x in self._pos])


class Intersection:
    ENCLOSES = 1
    ENCLOSED = 2
    OVERLAPS = 3
    LAST = 4
    FIRST = 5
    POINTS_ONLY = 6

    @staticmethod
    def random():
        return random.choice([#Intersection.ENCLOSES, #@todo bugged
                              Intersection.ENCLOSED,
                              Intersection.OVERLAPS,
                              Intersection.LAST, 
                              Intersection.FIRST,
                              Intersection.POINTS_ONLY
                              ])

    @staticmethod
    def to_name(intersection):
        return {Intersection.ENCLOSES: "ENCLOSES",
                Intersection.ENCLOSED: "ENCLOSED",
                Intersection.OVERLAPS: "OVERLAPS",
                Intersection.FIRST: "FIRST",
                Intersection.LAST: "LAST",
                Intersection.POINTS_ONLY: "POINTS_ONLY",
            }[intersection]


class Hyperrectangle:
    def __init__(self, from_pos, to_pos, num_orto):
        self._from_pos = from_pos
        self._to_pos = to_pos
        self._num_orto = num_orto

    def size(self):
        return [t - f for t, f in zip(self.to_pos(), self.from_pos())]

    def d(self):
        return self.from_pos().d()

    def from_pos(self):
        return self._from_pos

    def to_pos(self):
        return self._to_pos

    @staticmethod
    def square(d, size):
        return Hyperrectangle(Point([0]*d), Point([size]*d), d)

    @staticmethod
    def random(data_space, num_orto=None):
        if num_orto is None:
            num_orto = data_space.d()
        a = Point.random(data_space)
        b = Point.random(data_space)
        for idx in range(num_orto, data_space.d()):
            b[idx] = a[idx]
        return Hyperrectangle(a.min(b), a.max(b), num_orto)

    @staticmethod
    def point(data_space):
        a = Point.random(data_space)
        return Hyperrectangle(a, Point([x + 1 for x in a]), 0)

    def compare_dimension(self, other, intersection, d):
        if intersection == Intersection.ENCLOSES:
            return self.from_pos()[d] > other.from_pos()[d] and self.to_pos()[d] < other.to_pos()[d]
        if intersection == Intersection.ENCLOSED:
            return self.from_pos()[d] <= other.from_pos()[d] and self.to_pos()[d] > other.to_pos()[d]
        if intersection == Intersection.OVERLAPS:
            return self.to_pos()[d] > other.from_pos()[d] and self.from_pos()[d] <= other.to_pos()[d]
        if intersection == Intersection.LAST:
            return other.to_pos()[d] >= self.from_pos()[d] and other.to_pos()[d] < self.to_pos()[d]
        if intersection == Intersection.FIRST:
            return other.from_pos()[d] >= self.from_pos()[d] and other.from_pos()[d] < self.to_pos()[d]
        if intersection == Intersection.POINTS_ONLY:
            return all(x == 0 for x in other.size()) and self.to_pos()[d] > other.from_pos()[d] and \
                   self.from_pos()[d] <= other.from_pos()[d]

    def compare(self, other, intersection):
        ret = True
        for dx in range(self.d()):
            ret = ret and self.compare_dimension(other, intersection if dx < other.orto() else Intersection.ENCLOSED, dx)
        return ret

    def __str__(self):
        return "from: " + str(self.from_pos()) + " to: " + str(self.to_pos()) + " o: " + str(self.orto())

    def orto(self):
        return self._num_orto

    def get_corner(self, b_o_t):
        return Point([b if x else t for b, t, x in zip(self.from_pos(), self.to_pos(), b_o_t)])

class HyperrectangleValue(Hyperrectangle):
    def __init__(self, from_pos, to_pos, num_orto, value):
        super().__init__(from_pos, to_pos, num_orto)
        self._value = value

    def value(self):
        return self._value

    @staticmethod
    def gen(rect_gen, val_gen):
        h = rect_gen()
        return HyperrectangleValue(h.from_pos(), h.to_pos(), h.orto(), val_gen())

    def __str__(self):
        return super().__str__() + " value: " + str(self.value())

class SpsIndexWrapper:
    def __init__(self, index, d, o):
        try:
            self.index = sps.make_sps_index(num_dimensions=d, num_orthotope_dimensions=o)
        except Exception as e:
            print("requested num dimensions:", d)
            print("requested num orthotope dimensions:", o)
            raise e
        for x in index._data:
            if o > 0:
                self.index.add_point(x.from_pos().to_list(), x.to_pos().to_list(), x.value())
            else:
                self.index.add_point(x.from_pos().to_list(), x.value())
        self.idx = self.index.generate(verbosity=0)

    def _transl_inter(self, x):
        return {
            Intersection.ENCLOSES: sps.IntersectionType.encloses,
            Intersection.ENCLOSED: sps.IntersectionType.enclosed,
            Intersection.OVERLAPS: sps.IntersectionType.overlaps,
            Intersection.FIRST: sps.IntersectionType.first,
            Intersection.LAST: sps.IntersectionType.last,
            Intersection.POINTS_ONLY: sps.IntersectionType.points_only,
        }[x]

    def count(self, query, intersection=Intersection.ENCLOSES, verbosity=0):
        return self.index.count(self.idx, query.from_pos().to_list(), query.to_pos().to_list(), 
                         self._transl_inter(intersection), verbosity=verbosity)

    def __str__(self):
        return str(self.index)


class Index:
    def __init__(self, data):
        self._data = data

    def count(self, query, intersection=Intersection.ENCLOSES, verbosity=0):
        cnt = 0
        if verbosity > 0:
            print("query:", query, Intersection.to_name(intersection))
            print("data\tcount")
        for d in self._data:
            if query.compare(d, intersection):
                if verbosity > 0:
                    print(d, "\tcount:", d.value())
                cnt += d.value()
            else:
                if verbosity > 0:
                    print(d, "\tcount: -")
        if verbosity > 0:
            print("result:", cnt)
        return cnt

    @staticmethod
    def all_combinations(d):
        if d == 0:
            yield []
        else:
            for x in Index.all_combinations(d - 1):
                yield x + [True]
                yield x + [False]

    def prefix(self, query, intersection=Intersection.ENCLOSES):
        if not intersection in [Intersection.ENCLOSED, Intersection.FIRST, Intersection.LAST, Intersection.OVERLAPS]:
            print("WARNING: unimplemented for", Intersection.to_name(intersection))
        for b_o_t in Index.all_combinations(self.d()):
            cnt = 0
            q_pos = query.get_corner(b_o_t)
            for d in self._data:
                if intersection == Intersection.ENCLOSED:
                    d_pos = d.get_corner(b_o_t)
                if intersection == Intersection.FIRST:
                    d_pos = d.get_corner([True]*d.d())
                if intersection == Intersection.LAST:
                    d_pos = d.get_corner([False]*d.d())
                if intersection == Intersection.OVERLAPS:
                    d_pos = d.get_corner([not x for x in b_o_t])
                if all(d < q for q, d in zip(q_pos, d_pos)):
                    print("\t", d, "\tcount:", d.value())
                    cnt += d.value()
                else:
                    print("\t", d, "\tcount: -")
            if intersection == Intersection.ENCLOSED:
                print("prefix count for corner", q_pos.to_list() + query.size()[:self.orto()], "should be", cnt)
            if intersection in [Intersection.FIRST, Intersection.LAST, Intersection.OVERLAPS]:
                print("prefix count for corner", q_pos.to_list() + ["inf"]*self.orto(), "should be", cnt)


    @staticmethod
    def random(n, data_gen):
        data = []
        for _ in range(n):
            data.append(data_gen())
        return Index(data)

    
    def __str__(self):
        ret = ""
        for d in self._data:
            ret += str(d) + "\n"
        return ret

    def d(self):
        return self._data[0].d()
    
    def orto(self):
        return max(x.orto() for x in self._data)

    def to_sps_index(self, d, o):
        return SpsIndexWrapper(self, d, o)


class CountMultiple:
    def __init__(self, names, indices):
        self.names = names
        self.indices = indices

    def count(self, query, intersection=Intersection.ENCLOSES, count=0):
        reference_cnt = self.indices[0].count(query, intersection)
        for name, index in zip(self.names[1:], self.indices[1:]):
            cnt = index.count(query, intersection)

            if cnt != reference_cnt:
                print(name)
                index.count(query, intersection, verbosity=100)
                print()
                print()
                print(self.names[0])
                self.indices[0].count(query, intersection, verbosity=100)
                self.indices[0].prefix(query, intersection)
                print()
                print()
                print(name)
                print(index)
                print()
                print()
                print(self.names[0])
                print(self.indices[0])
                print()
                print()
                print(name, "got a different result than", self.names[0])
                print("expected:", reference_cnt, "but got:", cnt)
                print("query was", query, Intersection.to_name(intersection))
                print("d, o was:", self.indices[0].d(), self.indices[0].orto())
                print("seed was:", SEED)
                print("failure at attempt", count)
                exit()

def test_one(d=2, o=0, data_size=10, data_elements=3, query_elements=3, intersection=Intersection.random(), count=0):
    data_space = Hyperrectangle.random(Hyperrectangle.square(d, data_size))
    if min(data_space.size()) <= 1:
        data_space = Hyperrectangle.square(d, data_size)

    index = Index.random(data_elements,
                         lambda: HyperrectangleValue.gen(lambda: Hyperrectangle.random(data_space, num_orto=o), 
                                                         lambda: 1)
                            )
    #print("index", index, sep="\n")
    sps_index = index.to_sps_index(d, o)

    counter = CountMultiple(["py", "sps"], [index, sps_index])

    for _ in range(query_elements):
        query = Hyperrectangle.random(data_space)
        counter.count(query, intersection, count=count)

    print("success at attempt", count)


def random_n_from_s(d, o, s):
    area = s**d
    return [random.randrange(d, max(10, area))]# for _ in range(10)]

def intersection_from_o():
    return [Intersection.random() for _ in range(10)]

def data_sizes_exp():
    return [1,2,3,4,5,10,20,30,40,50,100,1000,10000]
    #return [x for i in range(6) for x in range(10**i, 6*(10**i), 10**i)]

def test_escalate(dos=[(2, 0)], 
                  data_sizes=data_sizes_exp,
                  data_elements=random_n_from_s,
                  query_elements=lambda d,o,s,n: [1],
                  intersections=intersection_from_o#lambda *x: [Intersection.ENCLOSED]
                  ):
    c = 1
    for s in data_sizes():
        for i in intersections():
            for d, o in dos:
                for n in data_elements(d, o, s):
                    for x in query_elements(d, o, s, n):
                        test_one(d, o, s, n, x, i, c)
                        c += 1


SEED = random.randrange(sys.maxsize)
SEED = 3232765277836967788
random.seed(SEED)

test_escalate([(2, 0), (1, 1)])