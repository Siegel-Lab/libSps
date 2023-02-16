import random
import sps

class Point:
    def __init__(self, pos):
        self._pos = pos

    def __getitem__(self, key):
        return self._pos[key]

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


class Intersection:
    ENCLOSES = 1
    ENCLOSED = 2
    OVERLAPS = 3
    LAST = 4
    FIRST = 5
    POINTS_ONLY = 6

class Hyperrectangle:
    def __init__(self, from_pos, to_pos):
        self._from_pos = from_pos
        self._to_pos = to_pos

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
        return Hyperrectangle([0]*d, [size]*d)

    @staticmethod
    def random(data_space):
        a = Point.random(data_space)
        b = Point.random(data_space)
        return Hyperrectangle(a.min(b), a.max(b))

    @staticmethod
    def point(data_space):
        a = Point.random(data_space)
        return Hyperrectangle(a, Point([x + 1 for x in a]))

    def compare_dimension(self, other, intersection, d):
        if intersection == Intersection.ENCLOSES:
            return self.from_pos()[d] <= other.from_pos()[d] and self.to_pos()[d] >= other.to_pos()[d]
        if intersection == Intersection.ENCLOSED:
            return self.from_pos()[d] >= other.from_pos()[d] and self.to_pos()[d] <= other.to_pos()[d]
        if intersection == Intersection.OVERLAPS:
            return self.to_pos()[d] >= other.from_pos()[d] and self.from_pos()[d] <= other.to_pos()[d]
        if intersection == Intersection.LAST:
            return self.to_pos()[d] >= other.to_pos()[d] and self.from_pos()[d] <= other.to_pos()[d]
        if intersection == Intersection.FIRST:
            return self.to_pos()[d] >= other.from_pos()[d] and self.from_pos()[d] <= other.from_pos()[d]
        if intersection == Intersection.POINTS_ONLY:
            return max(other.size()) == 1 and self.to_pos()[d] >= other.from_pos()[d] and self.from_pos()[d] <= other.from_pos()[d]

    def compare(self, other, intersection):
        ret = True
        for dx in range(self.d()):
            ret = ret and self.compare_dimension(other, intersection, dx)
        return ret

    def __str__(self):
        return "from: " + str(self.from_pos()) + " to: " + str(self.to_pos())

class HyperrectangleValue(Hyperrectangle):
    def __init__(self, from_pos, to_pos, value):
        super().__init__(from_pos, to_pos)
        self._value = value

    def value(self):
        return self._value

    @staticmethod
    def gen(rect_gen, val_gen):
        h = rect_gen()
        return HyperrectangleValue(h.from_pos(), h.to_pos(), val_gen())

    def __str__(self):
        return super().__str__() + " value: " + str(self.value())

class Index:
    def __init__(self, data):
        self._data = data

    def count(self, query, intersection=Intersection.ENCLOSES):
        cnt = 0
        for d in self._data:
            if query.compare(d, intersection):
                cnt += d.value()
        return cnt

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
        return self._data.d()

    def to_lib_sps_index(self):
        index = sps.make_sps_index(num_dimensions=)





def test(d=1, data_size=10, data_elements=3, query_elements=3):
    random.seed(6846854546132)

    data_space = Hyperrectangle.square(d, data_size)
    index = Index.random(data_elements,
                         lambda: HyperrectangleValue.gen(lambda: Hyperrectangle.random(data_space), 
                                                         lambda: 1)
                            )
    print("index", index, sep="\n")

    print("query", data_space)
    print("count", index.count(data_space))
    print()

    for _ in range(query_elements):
        query = Hyperrectangle.random(data_space)
        print("query", query)
        print("count", index.count(query, Intersection.OVERLAPS))
        print()

test()