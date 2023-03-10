import sys
import os
sys.path.append(os.getcwd())
from build_benchmark.libSps import VERSION, make_sps_index
import random
import time

def create_index(d, o, n, index=None):
    if index is None:
        index = make_sps_index("test/benchmark_index2", d, o, False, "Ram", True)
    for _ in range(n):
        pos_s = []
        for _ in range(0, d):
            pos_s.append(random.randrange(n))
        index.add_point(pos_s)
    id = index.generate(verbosity=0)
    return id, index

def query_index(id, index, d, n):
    qs = []
    for _ in range(n):
        pos_s = []
        pos_e = []
        for _ in range(0, d):
            a = random.randrange(n)
            b = random.randrange(n)
            pos_s.append(min(a,b))
            pos_e.append(max(a,b))
        qs.append((id, pos_s, pos_e))
    index.count_multiple(qs)

index = None
td = 0
while td < 2:
    t0 = time.perf_counter()
    for _ in range(10):
        id, index = create_index(3, 0, random.randrange(500), index)
        #index.reserve(5000, 500*3*1000, 1000*1000*1000, 2000, 1000)
    t1 = time.perf_counter()
    td = t1 - t0
    print(td)
#query_index(id, index, 3, 1000000)
t2 = time.perf_counter()

print("took", t1 - t0, "and", t2-t1)


# perf record -g -e cpu-clock python3 test/perf_benchmark.py
# perf script -F +pid | sed 's/sps::TypeDefs<unsigned long, unsigned int, 3ul, unsigned short, RamVecGenerator, RamVecGenerator, RamVecGenerator, RamVecGenerator, RamVecGenerator, RamVecGenerator, false, 0ul, false, StdOutProgressStream>/XXX/g' - > firefox.perf
# perf report --hierarchy --stdio --stdio-color always | sed 's/sps::TypeDefs<unsigned long, unsigned int, 3ul, unsigned short, RamVecGenerator, RamVecGenerator, RamVecGenerator, RamVecGenerator, RamVecGenerator, RamVecGenerator, false, 0ul, false, StdOutProgressStream>/XXX/g' - | less -r

# profiler.firefox.com