import sys
import os
sys.path.append(os.getcwd())
from build_benchmark.libSps import VERSION, make_sps_index
import random

def create_index(d, o, n):
    index = make_sps_index("test/benchmark_index2", d, o, "Ram", True)
    for _ in range(n):
        pos_s = []
        for _ in range(0, d):
            pos_s.append(random.randrange(n))
        index.add_point(pos_s)
    id = index.generate()
    return id, index


create_index(3, 0, 1000)


# perf record -g -e cpu-clock python3 test/perf_benchmark.py
# perf report --hierarchy --stdio