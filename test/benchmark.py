import sys
import os
sys.path.append(os.getcwd())
from build_benchmark.libSps import VERSION, make_sps_index, ValDiskSimpleVector_2D, ValCachedSimpleVector_2D
import random
import time
import os
import time


N_QUERY = 10000 # 10K
REPEAT_QUERY = 50

FILLS = [10000, 1000000] # 10K and 1G
AREAS = [10000, 1000000] # 10K and 1G

print("#VERSION", VERSION)
print("#N_QUERY", N_QUERY)
print("#REPEAT_QUERY", REPEAT_QUERY)
print("this will run for a while...")

files = [".prefix_sums", ".coords", ".overlays", ".datsets"]

def query(index, id, dims, genome_size, n=N_QUERY):
    ts = []
    for _ in range(REPEAT_QUERY):
        bins = []
        for _ in range(n):
            pos_s = []
            pos_e = []
            for _ in range(dims):
                x = random.randrange(genome_size)
                y = random.randrange(genome_size)
                pos_s.append(min(x, y))
                pos_e.append(max(x, y))
            bins.append((id, pos_s, pos_e))
        t1 = time.perf_counter()
        index.count_multiple(bins)
        t2 = time.perf_counter()
        ts.append(t2-t1)

    tm = sorted(ts)[len(ts)//2]
    # querytime
    print( ( n / tm ) / 1000, end="\t")

def fill(n, index, name, dims, is_ort, area, _area):
    print(*name, _area, n, sep="\t", end="\t")
    index.clear()
    t1 = time.perf_counter()
    for _ in range(n):
        if is_ort:
            pos_s = []
            pos_e = []
            for _ in range(0, 2):
                x = random.randrange(area)
                y = random.randrange(area)
                pos_s.append(min(x, y))
                pos_e.append(max(x, y))
            for _ in range(2, dims):
                pos_s.append(random.randrange(area))
            pos_e += pos_s[2:]
            index.add_point(pos_s, pos_e)
        else:
            pos_s = []
            for _ in range(0, dims):
                pos_s.append(random.randrange(area))
            index.add_point(pos_s)
    t2 = time.perf_counter()
    id = index.generate(verbosity=0)
    t3 = time.perf_counter()
    # index_name, filltime, generatetime
    print((t2-t1), (t3-t2), sep="\t", end="\t")
    return id

def make_indices():
    indices = []
    index_names = []
    params = []
    for dims in [2, 3]:
        for ort_dim in [False, True]:
            num_ort_dims = 2 if ort_dim else 0
            for bin_search in [False]:
            #for bin_search in [False, True]:
                #for storage in ["Disk"]:
                for storage in ["Disk", "Cached"]:
                    xs = [
                        str(dims) + "d", 
                        ("rect" if ort_dim else "point"), 
                        ("bin_search" if bin_search else "lookup_arr"), 
                        ("ram" if storage == "Disk" else "cached"), 
                    ]
                    index_names.append(xs)
                    indices.append(("test/benchmark_index", dims, num_ort_dims, bin_search, storage))
                    params.append((dims, ort_dim))
    return indices, index_names, params

def test():
    indices, index_names, params = make_indices()
    print("dims", "ort", "sparse_coords", 
          "storage", "area", "n", "filltime [s]", "generatetime [s]", "speed [queries / ms]", 
          "filesize [gb]", sep="\t")

    for index_params, name, param in zip(indices, index_names, params):
        index = make_sps_index(*index_params)
        for _area in AREAS:
            # d-th root of area
            area = int(_area ** ( 1 / (param[0] + param[1])))
            for n in FILLS:
                id = fill(n, index, name if n == FILLS[0] and _area == AREAS[0] else [""]*4, *param, area, _area)
                query(index, id, param[0], area)
                # filesize
                print(index.get_size(id)/ 10**9) # in GB
                index.clear()

        del index
        for file_suff in files:
            os.remove("test/benchmark_index" + file_suff)


def test_max_io(n_query=N_QUERY):
    print("checking maximal i/o of vector implementations")
    print("vector", "storage", "area", "speed [queries / ms]", sep="\t")
    for area in AREAS:
        vec_d = ValDiskSimpleVector_2D("test/benchmark_index_d", True)
        for _ in range(area):
            vec_d.add(random.choice(range(area)))
        t1 = 0
        t2 = 0
        for _ in range(REPEAT_QUERY):
            bins = []
            for _ in range(n_query):
                bins.append(random.choice(range(area)))
            t1 += time.perf_counter()
            vec_d.get_multiple(bins)
            t2 += time.perf_counter()
        print("vector", "ram", area, ( n_query * REPEAT_QUERY / (t2-t1) ) / 1000, sep="\t")

        del vec_d
        
        os.remove("test/benchmark_index_d.vals")

    for area in AREAS:
        vec_c = ValCachedSimpleVector_2D("test/benchmark_index_c", True)
        for _ in range(area):
            vec_c.add(random.choice(range(area)))
        ts = []
        for _ in range(REPEAT_QUERY):
            bins = []
            for _ in range(n_query):
                bins.append(random.choice(range(area)))
            t1 = time.perf_counter()
            vec_c.get_multiple(bins)
            t2 = time.perf_counter()
            ts.append(t2-t1)
            
        tm = sorted(ts)[len(ts)//2]
        print("", "cached", area, ( n_query / tm) / 1000, sep="\t")

        del vec_c
        
        os.remove("test/benchmark_index_c.vals")

random.seed(6846854546132)
test_max_io()
test()

# py test/benchmark.py