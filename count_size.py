import enum
import random
import math
from re import S
from statistics import median


def parse_heatmap(in_filename, compression):
    with open(in_filename, "r") as in_file_1:
        for line in in_file_1:
            # parse file columns
            read_name, strnd_1, chr_1, pos_1, _, strnd_2, chr_2, pos_2, _2, mapq_1, mapq_2 = line.split()
            # convert number values to ints
            pos_1, pos_2, mapq_1, mapq_2 = (
                int(x) for x in (pos_1, pos_2, mapq_1, mapq_2))
            pos_1 -= 1
            pos_2 -= 1

            yield chr_1, pos_1//compression, chr_2, pos_2//compression, max(mapq_1, mapq_2)


def input_to_points(len_file_name, heatmap_filenames, limit=float('inf'), compression=1):
    curr_start = 0
    chr_start_pos = {}
    with open(len_file_name, "r") as len_file:
        for line in len_file:
            chr_name, chr_len = line.split("\t")
            chr_start_pos[chr_name] = curr_start
            curr_start += int(chr_len)//compression
    points = []
    cnt = 0
    for idx, file_name in enumerate(heatmap_filenames):
        x = []
        for _ in range(11):
            x.append([])
        points.append(x)
        for chr_1, pos_1, chr_2, pos_2, map_q in parse_heatmap(file_name, compression):
            x = chr_start_pos[chr_1] + int(pos_1)
            y = chr_start_pos[chr_2] + int(pos_2)
            max_dist_betw_alignments = random.choice(list(range(10)))
            a = random.choice(range(10))
            points[-1][a+1].append((map_q, max_dist_betw_alignments, x, y))
            points[-1][0].append((map_q, max_dist_betw_alignments, x, y))
            cnt += 1
            if cnt > limit:
                return points, False
    return points, True


def product(l):
    ret = 1
    for x in l:
        ret *= x
    return ret


def best_piv_dim(points):
    best_d = 0
    split = 0
    for d in range(len(points[0])):
        s = len(set(p[d] for p in points))
        if s > split:
            split = s
            best_d = d
    return best_d


def pivot(points):
    best_d = 0
    split = set()
    for d in range(len(points[0])):
        s = set(p[d] for p in points)
        if len(s) > len(split):
            split = s
            best_d = d
    return sorted(split)[len(split)//2], best_d


def split_points(points, axes):
    piv, d = pivot(points)
    a = []
    b = []
    for p in points:
        if p[d] < piv:
            a.append(p)
        else:
            b.append(p)
    #print("split", len(points), "->", len(a), len(b))
    a_2 = []
    b_2 = []
    for i, ax in enumerate(axes):
        if i == d:
            a_2.append([])
            b_2.append([])
            for p in ax:
                if p < piv:
                    a_2[-1].append(p)
                else:
                    b_2[-1].append(p)
        else:
            a_2.append(ax)
            b_2.append(ax)
    return (a, a_2), (b, b_2)


def mem(points, axes):
    s = []
    for _ in range(len(axes)):
        s.append(set())
    for p in points:
        for i, x in enumerate(p):
            s[i].add(x)
    #print(points, axes)
    #assert (len(points) > 0) == (len(axes) > 0)
    # the position in the dimension and the continuous sum
    # also the bottom left position of the bin
    overlay_size = sum(len(x)-1 for x in axes) * 2 + len(axes) + 2 + 1
    #overlay_size_2 = sum(max(x)-min(x) if len(x) > 0 else 0 for x in axes) + len(axes) * 2 + 1
    overlay_fill = sum(len(x) for x in axes) / max(1,sum(max(x)-min(x) if len(x) > 0 else 0 for x in axes))
    # the position = d-long vec, the cont sum value and the pointer to the description
    intersections = product(len(x) for x in s)//2 + 1 + len(set(points))
    intersections_size = intersections * (len(axes) + 1)
    points_size = len(points) * 2
    return overlay_size, points_size, intersections_size, intersections, overlay_fill


def compute_mem(points, l=None, axes=None):
    if l is None:
        l = math.log2(len(points))*100
        print("max points per overlay", l)
        #l = math.pow(len(points), 1 / len(points[0]))
    if axes is None:
        axes = []
        for d in range(len(points[0])):
            axes.append(list(set(p[d] for p in points)))
    if len(set(points)) > 1 and l < len(points):
        split = split_points(points, axes)
        del points
        del axes
        split_overlay_size = 0
        split_points_size = 0
        split_points_sum = 0
        split_intersections_sum = 0
        num_overlays = 0
        split_intersections = 0
        split_overlay_fill = 0
        for points_split, split_axis in split:
            a, b, g, c, d, e, f = compute_mem(points_split, l, split_axis)
            split_overlay_size += a
            split_points_size += b
            split_intersections_sum += g
            split_points_sum += d
            num_overlays += c
            split_intersections += e
            split_overlay_fill += f
        return split_overlay_size, split_points_size, split_intersections_sum, num_overlays, \
               split_points_sum/len(split), split_intersections/len(split), split_overlay_fill/len(split)
    else:
        overlay_size, points_size, intersections_size, intersections, overlay_fill = mem(
            points, axes)
        return overlay_size, points_size, intersections_size, 1, len(points), intersections, overlay_fill

def sum_mem(points_l):
    l = []
    for ps in points_l:
        for p in ps:
            l.append(compute_mem(p))
    r = [sum(p[x] for p in l) for x in range(len(l[0]))]
    for x in range(4,7):
        r[x] = r[x] / len(l)
    return r

EXMAPLE = [[[
    (1, 1),
    (2, 2),
    (3, 3),
]]]
EXMAPLE_2 = [[[
    (1, 1),
    (3, 1),
    (2, 3),
    (3, 4),
    (4, 3),
]]]

BED_FOLDER = "/work/project/ladsie_012/ABS.2.2/2021-10-26_NS502-NS521_ABS_CR_RADICL_inputMicroC/bed_files"
BED_SUFFIX = "RNA.sorted.bed_K1K2.bed_K4.bed_R_D.bed_R_D_K1K2.bed_R_D_PRE1.bed"
FILES = [
    BED_FOLDER + "/NS504_P10_Total_3." + BED_SUFFIX,
    BED_FOLDER + "/NS505_N50_Total_1." + BED_SUFFIX,
    BED_FOLDER + "/NS508_P10_NPM_1." + BED_SUFFIX,
    BED_FOLDER + "/NS511_N50_NPM_1." + BED_SUFFIX,
]

if __name__ == "__main__":
    #sum_mem(EXMAPLE_2)
    #exit()
    for x in range(3, 100):
        print(10**x, "points:")
        points, have_all = input_to_points(
            "../out/Lister427.sizes", FILES, limit=10**x)
        a, b, g, c, d, e, f = sum_mem(points)
        print("\ttotal of", a, "overlay integers = ", 100*a/(a+b+g), "%")
        print("\ttotal of", b, "points integers = ", 100*b/(a+b+g), "%")
        print("\ttotal of", g, "intersection integers = ", 100*g/(a+b+g), "%")
        print("\ttotal of", c, "overlays ->" , math.log2(c), "search complexity")
        print("\taverage of", d, "points per overlay")
        print("\taverage of", e, "intersections per overlay")
        print("\taverage overlay fill: ", 100*f, "%")
        print("\ttotal of", 8 * (a+b+g) / 10**9, "gb")
        print("\tOR total of", 8 * (a+b) / 10**9, "gb with runtime penalty of factor", 1+d/math.log2(d+e),
                "memory reduced to:", 100*((a+b)/(a+b+g)), "%")
        if have_all:
            break

 # 2^(xd)(dN+1)+x^d(N/2^(xd)-1)^d
