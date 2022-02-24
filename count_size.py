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


def input_to_points(len_file_name, heatmap_filenames, limit=float('inf'), compression=100):
    curr_start = 0
    chr_start_pos = {}
    with open(len_file_name, "r") as len_file:
        for line in len_file:
            chr_name, chr_len = line.split("\t")
            chr_start_pos[chr_name] = curr_start
            curr_start += int(chr_len)//compression
    points_up = []
    points_right = []
    cnt = 0
    for file_name in heatmap_filenames:
        for chr_1, pos_1, chr_2, pos_2, map_q in parse_heatmap(file_name, compression):
            x = chr_start_pos[chr_1] + int(pos_1)
            y = chr_start_pos[chr_2] + int(pos_2)
            max_dist_betw_alignments = random.choice(list(range(10)))
            annotations = [True] + \
                [random.choice([True, False]) for _ in range(10)]
            p = points_up if y > x else points_right
            for a, _ in enumerate(annotations):
                p.append((map_q, max_dist_betw_alignments, a, x, y))
            cnt += 1
            if cnt > limit:
                return points_up, points_right, False
    return points_up, points_right, True


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
    # the position in the dimension and the continuous sum
    # also the bottom left position of the bin
    #overlay_size = sum(len(x) for x in axes) * 2 + len(axes)
    overlay_size = sum(max(x)-min(x) for x in axes) + len(axes) * 2
    overlay_fill = sum(len(x) for x in axes) / sum(max(x)-min(x) for x in axes)
    # the position = d-long vec, the cont sum value and the pointer to the description
    intersections = product(len(x) for x in s)//2 + 1
    cont_sum_size = (intersections + len(set(points))) * (len(axes) + 2)
    return overlay_size, cont_sum_size, intersections, overlay_fill


def compute_mem(points, axes=None):
    if axes is None:
        axes = []
        for d in range(len(points[0])):
            axes.append(list(set(p[d] for p in points)))
    if len(points) == 0:
        return 0, 0, 0
    overlay_size, cont_sum_size, intersections, overlay_fill = mem(
        points, axes)
    split = split_points(points, axes)
    split_size = 0
    for points_split, split_axis in split:
        split_overlay_size, split_cont_sum_size, _, _ = mem(
            points_split, split_axis)
        split_size += split_overlay_size + split_cont_sum_size
    if overlay_size + cont_sum_size <= split_size:
        #start = [min(p[d] for p in points) for d in range(len(points[0]))]
        #end = [max(p[d] for p in points) for d in range(len(points[0]))]
        #print("overlay with", len(points), "points")
        return overlay_size, cont_sum_size, 1, len(points), intersections, overlay_fill
    else:
        del points
        del axes
        split_overlay_size = 0
        split_cont_sum_size = 0
        split_points_sum = 0
        num_overlays = 0
        split_intersections = 0
        split_overlay_fill = 0
        for points_split, split_axis in split:
            a, b, c, d, e, f = compute_mem(points_split, split_axis)
            split_overlay_size += a
            split_cont_sum_size += b
            split_points_sum += d
            num_overlays += c
            split_intersections += e
            split_overlay_fill += f
        return split_overlay_size, split_cont_sum_size, num_overlays, split_points_sum/len(split), \
            split_intersections/len(split), split_overlay_fill/len(split)


EXMAPLE = [
    (1, 1),
    (2, 2),
    (3, 3),
]
EXMAPLE_2 = [
    (1, 1),
    (3, 1),
    (2, 3),
    (3, 4),
    (4, 3),
]

BED_FOLDER = "/work/project/ladsie_012/ABS.2.2/2021-10-26_NS502-NS521_ABS_CR_RADICL_inputMicroC/bed_files"
BED_SUFFIX = "RNA.sorted.bed_K1K2.bed_K4.bed_R_D.bed_R_D_K1K2.bed_R_D_PRE1.bed"
FILES = [
    BED_FOLDER + "/NS504_P10_Total_3." + BED_SUFFIX,
    BED_FOLDER + "/NS505_N50_Total_1." + BED_SUFFIX,
    BED_FOLDER + "/NS508_P10_NPM_1." + BED_SUFFIX,
    BED_FOLDER + "/NS511_N50_NPM_1." + BED_SUFFIX,
]

if __name__ == "__main__":
    for x in range(3, 100):
        print(10**x, "points:")
        points_up, points_right, have_all = input_to_points(
            "../out/Lister427.sizes", FILES, limit=10**x)
        a, b, c, d, e, f = compute_mem(points_up + points_right)
        print("\ttotal of", a, "overlay integers = ", 100*a/(a+b), "%")
        print("\ttotal of", b, "cont sum integers = ", 100*b/(a+b), "%")
        print("\ttotal of", c, "overlays")
        print("\taverage of", d, "points per overlay")
        print("\taverage of", e, "intersections per overlay")
        print("\taverage overlay fill: ", 100*f, "%")
        print("\ttotal of", 8 * (a+b) / 10**9, "gb")
        if have_all:
            break

 # 2^(xd)(dN+1)+x^d(N/2^(xd)-1)^d
