import os
import random

def parse_heatmap(in_filename):
    with open(in_filename, "r") as in_file_1:
        for line in in_file_1:
            # parse file columns
            read_name, strnd_1, chr_1, pos_1, _, strnd_2, chr_2, pos_2, _2, mapq_1, mapq_2 = line.split()
            # convert number values to ints
            pos_1, pos_2, mapq_1, mapq_2 = (
                int(x) for x in (pos_1, pos_2, mapq_1, mapq_2))
            pos_1 -= 1
            pos_2 -= 1

            yield chr_1, pos_1, chr_2, pos_2, max(mapq_1, mapq_2)

def input_to_points(len_file_name, heatmap_filenames, limit=float('inf')):
    curr_start = 0
    chr_start_pos = {}
    with open(len_file_name, "r") as len_file:
        for line in len_file:
            chr_name, chr_len = line.split("\t")
            chr_start_pos[chr_name] = curr_start
            curr_start += int(chr_len)
    points_up = []
    points_right = []
    cnt = 0
    for file_name in heatmap_filenames:
        for chr_1, pos_1, chr_2, pos_2, map_q in parse_heatmap(file_name):
            x = chr_start_pos[chr_1] + int(pos_1)
            y = chr_start_pos[chr_2] + int(pos_2)
            max_dist_betw_alignments = random.choice(list(range(10)))
            annotations = [True] + [random.choice([True, False]) for _ in range(10)]
            p = points_up if y > x else points_right
            for a, _ in enumerate(annotations):
                p.append((map_q, max_dist_betw_alignments, a, x, y))
            cnt += 1
            if cnt > limit:
                return points_up, points_right
    return points_up, points_right

def sort_file(file_name, c):
    # sort numerically, and by the column c
    os.system("sort -n -k " + str(c + 1) + " " + file_name + " -o " + file_name)

def write_bins_to_other_file(bins, num, c):
    ret = []
    last = None
    for bin in bins:
        x = bin[c]
        if x != last:
            num -= 1
            last = x
        if num < 0:
            return ret
        ret.append(bin)
    return ret


def compute_num_bins_helper(bins, dimensions, top=True):
    cnt = 0
    d = dimensions[0]
    bins.sort(key=lambda x: x[d])
    if len(dimensions) == 1:
        last = None
        for bin in bins:
            x = bin[d]
            if x != last:
                cnt += 1
                last = x
        return cnt, [cnt]
    else:
        num = 1
        summs = []
        while True:
            sub_bins = write_bins_to_other_file(bins, num, d)
            if top:
                print(len(sub_bins), "/", len(bins))
            _cnt, _summs = compute_num_bins_helper(sub_bins, dimensions[1:], False)
            cnt += _cnt
            summs.append(_summs)
            if len(sub_bins) == len(bins):
                break
            num += 1
        xx = []
        for idx in range(len(dimensions)-1):
            xx.append(sum(x[idx] for x in summs) / len(summs))
        return cnt, [num] + xx


def compute_num_bins(bins, dimensions):
    a, b = compute_num_bins_helper(bins, dimensions)
    print("bins total", a)
    for idx, x in enumerate(b):
        print("dimension", idx, "average number bins", x)

EXMAPLE = [
    (1,1),
    (2,2),
    (3,3),
]
EXMAPLE_2 = [
    (1,1),
    (3,1),
    (2,3),
    (3,4),
    (4,3),
]

BED_FOLDER="/work/project/ladsie_012/ABS.2.2/2021-10-26_NS502-NS521_ABS_CR_RADICL_inputMicroC/bed_files"
BED_SUFFIX="RNA.sorted.bed_K1K2.bed_K4.bed_R_D.bed_R_D_K1K2.bed_R_D_PRE1.bed"
FILES = [
    BED_FOLDER + "/NS504_P10_Total_3." + BED_SUFFIX,
    BED_FOLDER + "/NS505_N50_Total_1." + BED_SUFFIX,
    BED_FOLDER + "/NS508_P10_NPM_1." + BED_SUFFIX,
    BED_FOLDER + "/NS511_N50_NPM_1." + BED_SUFFIX,
]

if __name__ == "__main__":
    #compute_num_bins(EXMAPLE_2, [1, 0])
    #exit()

    points_up, points_right = input_to_points("../out/Lister427.sizes", FILES)
    compute_num_bins(points_up, [0, 1, 2, 3, 4])
    compute_num_bins(points_right, [0, 1, 2, 4, 3])