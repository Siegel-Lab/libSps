from build_dbg.libSps import *
#from build_rel.libSps import *
import random
from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool

print_all = False

def combinations(a, b):
    assert len(a) == len(b)
    if len(a) == 0:
        yield []
    else:
        for x in combinations(a[:-1], b[:-1]):
            yield x + [a[-1]]
            yield x + [b[-1]]


def render_overlays(tree, x, title):
    do = {
        "left": [],
        "right": [],
        "bottom": [],
        "top": [],
        "grid_pos": [],
        "index": [],
        "pred_indices": []
    }
    dp = {
        "x": [],
        "y": [],
        "index": [],
    }
    def m(x):
        return min(x, 100)

    for overlay_info in tree.get_overlay_info(x):
        do["left"].append(m(overlay_info.bottom_left[0]))
        do["bottom"].append(m(overlay_info.bottom_left[1]))
        do["right"].append(m(overlay_info.top_right[0]))
        do["top"].append(m(overlay_info.top_right[1]))
        do["grid_pos"].append(str(overlay_info.grid_pos))
        do["index"].append(str(overlay_info.index))
        do["pred_indices"].append(str(overlay_info.pred_indices))
        for p in overlay_info.points:
            dp["x"].append(p[0])
            dp["y"].append(p[1])
            dp["index"].append(str(overlay_info.index))
    output_file(filename="overlays.html", title=title)
    f = figure(sizing_mode="stretch_both", title=title)
    f.quad(source=ColumnDataSource(data=do), line_color="black", fill_color="white", fill_alpha=0.5)
    f.x(source=ColumnDataSource(data=dp), color="blue", size=10)
    tooltips = [
        ("index", "@index"),
        ("grid_pos", "@grid_pos"),
        ("pred_indices", "@pred_indices"),
    ]
    hover = HoverTool(tooltips=tooltips, point_policy="follow_mouse")
    f.add_tools(hover)
    save(f)



def fixed(tree, points, d=2, cont=0, area=False, enforce_wide_queries=True):
    for idx, pos in enumerate(points):
        if area:
            tree.add_point(*pos, "p" + str(idx))
        else:
            tree.add_point(pos, "p" + str(idx))
    x = tree.generate(0, len(tree), 5 if print_all else 1)
    if print_all:
        print("done generating")
        print(tree)
        print("generated")
    if area:
        p_start = points + [(list(range(d)), [max(p[i] for _, p in points)+1 for i in range(d)])]
        p_end = p_start + [(list(range(d)), [max(p[i] for _, p in p_start)+1 for i in range(d)])]
        max_w = [max(p2[i] - p1[i] for p1, p2 in points) for i in range(d)]
    else:
        p_start = points + [[max(p[i] for p in points)+1 for i in range(d)]]
        p_end = p_start + [[max(p[i] for p in p_start)+1 for i in range(d)]]
    for a1 in p_start:
        for a2 in p_end:
            if area:
                p1 = a1[0]
                p2 = a2[1]
                if enforce_wide_queries:
                    if not all(j-i >= w for w, i, j in zip(max_w, p1, p2)):
                        if print_all:
                            print("not wide enough")
                        continue
            else:
                p1 = a1
                p2 = a2
            if all(a < b for a, b in zip(p1, p2)):
                cnt = tree.count(x, p1, p2, 5 if print_all else 0)
                if area:
                    truth = sum(1 if all(i >= k and j < l for i, j, k, l in zip(ps, pe, p1, p2)) else 0 \
                                    for ps, pe in points)
                else:
                    truth = sum(1 if all(i >= j and i < k for i, j, k in zip(p, p1, p2)) else 0 for p in points)
                if not cnt == truth:
                    render_overlays(tree, x, str(d) + "-" + str(cont))
                    print("counts:", cnt, truth)
                    for pc in combinations(p1, p2):
                        if area:
                            if pc == p2:
                                corner_c = sum(1 if all(i < j for i, j in zip(pe, pc)) else 0 for _, pe in points)
                            else:
                                corner_c = sum(1 if all(i < j for i, j in zip(ps, pc)) else 0 for ps, _ in points)
                        else:
                            corner_c = sum(1 if all(i < j for i, j in zip(p, pc)) else 0 for p in points)
                        print("expected corner count", corner_c, "for", pc)
                    print("query", p1, p2)
                    print("points", points)
                    print(tree)
                    print("failure", d, cont)
                    exit()
                else:
                    if print_all:
                        print("OK")
            else:
                if print_all:
                    print("not valid")
    print("success", d, cont)


def test(tree, d, n=30, area=False, enforce_wide_queries=True):
    cont = 0
    for x in range(1, n):
        for _ in range(min(x*2, 100)):
            tree.clear()
            points = []
            for _ in range(x):
                pos_s = []
                pos_e = []
                for _ in range(d):
                    pos_s.append(random.choice(range(x)))
                    pos_e.append(pos_s[-1])
                    if d < 2:
                        pos_e[-1] += random.choice(range(x))
                if area:
                    points.append((pos_s, pos_e))
                else:
                    points.append(pos_s)
                if print_all:
                    print("adding", points[-1])
            fixed(tree, points, d, cont, area, enforce_wide_queries)
            cont += 1



random.seed(6846854546132)
#fixed(DependantDimSparsePrefixSum_2D("test/blub2"), 2, [[0,1], [1,0], [1,2], [0,3], [1,4]])

#test(CachedDependantDimPrefixSum_2D("test/blub1", True), 2)
#test(CachedDependantDimPrefixSum_4D("test/blub3", True), 4)


#test(DiskDependantDimPointsPrefixSum_3D("test/blub2", True), 3)
#test(DiskDependantDimPointsPrefixSum_5D("test/blub4", True), 5)

#test(DiskDependantDimPointsPrefixSum_2D("test/blub5", True), 2)

#test(DiskDependantDimRectanglesPrefixSum_2D("test/blub6", True), 2, area=True, enforce_wide_queries=False)
test(DiskDependantDimRectanglesPrefixSum_3D("test/blub8", True), 3, area=True, enforce_wide_queries=False)


#test(DiskDependantDimRectanglesPrefixSum_3D("test/blub7", True), 3, area=True)

## libSps.index("", 2, True, 0, "Disk", True)