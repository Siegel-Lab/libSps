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



def fixed(tree, points, d=2, cont=0):
    for idx, pos in enumerate(points):
        tree.add_point(pos, "p" + str(idx))
    x = tree.generate(0, len(points))
    if print_all:
        print("done generating")
        print(tree)
        print("generated")
    for p1 in points:
        for p2 in points + [[max(p[i] for p in points)+1 for i in range(d)]]:
            if all(a < b for a, b in zip(p1, p2)):
                cnt = tree.count(x, p1, p2)
                truth = sum(1 if all(i >= j and i < k for i, j, k in zip(p, p1, p2)) else 0 for p in points)
                itr = list(range(truth))#tree.get(x, p1, p2)
                if not cnt == len(itr) == truth:
                    render_overlays(tree, x, str(d) + "-" + str(cont))
                    print("counts:", cnt, len(itr), truth)
                    for pc in combinations(p1, p2):
                        corner_c = sum(1 if all(i < j for i, j in zip(p, pc)) else 0 for p in points)
                        print("expected corner count", corner_c, "for", pc)
                    print("query", p1, p2)
                    print("points", points)
                    print("iteration result", itr)
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


def test(tree, d, n=30):
    cont = 0
    for x in range(1, n):
        for _ in range(min(x*2, 100)):
            tree.clear()
            points = []
            for _ in range(x):
                pos = []
                for _ in range(d):
                    pos.append(random.choice(range(x)))
                points.append(pos)
                if print_all:
                    print("adding", points[-1])
            fixed(tree, points, d, cont)
            cont += 1



random.seed(6846854546132)
#fixed(DependantDimSparsePrefixSum_2D("test/blub2"), 2, [[0,1], [1,0], [1,2], [0,3], [1,4]])

#test(CachedDependantDimPrefixSum_2D("test/blub1", True), 2)
#test(CachedDependantDimPrefixSum_3D("test/blub2", True), 3)
#test(CachedDependantDimPrefixSum_4D("test/blub3", True), 4)
test(CachedDependantDimPrefixSum_5D("test/blub4", True), 5)
