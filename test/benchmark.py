import sys
import os
sys.path.append(os.getcwd())
from build_benchmark.libSps import VERSION, make_sps_index, ValDiskSimpleVector_2D#, ValCachedSimpleVector_2D
import random
import time
import os
import time
from bokeh.models import Legend, Span, Label, BoxAnnotation
from bokeh.plotting import figure, save, reset_output
from bokeh.layouts import gridplot
from bokeh.models import Div

MIN_FILL = 3# 1
MAX_FILL = 4# 8
N_QUERY = 10000 #100k
#N_QUERY = 10000000 #10,000k
QUAL_OVERLAY = 1000
QUAL_PTS = 1000

files = [".prefix_sums", ".coords", ".overlays", ".datsets"]
files_full = files + [".desc", ".points"]

def query(index, id, dims, genome_size, n=N_QUERY):
    bins = []
    for _ in range(n):
        pos_s = []
        pos_e = []
        for _ in range(dims):
            pos_s.append(random.randrange(genome_size))
            pos_e.append(pos_s[-1] + random.randrange(100))
        bins.append((pos_s, pos_e))
    t1 = time.perf_counter()
    index.count_multiple(id, bins)
    t2 = time.perf_counter()
    #querytime
    #print((t2-t1), n, sep="\t", end="\t")

def enumerate_all_combinations(a, b):
    if len(a) > 0:
        for c in enumerate_all_combinations(a[1:], b[1:]):
            yield [a[0]] + c
            yield [b[0]] + c
    else:
        yield []

def fill(n, index, name, dims, is_ort, density, distrib, asp_ratio, offset):
    index.clear()
    t1 = time.perf_counter()
    pts = []
    for _ in range(n):
        if distrib == "log_norm":
            x = int(random.lognormvariate(0, 1)* density / asp_ratio) + offset
        else:
            x = random.randrange(int(n * density / asp_ratio)) + offset

        if distrib == "diagonal":
            y = x
        elif distrib == "dist_dep_dec":
            y = min(max(0, int(abs(random.normalvariate(0, 10))) + x), int(n * density * asp_ratio))
        elif distrib == "log_norm":
            y = int(random.lognormvariate(0, 1) * density * asp_ratio) + offset
        else:
            y = random.randrange(int(n * density * asp_ratio)) + offset
        pos_s = [x, y]
        for _ in range(2, dims):
            pos_s.append(random.randrange(int(n * density * asp_ratio)))
        if is_ort:
            pos_s = [x//2 for x in pos_s]
            #pos_e = [max(pos_s[0], int(n * density / asp_ratio)) + random.randrange((int(n * density / asp_ratio))//2)]
            pos_e = []
            for s in pos_s: #pos_s[1:]:
                #pos_e.append(max(s, int(n * density * asp_ratio)) + random.randrange((int(n * density * asp_ratio))//2))
                pos_e.append(s)
            index.add_point(pos_s, pos_e, "")
            for pos_x in enumerate_all_combinations(pos_s, pos_e):
                pts.append(pos_x + [b-a for a, b in zip(pos_s, pos_e)])
        else:
            index.add_point(pos_s, "")
            pts.append(pos_s)
            if distrib == "dup":
                for _ in range(2**dims):
                    index.add_point(pos_s, "")
                    pts.append(pos_s)
    t2 = time.perf_counter()
    #id = index.generate(0, len(index), fac, verbosity=0)
    t3 = time.perf_counter()
    # index_name, filltime, generatetime
    #print(*name, n, (t2-t1), (t3-t2), sep="\t", end="\t")
    #return id
    if False:
        for d in range(dims + (dims if is_ort else 0)):
            print()
            print("dim", d)
            histo = {}
            for x in [x[d] for x in pts]:
                if not x in histo:
                    histo[x] = 0
                histo[x] = histo[x] + 1
            for x in range(max(histo.keys())+1):
                if not x in histo:
                    histo[x] = 0
                print(x, "\t", histo[x], "\t", "*"*histo[x])

    ret = 1
    for r in [max(p[d] for p in pts) - min(p[d] for p in pts) for d in range(dims)]:
        ret = ret * r
    return ret

def make_indices():
    indices = []
    index_names = []
    params = []
    for dims in [2]: #[2, 3]:
        for ort_dim in [False, True]:
            num_ort_dims = dims if ort_dim else 0
            for storage in ["Disk"]: #["Disk", "Cached"]:
                for uniform_overlay_grid in [True]: #[False, True]:
                    for dep_dim in ([False] if uniform_overlay_grid else [False, True]):
                        xs = [
                            str(dims) + "d", 
                            ("rect" if ort_dim else "point"), 
                            ("ram" if storage == "Disk" else "cached"), 
                            ("uniform" if uniform_overlay_grid else "adaptive"), 
                            ("dep" if dep_dim else "indep")
                        ]
                        index_names.append(xs)
                        indices.append(("test/benchmark_index", dims, dep_dim, uniform_overlay_grid, 
                                        num_ort_dims, storage))
                        params.append((dims, ort_dim))
    return indices, index_names, params

def test(plot=True, max_pred_file_size=1, fac_base=2):
    print("#VERSION", VERSION)
    print("#N_QUERY", N_QUERY)
    indices, index_names, params = make_indices()
    #print("dims", "ort", "storage", "uniform overlay grid", "dep", "fac", "n", "filltime", "generatetime [s]", "querytime [s]", "num queried bins", "filesize actual [gb]", "filesize predicted [gb]", "filesize difference [gb]", sep="\t")
    
    plots = []
    plots.append([None, 
                  Div(text="<b>internal grid</b>", align="center"), 
                  Div(text="<b>overlay grid</b>", align="center"), 
                  Div(text="<b>internal lookup table</b>", align="center"), 
                  Div(text="<b>overlay lookup table</b>", align="center"), 
                  Div(text="<b>global lookup table</b>", align="center"), 
                  Div(text="<b>total file size</b>", align="center"),
                  Div(text="<b>total file size error</b>", align="center")])
    for _n in range(MIN_FILL, MAX_FILL):
        n = 10**_n
        for offset in [0]: #[0, 100]:
            for density in [1]: #[0.1, 1, 10]:
                for asp_ratio in [1]:# [1, 10]:
                    for distrib in ["even", "dup"]: #["even", "dist_dep_dec", "diagonal", "log_norm", "dup"]:
                        for index_params, name, param in zip(indices, index_names, params):
                            ipsas,opsas,iscas,oscas,ipsps,opsps,iscps,oscps,fsa,fsp,las,lps,eps,gcs = ([], [], [], 
                                    [], [], [], [], [], [], [], [], [], [], [])
                            xs = []
                            actual_left = float('inf')
                            actual_right = 0
                            print(*name, "n=", n, "s=", int(n*density), distrib, "a=", asp_ratio, "o=", offset)
                            for _fac in range(0, 10):
                                fac = fac_base ** _fac - 1
                                index = make_sps_index(*index_params)
                                num_cells = fill(n, index, name if _n == MIN_FILL else [""]*6, *param, 
                                                 density, distrib, asp_ratio, offset)

                                if _fac == 0:
                                    picked_num = index.pick_num_overlays(0, len(index), QUAL_OVERLAY, QUAL_PTS)
                                    picked_size = index.estimate_size( 
                                            index.to_factor(picked_num, 0, len(index)),
                                            0, len(index), QUAL_OVERLAY, QUAL_PTS) / 10**9 # in GB
                                    picked_size_s = index.estimate_size( 
                                            index.to_factor(int(num_cells ** ( 1/(param[0]+param[1])) ), 0, len(index)),
                                            0, len(index), QUAL_OVERLAY, QUAL_PTS) / 10**9 # in GB
                                    picked_size_g = index.estimate_size( 
                                            index.to_factor(int(num_cells ** ( 1/2 )), 0, len(index)),
                                            0, len(index), QUAL_OVERLAY, QUAL_PTS) / 10**9 # in GB

                                file_size_estimate = index.estimate_size(fac, 0, len(index), 
                                                                         QUAL_OVERLAY, QUAL_PTS) / 10**9 # in GB
                                fsp.append(file_size_estimate)

                                if file_size_estimate < max_pred_file_size:
                                    id = index.generate(0, len(index), fac, 0, QUAL_OVERLAY, QUAL_PTS)
                                    query(index, id, param[0], n)
                                    ipsa = index.get_num_internal_prefix_sums(id)
                                    ipsas.append(ipsa)
                                    opsa = index.get_num_overlay_prefix_sums(id)
                                    opsas.append(opsa)
                                    isca = index.get_num_internal_sparse_coords(id)
                                    iscas.append(isca)
                                    osca = index.get_num_overlay_sparse_coords(id)
                                    oscas.append(osca)
                                    la = index.get_num_global_sparse_coords(id)
                                    las.append(la)
                                    gcs.append(index.get_size(id)/ 10**9) # in GB
                                else:
                                    ipsas.append(float('NaN'))
                                    opsas.append(float('NaN'))
                                    iscas.append(float('NaN'))
                                    oscas.append(float('NaN'))
                                    las.append(float('NaN'))
                                    gcs.append(float('NaN'))
                                
                                opsp, ipsp, oscp, iscp, num_overlays, lp = index.estimate_num_elements(
                                                                                    fac, 0, len(index), 
                                                                                    QUAL_OVERLAY, QUAL_PTS)
                                opsps.append(opsp)
                                ipsps.append(ipsp)
                                oscps.append(oscp)
                                iscps.append(iscp)
                                lps.append(lp)
                                xs.append(num_overlays)
                                del index
                                if file_size_estimate < max_pred_file_size:
                                    actual_left = min(actual_left, num_overlays)
                                    actual_right = max(actual_right, num_overlays)
                                    file_size = sum(os.path.getsize("test/benchmark_index" + f) for f in files)  / 10**9 # in GB
                                    fsa.append(file_size)
                                    eps.append(file_size_estimate/file_size)
                                else:
                                    fsa.append(float('NaN'))
                                    eps.append(float('NaN'))
                                
                                #print(file_size / 10**9, file_size_estimate / 10**9, (file_size - file_size_estimate) / 10**9, sep="\t")
                                for f in files_full:
                                    os.remove("test/benchmark_index" + f)
                                if file_size_estimate < max_pred_file_size:
                                    print(".", end="", flush=True)
                                else:
                                    print("!", end="", flush=True)
                            print()
                            if plot:
                                for axis_type in ["log"]: #["linear", "log"]:
                                    plot_ips = figure(x_axis_type=axis_type)
                                    plot_ips.circle(x=xs, y=ipsas, legend_label="actual", fill_color="green", 
                                                    line_color=None, size=8)
                                    plot_ips.circle(x=xs, y=ipsps, legend_label="predicted", fill_color=None, 
                                                    line_color="red", size=8)
                                    

                                    plot_ops = figure(x_axis_type=axis_type)
                                    plot_ops.circle(x=xs, y=opsas, fill_color="green", line_color=None, size=8)
                                    plot_ops.circle(x=xs, y=opsps, fill_color=None, line_color="red", size=8)

                                    plot_isc = figure(x_axis_type=axis_type)
                                    plot_isc.circle(x=xs, y=iscas, fill_color="green", line_color=None, size=8)
                                    plot_isc.circle(x=xs, y=iscps, fill_color=None, line_color="red", size=8)

                                    plot_osc = figure(x_axis_type=axis_type)
                                    plot_osc.circle(x=xs, y=oscas, fill_color="green", line_color=None, size=8)
                                    plot_osc.circle(x=xs, y=oscps, fill_color=None, line_color="red", size=8)

                                    plot_ls = figure(x_axis_type=axis_type)
                                    plot_ls.circle(x=xs, y=las, fill_color="green", line_color=None, size=8)
                                    plot_ls.circle(x=xs, y=lps, fill_color=None, line_color="red", size=8)
                                
                                    plot_fs = figure(x_axis_type=axis_type, y_axis_type="log")
                                    act = plot_fs.circle(x=xs, y=fsa, fill_color="green", line_color=None, size=8)
                                    plot_fs.circle(x=xs, y=gcs, fill_color=None, line_color="blue", size=9)
                                    pred = plot_fs.circle(x=xs, y=fsp, fill_color=None, line_color="red", size=8)
                                    
                                    plot_fs.add_layout(BoxAnnotation(right=picked_num, bottom=picked_size,
                                                                    line_color='black', fill_color=None,
                                                                    line_dash='solid', line_width=1, line_alpha=1))
                                    plot_fs.add_layout(Label(x=picked_num + 1, y=230, y_units='screen', text="libSps",
                                                            background_fill_color="white", background_fill_alpha=1.0,
                                                            level="overlay"))
                                    x_pos = num_cells ** ( 1/(param[0]+param[1]) )
                                    plot_fs.add_layout(BoxAnnotation(right=x_pos, bottom=picked_size_s,
                                                                    line_color='black', fill_color=None,
                                                                    line_dash='dotted', line_width=1, line_alpha=1))
                                    plot_fs.add_layout(Label(x=x_pos + 1, y=210, y_units='screen', text="Shekelyan et al.",
                                                            background_fill_color="white", background_fill_alpha=1.0,
                                                                                    level="overlay"))
                                    x_pos = num_cells ** ( 1/2 )
                                    plot_fs.add_layout(BoxAnnotation(right=x_pos, bottom=picked_size_g,
                                                                    line_color='black', fill_color=None,
                                                                    line_dash='dashed', line_width=1, line_alpha=1))
                                    plot_fs.add_layout(Label(x=x_pos + 1, y=190, y_units='screen', text="Geffner et al.",
                                                            background_fill_color="white", background_fill_alpha=1.0,
                                                                                    level="overlay"))
                                                                                
                                    plot_ep = figure(x_axis_type=axis_type, y_axis_type="log", y_range=(0.09,11))
                                    plot_ep.line(x=xs, y=eps, line_color="black")
                                    plot_ep.circle(x=xs, y=eps, fill_color="black", line_color=None, size=8)
                                    plot_ep.yaxis[0].ticker.base = 10
                                    plot_ep.yaxis[0].formatter.ticker = plot_ep.yaxis[0].ticker

                                    for p in [plot_ips, plot_ops, plot_isc, plot_osc, plot_ls, plot_fs, plot_ep]:
                                        p.yaxis.axis_label = "amount"
                                        p.xaxis.axis_label = "num overlays"
                                        p.background_fill_color = "lightgrey"
                                        if actual_right > actual_left:
                                            actual_box = BoxAnnotation(left=actual_left, right=actual_right, fill_color='white',
                                                                    level='image', fill_alpha=1)
                                            p.add_layout(actual_box)
                                    plot_fs.yaxis.axis_label = "size [GB]"
                                    plot_ep.yaxis.axis_label = "Log10(predicted / actual) filesize"


                                    plots.append([Div(text="<b>" + "<br>".join(name) + "<br>n=" + str(n) + "<br>s=" + 
                                                    str(int(density * n)) + "<br>" + distrib + "<br>a=" + str(asp_ratio) + "<br>o=" + str(offset) + 
                                                    "</b>", 
                                                    width=50, sizing_mode="stretch_height",
                                                    align="center"),
                                                plot_ips, plot_ops, plot_isc, plot_osc, plot_ls, plot_fs, plot_ep])
                                    #plots.append([plot_fs])

    if plot:
        print("saving plots")
        save(gridplot(plots, sizing_mode="scale_width", merge_tools=False))


def test_max_io(n_query=N_QUERY):
    print("#VERSION", VERSION)
    print("#N_QUERY", N_QUERY)
    print("vector", "storage", "n", "querytime [s]", sep="\t")
    for _n in range(MIN_FILL, MAX_FILL):
        vec_d = ValDiskSimpleVector_2D("test/benchmark_index_d", True)
        n = 10**_n
        for _ in range(n):
            vec_d.add(random.choice(range(n)))
        bins = []
        for _ in range(n_query):
            bins.append(random.choice(range(n)))
        t1 = time.perf_counter()
        vec_d.get_multiple(bins)
        t2 = time.perf_counter()
        print("vector" if _n == 1 else "", "ram" if _n == 1 else "", n, (t2-t1), sep="\t")

        del vec_d
        
        os.remove("test/benchmark_index_d.vals")
    for _n in range(MIN_FILL, MAX_FILL):
        vec_c = ValCachedSimpleVector_2D("test/benchmark_index_c", True)
        n = 10**_n
        for _ in range(n):
            vec_c.add(random.choice(range(10000)))
        bins = []
        for _ in range(n_query):
            bins.append(random.choice(range(n)))
        t1 = time.perf_counter()
        vec_c.get_multiple(bins)
        t2 = time.perf_counter()
        print("", "cached" if _n == 1 else "", n, (t2-t1), sep="\t")

        del vec_c
        
        os.remove("test/benchmark_index_c.vals")

random.seed(6846854546132)
#test_max_io()
test()
