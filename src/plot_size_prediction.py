import sys,os
import argparse
from bokeh.models import Legend, Span, Label, BoxAnnotation
from bokeh.plotting import figure, save, reset_output
from bokeh.layouts import gridplot
from bokeh.models import Div
from bokeh.io import output_file

parser = argparse.ArgumentParser(description='Create indices for the smoother Hi-C data viewer.')

parser.add_argument('filepath_prefix', help="Prefix path of the index on the filesystem")
parser.add_argument('-f', '--from_points', help="Index of first point that shall be included in the dataset, defaults to %(default)s", default=0, type=int)
parser.add_argument('-t', '--to_points', help="Index of first point that shall not be included in the dataset anymore. Negative values pick all points. Defaults to %(default)s",
                    default=-1, type=int)
parser.add_argument('-d','--num_dimension', metavar='D', 
                    help="Number of dimensions for datapoints in the index, defaults to %(default)s.", default=2,
                    type=int)
parser.add_argument('-u','--with_uniform_overlay_boxes', 
                    help="Whether the overlay grid is completely uniform.", 
                    action='store_true', default=False)
parser.add_argument('-e','--with_dependent_dimension', 
                    help="Whether the overlay grid in dimension 1 is dependent on the grid in dimension 0.", 
                    action='store_true', default=False)
parser.add_argument('-x','--factor_as_x_axis', 
                    help="Display the factor as x instead of the number of overlays.", 
                    action='store_true', default=False)
parser.add_argument('-o','--num_orthotope_dimensions', metavar='O', 
                    help="Number of orthotope dimensions (set this to 1, 2, 3, … to add intervals, rectangles, cubes, … instead of points), defaults to %(default)s.", 
                    default=0, type=int)
parser.add_argument('-s','--storage_type', metavar='S',
                    help="the way the datastructure interacts with the filesystem, defaults to %(default)s.", 
                    default="Disk")
parser.add_argument('-F','--exp_from', metavar='F',
                    help="Start exponent for the plot's x values, defaults to %(default)s.", 
                    default=0, type=int)
parser.add_argument('-T','--exp_to', metavar='T',
                    help="End exponent for the plot's x values, defaults to %(default)s.", 
                    default=10, type=int)
parser.add_argument('-b','--exp_base', metavar='B',
                    help="Base value for the plot's x values, defaults to %(default)s.", 
                    default=2, type=float)
parser.add_argument('-a','--actual', metavar='A',
                    help="x positions for which to compute the actual size.", nargs='*', type=int, default=[])
parser.add_argument('-m','--max_size', metavar='M',
                    help="If the prediction predicts the datastructure to be smaller than M, compute the datastructures actual size, defaults to %(default)s.", 
                    default=0, type=float)
parser.add_argument('-q','--quality_overlays', metavar='M',
                    help="If the prediction predicts the datastructure to be smaller than M, compute the datastructures actual size, defaults to %(default)s.", 
                    default=10000, type=int)
parser.add_argument('-Q','--quality_points', metavar='M',
                    help="If the prediction predicts the datastructure to be smaller than M, compute the datastructures actual size, defaults to %(default)s.", 
                    default=1000, type=int)
parser.add_argument('-l','--lib_sps_path', metavar='P',
                    help="Path to libSps' folder in case it is located outside of the current working directory, defaults to %(default)s.", default=".")
parser.add_argument('-O','--skip_optimum', 
                    help="Skip computing and displaying the automatically determined optimum", 
                    action='store_true', default=False)


args = parser.parse_args()

sys.path.append(args.lib_sps_path)
try:
    from libSps import VERSION, make_sps_index
except:
    print("Could not load libSps.")
    print("Please provide a path to libSps via the -l <path/to/folder> parameter.")


for name, param in [
            ("LibSps Version", VERSION),
            ("prefix", args.filepath_prefix),
            ("no. dimensions", args.num_dimension),
            ("w/ dependent dimension", args.with_dependent_dimension),
            ("w/ uniform overlay boxes", args.with_uniform_overlay_boxes),
            ("no. orthotope dimensions", args.num_orthotope_dimensions),
            ("storage type", args.storage_type),
        ]:
    print(name + ":", param, sep=" "*(30-len(name)))

FILES = [".prefix_sums", ".coords", ".overlays", ".datsets"]
# load index
index = make_sps_index(args.filepath_prefix, args.num_dimension, args.with_dependent_dimension, 
                       args.with_uniform_overlay_boxes, args.num_orthotope_dimensions, args.storage_type, False )
if args.to_points < 0:
    args.to_points = len(index)
fac_list = sorted([args.exp_base ** _fac for _fac in range(args.exp_from, args.exp_to)] + args.actual)

print("predicting index size for the points", args.from_points, "-", args.to_points,
      "and factors", *fac_list, ".")

ipsas,opsas,iscas,oscas,ipsps,opsps,iscps,oscps,fsa,fsp,las,lps,eps = ([], [], [], [], [], [], 
                                                                       [], [], [], [], [], [], [])
xs = []
actual_left = float('inf')
actual_right = 0

print("predicting sizes...")
size_estimates = index.estimate_num_elements(fac_list, args.from_points, args.to_points, args.quality_overlays, 
                                             args.quality_points)

print("checking sizes", end="", flush=True)
for fac, (opsp, ipsp, oscp, iscp, num_overlays, lp, file_size_estimate) in zip(fac_list, size_estimates):
    file_size_estimate = file_size_estimate / 10**9 # in GB
    opsps.append(opsp)
    ipsps.append(ipsp)
    oscps.append(oscp)
    iscps.append(iscp)
    lps.append(lp)
    if args.factor_as_x_axis:
        xs.append(fac)
    else:
        xs.append(num_overlays)

    fsp.append(file_size_estimate)

    if fac in args.actual or file_size_estimate < args.max_size:
        id = index.generate(args.from_points, args.to_points, fac, 0, args.quality_overlays, args.quality_points)
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

        if args.factor_as_x_axis:
            actual_left = min(actual_left, fac)
            actual_right = max(actual_right, fac)
        else:
            actual_left = min(actual_left, num_overlays)
            actual_right = max(actual_right, num_overlays)
        file_size = index.get_size(id)/ 10**9 # in GB
        fsa.append(file_size)
        eps.append(file_size_estimate/file_size)

        index.pop_dataset()
    else:
        ipsas.append(float('NaN'))
        opsas.append(float('NaN'))
        iscas.append(float('NaN'))
        oscas.append(float('NaN'))
        las.append(float('NaN'))
        fsa.append(float('NaN'))
        eps.append(float('NaN'))
    
    print(".", end="", flush=True)
print()

if not args.skip_optimum:
    print("searching for optimum...")

    picked_num = index.pick_num_overlays(args.from_points, args.to_points)
    picked_size = index.estimate_num_elements([index.to_factor(picked_num, args.from_points, args.to_points)], 
                                    args.from_points, args.to_points)[0][-1] / 10**9 # in GB
    if args.factor_as_x_axis:
        picked_num = index.to_factor(picked_num, args.from_points, args.to_points)

print("generating and saving plot...")


axis_type = "log"

plot_ips = figure(x_axis_type=axis_type, title="internal grid")
plot_ips.circle(x=xs, y=ipsas, legend_label="actual", fill_color="green", 
                line_color=None, size=8)
plot_ips.circle(x=xs, y=ipsps, legend_label="predicted", fill_color=None, 
                line_color="red", size=8)


plot_ops = figure(x_axis_type=axis_type, title="overlay grid")
plot_ops.circle(x=xs, y=opsas, fill_color="green", line_color=None, size=8)
plot_ops.circle(x=xs, y=opsps, fill_color=None, line_color="red", size=8)

plot_isc = figure(x_axis_type=axis_type, title="internal lookup table")
plot_isc.circle(x=xs, y=iscas, fill_color="green", line_color=None, size=8)
plot_isc.circle(x=xs, y=iscps, fill_color=None, line_color="red", size=8)

plot_osc = figure(x_axis_type=axis_type, title="overlay lookup table")
plot_osc.circle(x=xs, y=oscas, fill_color="green", line_color=None, size=8)
plot_osc.circle(x=xs, y=oscps, fill_color=None, line_color="red", size=8)

plot_ls = figure(x_axis_type=axis_type, title="global lookup table")
plot_ls.circle(x=xs, y=las, fill_color="green", line_color=None, size=8)
plot_ls.circle(x=xs, y=lps, fill_color=None, line_color="red", size=8)

plot_fs = figure(x_axis_type=axis_type, y_axis_type="log", title="total file size")
act = plot_fs.circle(x=xs, y=fsa, fill_color="green", line_color=None, size=8)
pred = plot_fs.circle(x=xs, y=fsp, fill_color=None, line_color="red", size=8)

if not args.skip_optimum:
    plot_fs.add_layout(BoxAnnotation(right=picked_num, top=picked_size,
                                    line_color='black', fill_color=None,
                                    line_dash='solid', line_width=1, line_alpha=1))
    plot_fs.add_layout(Label(x=picked_num + 1, y=0, y_units='screen', text="predicted minimum",
                            background_fill_color="white", background_fill_alpha=1.0,
                            level="overlay"))

plot_ep = figure(x_axis_type=axis_type, y_axis_type="log", title="total file size error")
plot_ep.line(x=xs, y=eps, line_color="black")
plot_ep.circle(x=xs, y=eps, fill_color="black", line_color=None, size=8)
plot_ep.yaxis[0].ticker.base = 10
plot_ep.yaxis[0].formatter.ticker = plot_ep.yaxis[0].ticker

for p in [plot_ips, plot_ops, plot_isc, plot_osc, plot_ls, plot_fs, plot_ep]:
    p.yaxis.axis_label = "amount"
    if args.factor_as_x_axis:
        p.xaxis.axis_label = "factor"
    else:
        p.xaxis.axis_label = "No. overlay boxes"
    p.background_fill_color = "whitesmoke"
    if actual_right > actual_left:
        actual_box = BoxAnnotation(left=actual_left, right=actual_right, fill_color='white',
                                   level='image', fill_alpha=1)
        p.add_layout(actual_box)
plot_fs.yaxis.axis_label = "size [GB]"
plot_ep.yaxis.axis_label = "Log10(predicted filesize / actual filesize)"

output_file(args.filepath_prefix + ".size_prediction.html", title="Size prediction for " + args.filepath_prefix)

save(gridplot([[plot_ips, plot_ops, plot_isc, plot_osc, plot_ls], [plot_fs, plot_ep, None, None, None]], 
            sizing_mode="scale_width", merge_tools=False))

print("saved plot at", args.filepath_prefix + ".size_prediction.html")
print("closing index...")