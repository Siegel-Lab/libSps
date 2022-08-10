from libSps import VERSION, make_sps_index
import argparse
from bokeh.models import Legend, Span, Label, BoxAnnotation
from bokeh.plotting import figure, save, reset_output
from bokeh.layouts import gridplot
from bokeh.models import Div

parser = argparse.ArgumentParser(description='Create indices for the smoother Hi-C data viewer.')

parser.add_argument('filepath_prefix', help="Prefix path of the index on the filesystem")
parser.add_argument('from_points', help="Index of first point that shall be included in the dataset, defaults to %(default)s", default=0)
parser.add_argument('to_points', help="Index of first point that shall not be included in the dataset anymore. Negative values pick all points. Defaults to %(default)s",
                    default=-1)
parser.add_argument('-d','--num_dimension', metavar='D', 
                    help="Number of dimensions for datapoints in the index, defaults to %(default)s.", default=2)
parser.add_argument('-u','--with_uniform_overlay_boxes', 
                    help="Whether the overlay grid is completely uniform.", 
                    action='store_true')
parser.add_argument('-e','--with_dependent_dimension', 
                    help="Whether the overlay grid in dimension 1 is dependent on the grid in dimension 0.", 
                    action='store_true')
parser.add_argument('-o','--num_orthotope_dimensions', metavar='O', 
                    help="Number of orthotope dimensions (set this to 1, 2, 3, … to add intervals, rectangles, cubes, … instead of points), defaults to %(default)s.", 
                    default=0)
parser.add_argument('-s','--storage_type', metavar='S',
                    help="the way the datastructure interacts with the filesystem, defaults to %(default)s.", 
                    default="Ram")
parser.add_argument('-f','--exp_from', metavar='F',
                    help="Start exponent for the plot's x values, defaults to %(default)s.", 
                    default=0)
parser.add_argument('-t','--exp_to', metavar='T',
                    help="End exponent for the plot's x values, defaults to %(default)s.", 
                    default=10)
parser.add_argument('-b','--exp_base', metavar='B',
                    help="Base value for the plot's x values, defaults to %(default)s.", 
                    default=2)
parser.add_argument('--exp_sub', metavar='SUB',
                    help="Subtrahend for the plot's x values, defaults to %(default)s.", 
                    default=1)
parser.add_argument('-a','--actual', metavar='A',
                    help="x positions for which to compute the actual size.", nargs='*', type=int, default=[])
parser.add_argument('-m','--max_size', metavar='M',
                    help="If the prediction predicts the datastructure to be smaller than M, compute the datastructures actual size, defaults to %(default)s.", 
                    default=0)


args = parser.parse_args()


print("LibSps Version:", VERSION)

FILES = [".prefix_sums", ".coords", ".overlays", ".datsets"]
# load index
index = make_sps_index(args.filepath_prefix, args.num_dimension, args.with_uniform_overlay_boxes, 
                       args.with_dependent_dimension, args.num_orthotope_dimensions, args.storage_type, False )
if args.to_points < 0:
    args.to_points = len(index)

ipsas,opsas,iscas,oscas,ipsps,opsps,iscps,oscps,fsa,fsp,las,lps,eps = ([], [], [], [], [], [], 
                                                                       [], [], [], [], [], [], [])
xs = []
actual_left = float('inf')
actual_right = 0
for fac in sorted([args.exp_base ** _fac - args.exp_sub for _fac in range(args.exp_from, args.exp_to)] + args.actual):
    file_size_estimate = index.estimate_size(fac, args.from_points, args.to_points) / 10**9 # in GB
    fsp.append(file_size_estimate)

    if fac in args.actual or file_size_estimate < args.max_size:
        id = index.generate(0, len(index), fac, verbosity=0)
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

        actual_left = min(actual_left, num_overlays)
        actual_right = max(actual_right, num_overlays)
        file_size = sum(os.path.getsize("test/benchmark_index" + f) for f in FILES)  / 10**9 # in GB
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

    opsp, ipsp, oscp, iscp, num_overlays, lp = index.estimate_num_elements(fac, args.from_points, args.to_points)
    opsps.append(opsp)
    ipsps.append(ipsp)
    oscps.append(oscp)
    iscps.append(iscp)
    lps.append(lp)
    xs.append(num_overlays)