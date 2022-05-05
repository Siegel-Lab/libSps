# libSps - Sparse Prefix Sums C++ Library

libSps is a C++ Library that allows counting points in n-dimensional space in constant time O(1). 
It was written to analyze [Hi-C and RADICL-seq data](https://en.wikipedia.org/wiki/Chromosome_conformation_capture "Wikipedia") 
and is based on the algorithms in Shekelyan et al.<sup>1</sup> as well as Schmidt et al.<sup>2</sup>.

The basic principle is to compute and store the prefix sums of all datapoints.
For example: The points 1, 3, 3, & 5

            X       XX      X
        |   |   |   |   |   |
        0   1   2   3   4   5

would have the following prefix sums (the prefix sum is incremented at each point): 0, 1, 1, 3, 3, 4

    4 -                     ----
    3 -             --------
    2 -
    1 -     --------
    0 - ----
        |   |   |   |   |   |
        0   1   2   3   4   5

If we now want to count the number of points in any given interval (a, b], 
we merely have to subtract the prefix sum at position a from the prefix sum at position b.
E.g. count( (1, 4] ) = 3 - 1 = 2; i.e. the two points at postion 3 but no other point are within the interval (1, 4].
No matter the number of points in our index nor the size of the queried interval, we need two lookups and one substraction operation to count the points.

## Getting started

The easiest way to install libSps is Bioconda

    conda -c bioconda install libSps

First, we create an index and fill it with points.
The following is example python code:

    from libSps import make_sps_index

    # create a 2D index
    index = make_sps_index()

    # add points
    initial_size = len(index)
    index.add_point((5, 7), "description of point1")
    index.add_point((20, 1), "description of point2")

    # preprocess the index
    dataset_id = index.generate(initial_size, len(index))

Once we have a preprocessed dataset, we can query it:

    a = index.count(dataset_id, (0, 0), (0, 0))
    # a == 0

    b = index.count(dataset_id, (3, 0), (10, 10))
    # b == 1

Note that, index.count() takes constant O(1)! time, no matter the amount of points in the index or size of the area that is queried. 


## Installation

There are two ways to install libSps: via bioconda and via github.
The bioconda installation is easier, but restricts you to using <=5 dimensions.

### BioConda

    conda -c bioconda install libSps

### CMake

To compile libSps from the githug source use the following commands:

    git clone https://github.com/MarkusRainerSchmidt/libSps
    mkdir build
    cd build
    cmake ../libSps/ # see below to config correctly
    make

Various parameters of libSps are set during compile time as they affect the underlying memory layout and as this improves runtime performance. 
Hence, With CMake, you can configure libSps in various ways:

add `-DWITH_PYTHON=ON` to create a shared object file that can be imported as a python module. If this parameter is set:
- set `-DNUM_DIMENSIONS_A=X -DNUM_DIMENSIONS_B=Y ... -DNUM_DIMENSIONS_D=Z` to create indices with X, Y, ..., Z dimensions.
- set `-DW_DEPENDENT_DIM=ON` to create indices where the overlay grids' dimension 1 is dependent on its dimension 0.
- set `-DWO_DEPENDENT_DIM=ON` to create indices where the overlay grids' dimensions are independent of each other.
- set `-DW_CUBES=ON` to create indices where entries are not point-like but interval-like on the first 3 dimensions.
- set `-DW_RECTANGLES=ON` to create indices where entries are not point-like but interval-like on the first 2 dimensions.
- set `-DW_INTERVALS=ON` to create indices where entries are not point-like but interval-like on the first dimension.
- set `-DW_POINTS=ON` to create indices where entries are fully point-like.
- set `-DDISK=ON` to create indices that load all data from a file on startup and store it back to the file on shutdown. Expect them to consume as much RAM as the filesize.
- set `-DCACHED=ON` to create indices that use a cache to load data from and store data to a file dynamically as needed during runtime. Expect this storage type to be slightly slower than the other two options. For large datasets this storage is necessary, as it allows the RAM usage to be independent of the amount of data stored.
- set `-DRAM=ON` to create indices stores all information in RAM and never interacts with the filesystem.

If multiple contradicting options are turned on as e.g. `-DNUM_DIMENSIONS_A=7 -DNUM_DIMENSIONS_B=3` and `-DRAM=ON -DDISK=ON` indices with all valid combinations will be compiled and made available via the `make_sps_index()` factory function.
Hence, you may want to consider turning options off to save compiletime: E.g. `-DCACHED=OFF -DNUM_DIMENSIONS_A=0`.

If you plan using libSps from C++ you can set all these options via template parameters (see "Creating indices" below).

## Usage

The libSps workflow is seperated in two phases: Creating an index and querying the index. 

### Creating indices

#### Python

#### C++

### Querying indices


### Documentation

[Technical Documentation for the python module can be found here.](https://github.com/MarkusRainerSchmidt/libSps/docs/py/index.html "Documentation")


[Technical Documentation for the c++ code can be found here.](https://github.com/MarkusRainerSchmidt/libSps/docs/cpp/index.html "Documentation") 

#### Dependent Dimensions

#### Intervals / Cubes / Rectangles / Orthotopes

## Limitations

- libSps does not work with negative coordinates.
- dimension 1 is the only dimension that can be made dependent on another dimension and it can only be made dependent on dimenison 0.

## Thanks

libSps uses several other projects. These are:
- [Pybind11](https://github.com/pybind/pybind11 "GitHub"): Seamless operability between C++11 and Python
- [STXXL](https://github.com/stxxl/stxxl "GitHub"): Standard Template Library for Extra Large Data Sets

## Citing libSps

For citing libSps, please use:
@todo

## References

1. Shekelyan, M., Dignös, A. & Gamper, J. Sparse prefix sums: Constant-time range sum queries over sparse multidimensional data cubes. Information Systems 82, 136–147 (2019).
2. Schmidt et al. @todo
