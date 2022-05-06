# libSps - The Sparse Prefix Sums C++ Library

## Introduction

libSps is a C++ Library that allows counting points in n-dimensional space in constant time O(1). 
It was written to analyze [Hi-C and RADICL-seq data](https://en.wikipedia.org/wiki/Chromosome_conformation_capture "Wikipedia") 
and is based on the algorithms in Shekelyan et al. [1] as well as Schmidt et al. [2].

The basic principle is to compute and store the prefix sums of all datapoints.
For example: The points 1, 3, 3 & 5

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
    index.add_point((5, 7), "description of point1")
    index.add_point((20, 1), "description of point2")

    # preprocess the index
    dataset_id = index.generate()

Once we have a preprocessed dataset, we can query it:

    a = index.count(dataset_id, (0, 0), (0, 0))
    # a == 0

    b = index.count(dataset_id, (3, 0), (10, 10))
    # b == 1

Note that, index.count() takes constant O(1)! time, no matter the amount of points in the index or size of the area that is queried. 


## Installation

There are two ways to install libSps: via bioconda and via github/cmake.
The bioconda installation is easier, but restricts you to using <=5 dimensions.

### BioConda

    conda -c bioconda install libSps

### Github/CMake

To compile libSps from the githug source use the following commands:

    # clone repository
    git clone https://github.com/MarkusRainerSchmidt/libSps
    cd libSps

    # create and activate conda environment (alternativeley, install the software in conda_env/libSps.yaml yourself)
    ./conda_env/create_env.sh
    conda activate libSps

    # configure build 
    mkdir build
    cd build
    cmake ../libSps/ # Here you can configure various options (see below)

    # build
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

### Python

This section contains the Documentation for using libSps from python3.

#### Creating indices

First you need to pick the correct implementation of the sparse prefix sum index for your use-case.
For this, it is easiest to call the `make_sps_index` factory function with the appropriate parameters.

For example creating an Index for 4-dimensional points, that can deal with a large amount of points saving all data in files and using a cache, use the following parameters:

    index = make_sps_index("example_index_prefix", 4, storage_type="Cached")

If instead of Points, you want to query Intervals, Rectangles, Cubes or Orthotopes of higher dimensionality, adjust the num_orthotope_dimensions parameter.
Note that e.g. Rectangles can still be placed in more than 2-Dimensional space. 
For Example Rectangles that are placed in 3-Dimensional space will extend in the first two dimensions and be completely flat in the third.

    # Rectangles in 3-dimensional space
    index_2 = make_sps_index("example_index_prefix_2", 3, num_orthotope_dimensions=2)

Next you need to fill the index with points. 
At the moment a description for each point can be stored but not retrieved again (hence you can also just leave the description empty).
The number of coordinates you provide for each point must match the dimensionality of the index you picked in step 1.

    index.add_point((5, 7, 1, 20), "description of point for index")

If your index contains orthotope dimensions addind one datapoint requires its bottom-left-front-... point and its top-right-back-... point.
Again the dimensionality of the given points must match the picked number of dimensions.
For dimensions that are not part of the orthotope the coordinates of both given points must be equal.

For example in index_2 the third dimension is not part of the rectangle and hence equal (=1) for both points.

    index_2.add_point((5, 7, 1), (10, 8, 1), "description of point for index_2")

Once all points have been added you call `generate` to create the sparse prefix sum matrix. 
This may take a long time to compute.

    dataset_id = index.generate()

Generate returns the id of the generated dataset. 
You will need this id for querying the index later.

If you want to store multiple datsets in the same index you can provide the generate function with the index of the first point that shall be part of the dataset and the index of one past the last point that shall be part of the dataset.
Note that if your index is orthotope a single datapoint may take up multiple indices.
You should hence use `len(index)` for determining the start and end indices of your datasets.

    index = make_sps_index()

    # create dataset a
    start_a = len(index)
    index.add_point((5, 7))
    index.add_point((20, 1))
    end_a = len(index)
    dataset_a_id = index.generate(start_a, end_a)

    # create dataset b
    start_b = len(index)
    index.add_point((1, 1))
    index.add_point((9, 3))
    end_b = len(index)
    dataset_b_id = index.generate(start_b, end_b)


#### Querying indices

For querying indices you provide the bottom-left-front-.... position and the top-right-back-... position of the region you want tou count the number of points in. You also need the id of your dataset you get when generating.

    a = index.count(dataset_a_id, (0, 0), (10, 7))

The bottom-left-front-.... position is inclusive and the top-right-back-... position is exclusive.

#### Documentation

[Technical Documentation for the python module can be found here.](https://github.com/MarkusRainerSchmidt/libSps/docs/Python.html "Python Documentation")

### C++

This section contains the Documentation for using libSps from C++17.
The C++ library is header-only, so you merely need to include the appropriate header files to get it working.

#### Creating indices

For creating an instance of the Index class you need to provide the `type_defs` template parameter.
As this parameter you should provide and Instance of the `TypeDefs` class.
You can either instantiate this class by itself or use the more convenient `InMemTypeDef`, `CachedTypeDef` or `DiskTypeDef`, which are defined in sps/default.h.

Here is an example:

    #include "sps/default.h"
    #include "sps/index.h"

    typename sps::Index<InMemTypeDef<2, true, true>> xIndex( );

Functionally, this does the same as the python `make_sps_index` factory function. 

Then cou can call `addPoint` and `generate` to set up your datasets as needed.

#### Querying indices

Querying indices works the same as for the python Module.

#### Documentation

[Technical Documentation for the c++ library can be found here.](https://github.com/MarkusRainerSchmidt/libSps/docs/Cpp.html "C++ Documentation") 

### Dependent Dimensions

libSps can make that distribution of overlays in the 2nd dimension dependent on the first dimension.

**What are overlays?:**
In sparse prefix sum datastructures overlays are used to reduce storage requirements. 
To achieve this, they span a grid across the entire dataset.
For overlays to work best, points should be evenly distributed amongst the overlays.

**What does making one dimension dependent on another do?:**
Usually the overlay grid is created, so that overlays are larger in regions where points are sparse and smaller in dense regions.
This is done independently for each dimension, i.e. by sorting all points by their position in this dimension and then creating a grid division after a fixed number of points.
If dimension 2 is dependent on dimension 1, these divisions for dimension 2 are not created globally, but locally for each slice of dimension 1.
This breaks the grid in this dimension, i.e. divisions will be placed at different heights for each slice. 


**When is this usefull?:**
This is usefull if your dataset looks evenly distributed, when looking at every dimension individually, but is strongly clustered when considering the first two dimensions together.
In this case the usual approach of stretching the grid will fail to separate points evenly among all overlays and therefore produce an extremely bloated datastructure.

One example for such data is Hi-C data, where most points (~90%) lie close to the 45&deg; diagonal.
Looking seperatately at dimesion 1 or 2, all points look like they are evenly distribute, while actually they are not.

For more details you should read Schmidt et al. [2]. 

### Intervals / Rectangles / Cubes / Orthotopes

libSps can be used to store and count n-dimensional orthotopes.
1, 2, 3, ... dimensional orthotopes are intervals, rectangles, cubes, ... .
libSps is able to place an n-orthotope in a m-dimensional space (m>=n), where the orthotope is simply a flat plane in the last m-n dimensions.

`count` will only return orthotopes that are fully contained in the given region.

## Limitations

- libSps does not work with negative coordinates.
- dimension 2 is the only dimension that can be made dependent on another dimension and it can only be made dependent on dimension 1.

## Thanks

libSps uses several other projects. These are:
- [Pybind11](https://github.com/pybind/pybind11 "GitHub"): Seamless operability between C++11 and Python
- [STXXL](https://github.com/stxxl/stxxl "GitHub"): Standard Template Library for Extra Large Data Sets

## Citing libSps

For citing libSps, please use:
@todo

## References

[1] Shekelyan, M., Dignös, A. & Gamper, J. Sparse prefix sums: Constant-time range sum queries over sparse multidimensional data cubes. Information Systems 82, 136–147 (2019).

[2] Schmidt et al. @todo
