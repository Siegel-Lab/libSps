# libSps - Manual

## Installation

There are two ways to install libSps: Via bioconda and via GitHub/CMake.
The bioconda installation is easier, but restricts you to using the default compilation parameters (see below).

### PIP

The easiest way to install libSps is via pip

    pip install libsps

### Default Compilation Parameters

By default libSps is compiled to handle the following types of data:

| Parameter | Default |
|-|-|
| Dimensions | 2 and 3 dimensional |
| Intervals / Rectangles / Cubes | Intervals |
| Storage type | RAM |

Note that these parameters can be changed using the GitHub/CMake installation method.

### GitHub/CMake

To compile libSps from the githug source use the following commands:

    # clone repository
    git clone https://github.com/Siegel-Lab/libSps.git
    cd libSps

    # create and activate conda environment
    conda env create -f conda_conf/dev_env_linux.yml
    conda activate sm_dev

Then you can choose to compile using pip or CMake.

Pip:

    pip install -e .

CMake:

    mkdir build
    cd build
    cmake ..
    make

libSps compiles under Windows, Linux, and MacOS.
Various parameters of libSps are set during compile time.
We chose to use compile time parameters as the underlying memory layout is affected and knowing the layout ahead of time improves runtime performance.
If you plan using libSps from C++ you can set all these options via template parameters (see "Creating indices" below).
However, for python you have to define the parameters at compile time.

To set parameters in cmake use the following syntax while initializing the build directory:

    cmake .. -D<parameter>=<value>

To set parameters for pip export the parameters as environment variables before building:

    export SPS_<parameter>=<value>

You can configure several parameters individually for several indices (the indices are distinguished by appending _A, _B, _C, ... to the parameter name).:

- set `DIMENSIONS_A=X` to create an index with X dimensions.
- set `ORTHOTOPE_A=X` to create indices the first X dimensions store hyperrectangles data instead of points. I.e. with DIMENSIONS_A=3 and ORTHOTOPE_A=2 rectangles are stored in 3D space.
- set `STORAGE_A=X` to pick the storage system for the index, where 0=Ram, 1=Disk, and 2=Cached. Ram indices are stored in main memory only, disk indices load their data from the disk into ram on startup, and cached indices are stored on disk but keep a cache in main memory.

If multiple indices are created, they are made available via the `make_sps_index()` factory function.


## Usage

The libSps workflow is seperated in two phases: Creating an index and querying the index. 

### Python

This section contains the documentation for using libSps from Python3.
libSps should be imported as follows:

    from sps import make_sps_index

#### Creating indices

First you need to pick the correct implementation of the sparse prefix sum index for your use-case.
For this, it is easiest to call the `make_sps_index` factory function with the appropriate parameters.

For example creating an Index for 4-dimensional points, that 
- can deal with a large amount of points, and
- saves all data in files using a cache.

use the following parameters:

    index = make_sps_index("example_index_prefix", 4, storage_type="Cached")

If instead of points, you want to query intervals, rectangles, cubes or orthotopes of higher dimensionality, adjust the num_orthotope_dimensions parameter.
Note that e.g. rectangles can still be placed in more than 2-dimensional space. 
For example rectangles that are placed in 3-dimensional space will extend in the first two dimensions and be completely flat in the third.

    # Rectangles in 3-dimensional space
    index_2 = make_sps_index("example_index_prefix_2", 3, num_orthotope_dimensions=2)

Next you need to fill the index with points. 
At the moment a description for each point can be stored but not retrieved again (hence you can also just leave the description empty).
The number of coordinates you provide for each point must match the dimensionality of the index you picked in step 1.

    index.add_point((5, 7, 1, 20)) # point at position 5/7/1/20

If your index contains orthotope dimensions (see Manual->Usage->Intervals) adding one datapoint requires its bottom-left-front-... point and its top-right-back-... point.
Again the dimensionality of the given points must match the picked number of dimensions.
For dimensions that are not orthotope the coordinates of both given points must be equal.

For example in index_2 the third dimension is not part of the rectangle and hence equal (in this case: 1) for both points.

    index_2.add_point((5, 7, 1), (10, 8, 1)) # rectangle from 5/7 to 10/8 at z=1

Once all points have been added you call `generate` to create the sparse prefix sum matrix. 
This may take a long time to compute.

    dataset_id = index.generate()

Generate returns the id of the generated dataset. 
You will need this id for querying the index later.

If you want to store multiple datsets in the same index you can call generate multiple times times. Each time it will use the points that have been added since the last call to generate and return a new ID that identifies the dataset.

    index = make_sps_index()

    # create dataset a
    index.add_point((5, 7), 1) # point at position 5/7 with value 1
    index.add_point((20, 1), 5) # point at position 20/1 with value 5
    dataset_a_id = index.generate()

    # create dataset b
    index.add_point((1, 1), 2) # point at position 1/1 with value 2
    index.add_point((9, 3), 7) # point at position 9/3 with value 7
    dataset_b_id = index.generate()


#### Querying indices

For querying indices you provide the bottom-left-front-.... position and the top-right-back-... position of the region you want to count the number of points in. You also need the id of your dataset you got while calling generate.

    a = index.count(dataset_a_id, (0, 0), (10, 7))

The bottom-left-front-.... position is inclusive and the top-right-back-... position is exclusive.

#### Documentation

Technical Documentation for the python module can be found [here](https://libsps.readthedocs.io/en/stable/Python.html "Python Documentation").

### C++

This section contains the Documentation for using libSps from C++17.
The C++ library is header-only, so you merely need to include the appropriate header files to get it working.

@todo make a small example project that includes libSps

@todo give example code for the vec generator and sorter

#### Creating indices

For creating an instance of the Index class you need to provide the `type_defs` template parameter.
As this parameter you should provide an instance of the `TypeDefs` class.
You can either instantiate this class by itself or use the more convenient `InMemTypeDef`, `CachedTypeDef` or `DiskTypeDef`, which are defined in sps/default.h.

Here is an example:

    #include "sps/default.h"
    #include "sps/index.h"

    typename sps::Index<InMemTypeDef<2, 0>> xIndex( );

Functionally, this does the same as the python `make_sps_index` factory function. 

Then cou can call `addPoint` and `generate` to set up your datasets as needed.

#### Querying indices

Querying indices works the same as for the python Module.

#### Documentation

Technical Documentation for the c++ library can be found [here](https://libsps.readthedocs.io/en/stable/Cpp.html "C++ Documentation").

### Intervals

libSps can be used to store and count n-dimensional orthotopes.
1, 2, 3, ... dimensional orthotopes are intervals, rectangles, cubes, ... .
libSps is able to place an n-orthotope in m-dimensional space (m>=n), where the orthotope is simply a flat plane in the last m-n dimensions.

`count` will only return orthotopes that are fully contained in the given region.

In order count orthotopes, libSps will store the prefix sums of the orthotopes corners seperately.
I.e. (assuming we store rectangles) for each position there will be a prefix sum for the bottom-left, bottom-right, top-left, and top-right corners of the stored rectangles.
While querying in a given rectangle, we can then use one of these 4 values in each corner of the queried rectangle to only count the stored rectangles that are fully within the queried rectangle.
This does work if and only if there is no stored rectangle that fully encloses the queried rectangle.
For being able to work with arbitrary data libSps hence adds one "helper"-dimension for each orthotope dimension, where orthrotopes are placed according to their size in the respective orthotope dimension.
While querying we can then use this helper-dimension to exclude rectangles that are big enough to enclose the queried area.
It is always fine to exclude all rectangles that are big enough to enclose the queried area from the get-go, since these rectangles are too large to be fully contained by the queried area and therefore never need to be counted.

For more details you should read Schmidt et al. [2].

### Format specification

libSps can store indices on disk permanentely.
To achieve this, several files are created.
All these files are essentially vectors of varying length, where the individual elements are objects with sizes depending on the template parameters of the used type_defs (c++) or parameters of the make_sps_index function (Python). 
These object sizes are not checked at runtime, so it is up to you, to ensure that the files created by an index with parameter set X,
are always opened with this parameter set X. 
Merely the storage type (i.e. Cached or Disk) can be changed for each open operation.

In the following, we will list the objects stored in each index. 
The size of the vectors is not stored and needs to be inferred from the file size and object size. in the following td::... referst to one of these template parameters.


| file | object | desc |
|------|--------|------|
| .points | **point** | The individual points that were added. Contains the points of all datasets. The order of points may not reflect the order they have been added in. |
| .prefix_sums | **prefix_sum_array** | The prefix sum for one position in space. |
| .coords | **sparse_coordinate** | The translation from real to sparse coordinates. Each object is of type td::coordinate_t. |
| .overlays | **overlay** | The overlay grid. |
| .datasets | **dataset** | The individual datasets. |


A **point** is composed of:
- its position: a td::D sized array of values of type td::coordinate_t.
- a pointer to it's descritpion of type size_t.
- if td::IS_ORTHOTOPE, each add_point call adds all corners of the described orthotope. \
  In this case each of the added points also stores a value of type uint8_t that indicates which corner it refers to. \
  Corners are counted binary-like, i.e. in inverse order and from start to end position. \
  E.g. for 3 dimensions (left-right, down-up, front-back): left-down-front (000), left-down-back (001), left-up-front (010), left-up-back (011), right-down-front (100), ... \
  See the Intervals section for more details.


A **prefix_sum_array** is composed of:
- 2^td::ORTHOTOPE_DIMS many values of type td::val_t. \
  If there are orthotope dimensions, the prefix sums for all corners of the orthotopes are stored seperately. \
  See the Intervals section for more details.


An **overlay** is composed of:
- td::D * (td::D - 1) many **sparse_coordinate_vector_interval** s. \
  These store the translation from real to sparse space for the overlay cells. \
  Each block of td::D - 1 intervals belongs to one dimension *n* of the overlay and stores the translation to sparse space for all other *n-1* dimensions.
- td::D many **sparse_coordinate_vector_interval**. \
  These store the translation from real to sparse space, individually for each dimension, for prefix sums of points within this overlay.
- td::D many **prefix_sum_vector_interval** s, where N = td::D - 1. \
  These store the overlay cells' prefix sums for each dimension *n* of the overlay.
- one **prefix_sum_vector_interval**, where N = td::D. \
  Stores the prefix sums for the points that belong to this overlay.
- one **point_vector_interval**. \
  Stores the points that belong to this overlay.


A **sparse_coordinate_vector_interval** is composed of:
- index of type td::coordinate_t.
- start of type td::coordinate_t.
- end of type td::coordinate_t.

Here start and end give the coordinates of the interval in real space that is mapped to sparse space.
Index gives the start index of the interval in .coords.
The size of the interval is end-start.
To map a coordinate X from real to sparse space one has to lookup the value in .coords at position index+X-start. 

A **prefix_sum_vector_interval**, with a given N, is composed of:
- N many axis sizes of type td::coordinate_t.
- index of type td::coordinate_t.

Axis sizes give the size of each dimension of an N-dimensional grid.
the index i of a point p is computed as follows: 

    set i = 0.
    for n from 0 to N:
        set i = i * axis_size[n] + p[n].

Index gives the start index of the interval in .prefix_sums.


A **point_vector_interval** is composed of:
- a start index of type td::coordinate_t. \
  Points to the first point in .points of this interval.
- a end index of type td::coordinate_t. \
  Points to one past the last point in .points of this interval.

A **dataset** is composed of:
- td::D many **sparse_coordinate_vector_interval**. \
  These store the translation from real to sparse space, individually for each dimension, for the overlays of this dataset.
- one **sparse_coordinate_vector_interval_array**. \
  Stores the sparse coordinates of dimension 2 in case this dimension is dependent on dimension 1. \
  Otherwise this object may be uninitialized.
- one **overlay_vector_interval**, where N = td::D.

A **sparse_coordinate_vector_interval_array** is composed of:
- index of type td::coordinate_t.
- start of type td::coordinate_t.
- end of type td::coordinate_t.
- amount of type td::coordinate_t.

Works the same as sparse_coordinate_vector_interval, expect for storing amount many coordinate translations.
All stored coordinate translation must have the same start and end indices.
Further, they need to be consecutive in .coords.

A **overlay_vector_interval**, with a given N, is composed of:
- N many axis sizes of type td::coordinate_t.
- index of type td::coordinate_t.

Axis sizes give the size of each dimension of an N-dimensional grid.
the index i of a point p is computed as follows: 

    set i = 0.
    for n from 0 to N:
        set i = i * axis_size[n] + p[n].

Index gives the start index of the interval in .overlays.

#### Default sizes

| type | default size in bytes | desc |
|-|-|-|
| td::D | n/a | Number of dimensions. Note: Making one dimension orthotope always also adds a non orthotope dimension to the index. These extra dimensions are also included in td::D |
| td::coordinate_t | 8 bytes (uint64_t) | An individual coordinate of a point. |
| td::ORTHOTOPE_DIMS | n/a | Number of orthotope dimensions. Note: Making one dimension orthotope always also adds a non orthotope dimension to the index. These extra dimensions are also included in td::D |
| td::val_t | 4 bytes (uint32_t) | A single prefix sum. |


## Limitations

- libSps does not work with negative coordinates.

## Thanks

libSps uses several other projects. These are:
- [Pybind11](https://github.com/pybind/pybind11 "GitHub"): Seamless operability between C++11 and Python
- [STXXL](https://github.com/stxxl/stxxl "GitHub"): Standard Template Library for Extra Large Data Sets

## Citing libSps

For citing libSps, please use:

Markus R Schmidt, Anna Barcons-Simon, Claudia Rabuffo, T Nicolai Siegel, Smoother: on-the-fly processing of interactome data using prefix sums, Nucleic Acids Research, 2024; gkae008, https://doi.org/10.1093/nar/gkae008
## References

[1] Shekelyan, M., Dignös, A. & Gamper, J. Sparse prefix sums: Constant-time range sum queries over sparse multidimensional data cubes. Information Systems 82, 136–147 (2019).

[2] Markus R Schmidt, Anna Barcons-Simon, Claudia Rabuffo, T Nicolai Siegel, Smoother: on-the-fly processing of interactome data using prefix sums, Nucleic Acids Research, 2024; gkae008, https://doi.org/10.1093/nar/gkae008

## Version

This Documentation was generated for libSps |libSpsVersion|.