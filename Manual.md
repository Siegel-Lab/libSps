# libSps - Manual

## Installation

There are two ways to install libSps: Via bioconda and via GitHub/CMake.
The bioconda installation is easier, but restricts you to using the default compilation parameters (see below).

### BioConda

    conda -c bioconda install libSps

### Default Compilation Parameters

By default libSps is compiled to handle the following types of data:

| Parameter | Default |
|-|-|
| Dimensions | 2 and 3 dimensional |
| Dependent dimension | off |
| Intervals / Rectangles / Cubes | off |
| Storage type | main memory and cached |

Note that these parameters can be changed using the GitHub/CMake installation method.

### GitHub/CMake

To compile libSps from the githug source use the following commands:

    # clone repository
    git clone https://github.com/MarkusRainerSchmidt/libSps
    cd libSps

    # create and activate conda environment (alternativeley, install the packages in conda_env/libSps.yaml yourself)
    ./conda_env/create_env.sh
    conda activate libSps

    # configure build 
    mkdir build
    cd build
    cmake ../libSps/ # Here you can configure various options (see below)

    # build
    make

Various parameters of libSps are set during compile time.
We chose to use compile time parameters as the underlying memory layout is affected and as this improves runtime performance.
If you plan using libSps from C++ you can set all these options via template parameters (see "Creating indices" below).
However, for python you have to define the parameters at compile time.
With CMake, you can configure libSps in various ways:

Add `-DWITH_PYTHON=ON` to create a shared object file that can be imported as a python module. If this parameter is set:
- set `-DNUM_DIMENSIONS_A=X -DNUM_DIMENSIONS_B=Y ... -DNUM_DIMENSIONS_D=Z` to create indices with X, Y, ..., Z dimensions.
- set `-DW_DEPENDENT_DIM=ON` to create indices where dimension 1 of the overlay grid is dependent on its dimension 0.
- set `-DWO_DEPENDENT_DIM=ON` to create indices where the overlay grids' dimensions are independent of each other.
- set `-DW_CUBES=ON` to create indices where entries are not point-like but interval-like on the first 3 dimensions. I.e. cubes placed in >=3-dimensional space.
- set `-DW_RECTANGLES=ON` to create indices where entries are not point-like but interval-like on the first 2 dimensions. I.e. rectangles placed in >=2-dimensional space.
- set `-DW_INTERVALS=ON` to create indices where entries are not point-like but interval-like on the first dimension. I.e. intervals placed in n-dimensional space.
- set `-DW_POINTS=ON` to create indices where entries are fully point-like.
- set `-DDISK=ON` to create indices that load all data from a file on startup and store it back to the file on shutdown. Expect them to consume as much RAM as the filesize.
- set `-DCACHED=ON` to create indices that use a cache to load data from and store data to a file dynamically as needed during runtime. Expect this storage type to be slightly slower than the other two options. For large datasets this storage is necessary, as it allows the RAM usage to be independent of the amount of data stored.
- set `-DRAM=ON` to create indices that store all information in RAM and never interact with the filesystem.

If multiple contradicting options are turned on as e.g. `-DNUM_DIMENSIONS_A=7 -DNUM_DIMENSIONS_B=3` and `-DRAM=ON -DDISK=ON` indices with all valid combinations will be compiled and made available via the `make_sps_index()` factory function.
You may want to consider turning options off to save compiletime: E.g. `-DCACHED=OFF -DNUM_DIMENSIONS_A=0`.



## Usage

The libSps workflow is seperated in two phases: Creating an index and querying the index. 

### Python

This section contains the documentation for using libSps from Python3.

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

    index.add_point((5, 7, 1, 20), "description of point for index")

If your index contains orthotope dimensions (see Manual->Usage->Intervals) adding one datapoint requires its bottom-left-front-... point and its top-right-back-... point.
Again the dimensionality of the given points must match the picked number of dimensions.
For dimensions that are not orthotope the coordinates of both given points must be equal.

For example in index_2 the third dimension is not part of the rectangle and hence equal (in this case: 1) for both points.

    index_2.add_point((5, 7, 1), (10, 8, 1), "description of point for index_2")

Once all points have been added you call `generate` to create the sparse prefix sum matrix. 
This may take a long time to compute.

    dataset_id = index.generate()

Generate returns the id of the generated dataset. 
You will need this id for querying the index later.

If you want to store multiple datsets in the same index you can provide the generate function with the index of the first point that shall be part of the dataset and the index of one past the last point that shall be part of the dataset.
Note that if your index is orthotope a single datapoint may take up multiple indices.
You should hence use `len(index)` for determining the start and end indices of your datasets before and after adding points.

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

For querying indices you provide the bottom-left-front-.... position and the top-right-back-... position of the region you want to count the number of points in. You also need the id of your dataset you got while calling generate.

    a = index.count(dataset_a_id, (0, 0), (10, 7))

The bottom-left-front-.... position is inclusive and the top-right-back-... position is exclusive.

#### Documentation

Technical Documentation for the python module can be found [here](https://github.com/MarkusRainerSchmidt/libSps/docs/Python.html "Python Documentation").

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

    typename sps::Index<InMemTypeDef<2, true, true>> xIndex( );

Functionally, this does the same as the python `make_sps_index` factory function. 

Then cou can call `addPoint` and `generate` to set up your datasets as needed.

#### Querying indices

Querying indices works the same as for the python Module.

#### Documentation

Technical Documentation for the c++ library can be found [here](https://github.com/MarkusRainerSchmidt/libSps/docs/Cpp.html "C++ Documentation").

### Dependent Dimensions

libSps can distribute overlays in dimension 1 dependent on dimension 0.


<img src="../_static/distributing-overlays.png" />

*Making dimension 1 dependent on dimension 0 breaks the uniform grid. Instread, rows are placed differently in each column.*

**What are overlays?:**
Overlays are used to reduce storage requirements in sparse prefix sum datastructures. 
To achieve this, they span in a grid across the entire dataset.
For this to work best, points should be evenly distributed amongst the overlays.

**What does making one dimension dependent on another do?:**
Usually the overlay grid is created, so that overlays are larger in regions where points are sparse and smaller in dense regions.
This is done independently for each dimension, i.e. by sorting all points by their position in this dimension and then creating a grid division after a fixed number of points.
If dimension 1 is dependent on dimension 0, these divisions for dimension 1 are not created globally, but locally for each slice of dimension 0.
This breaks the grid in this dimension, i.e. divisions will be placed at different heights for each slice. 


**When is this usefull?:**
This is usefull if your dataset looks evenly distributed, when looking at every dimension individually, but is strongly clustered when considering the first two dimensions together.
In this case the usual approach of stretching the grid will fail to separate points evenly among all overlays and therefore produce an extremely bloated datastructure.

One example for such data is Hi-C data, where most points (>90%) lie close to the 45-degree diagonal.
Looking seperatately at dimesion 1 or 2, all points look like they are evenly distributed, while actually they are not.

For more details you should read Shekelyan et al. [1] and Schmidt et al. [2] or have a look at the [documentation of Smoother](https://github.com/MarkusRainerSchmidt/smoother/docs/index.html "Go to the Smoother GitHub page.").

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
| .desc | char | The description strings of all points. Each description is EOF terminated. Hence, points only need to store pointers to the start of their description. Strings are stored as individual characters (1 byte a piece). |
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
  Corers are counted binary-like, i.e. in inverse order and from start to end position. \
  E.g. for 3 dimensions (left-right, down-up, front-back): left-down-front (000), left-down-back (001), left-up-front (010), left-up-back (011), right-down-front (100), ... \
  See the [Intervals section](#intervals) for more details.


A **prefix_sum_array** is composed of:
- 2^td::ORTHOTOPE_DIMS many values of type td::val_t. \
  If there are orthotope dimensions, the prefix sums for all corners of the orthotopes are stored seperately. \
  See the [Intervals section](#intervals) for more details.


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
| td::coordinate_t | 4 bytes (uint32_t) | An individual coordinate of a point. |
| td::ORTHOTOPE_DIMS | n/a | Number of orthotope dimensions. Note: Making one dimension orthotope always also adds a non orthotope dimension to the index. These extra dimensions are also included in td::D |
| td::val_t | 4 bytes (uint32_t) | A single prefix sum. |


## Limitations

- libSps does not work with negative coordinates.
- dimension 1 is the only dimension that can be made dependent on another dimension and it can only be made dependent on dimension 0.

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
