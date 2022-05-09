# libSps - The Sparse Prefix Sums Library

## Introduction

libSps is a C++ Library that allows counting points in n-dimensional space in constant time O(1). 
It was written to analyze [Hi-C and RADICL-seq data](https://en.wikipedia.org/wiki/Chromosome_conformation_capture "Wikipedia") 
and is based on the algorithms in Shekelyan et al. [1] as well as Schmidt et al. [2].

It can be loaded and used as a python3 module.

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

## Getting started quickly

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

## Manual

The full manual can be found [here](https://github.com/MarkusRainerSchmidt/libSps/Manual.md "Full Manual").

## Citing libSps

For citing libSps, please use:
@todo

## References

[1] Shekelyan, M., Dignös, A. & Gamper, J. Sparse prefix sums: Constant-time range sum queries over sparse multidimensional data cubes. Information Systems 82, 136–147 (2019).

[2] Schmidt et al. @todo
