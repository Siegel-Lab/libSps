# libSps - The Sparse Prefix Sums Library

## Overview

libSps is a versatile C++ library designed for analyzing n-dimensional data in constant time (O(1)). 
Specifically, it implements (hyper-)rectangle count queries.
This library is ideal for processing [Hi-C and RADICL-seq data](https://en.wikipedia.org/wiki/Chromosome_conformation_capture "Wikipedia").
libSps is based on algorithms developed by Shekelyan et al. [1] and Schmidt et al. [2].

This library offers seamless integration with Python3 and can be used as a header-only library for C++-17 projects.

## How it Works

The core of libSps is its ability to compute and store the prefix sums of all data points. 
Consider the example of the points 1, 3, 3, and 5:

            X       XX      X
        |   |   |   |   |   |
        0   1   2   3   4   5

The prefix sums of these points would be: 0, 1, 1, 3, 3, 4.

    4 -                     ----
    3 -             --------
    2 -
    1 -     --------
    0 - ----
        |   |   |   |   |   |
        0   1   2   3   4   5

To count the number of points in any interval (a, b], we subtract the prefix sum at position a from the prefix sum at position b. 
For example: count((1, 4]) = 3 - 1 = 2.
I.e. the two points at postion 3 but no other point are within the interval (1, 4].

The beauty of this approach is that it takes constant O(1) time, no matter the number of points in the index or the size of the queried interval.

## Getting Started

The easiest way to install libSps is via pip

    pip install libsps

libSps runs under Windows, Linux, and MacOS. Once installed, you can create an index and add points to it:

    from libSps import make_sps_index

    # Create a 2D index
    index = make_sps_index()

    # Add points
    index.add_point((5, 7), 1) # x: 5, y: 7, value: 1
    index.add_point((20, 1), 1) # x: 20, y: 1, value: 1

    # Preprocess the index
    dataset_id = index.generate()

Querying the index is just as simple:

    a = index.count(dataset_id, (0, 0), (0, 0))
    # a == 0

    b = index.count(dataset_id, (3, 0), (10, 10))
    # b == 1

The count function takes constant O(1) time, no matter the size of the queried rectangle or the number of points in the index.

## Complete Documentation

For more information and in-depth instructions, check out the [manual](https://github.com/Siegel-Lab/libSps/blob/master/Manual.md "Manual") or the [technical documentation](https://libsps.readthedocs.io/ "Technical Documentation").

## Citing libSps

For citing libSps, please use:
@todo

## References

[1] Shekelyan, M., Dignös, A. & Gamper, J. Sparse prefix sums: Constant-time range sum queries over sparse multidimensional data cubes. Information Systems 82, 136–147 (2019).

[2] Schmidt et al. @todo

