Documentation for the C++ Library
*********************************

.. list-table::
    :header-rows: 0

    * - :doc:`AbstractIndex<./AbstractIndex>`
      - Abstract Superclass for all Indices.
    * - :doc:`IntersectionType<./IntersectionType>`
      - An enum for querying the index.
    * - :doc:`Index<./CppIndex>`
      - The main sparse prefix sum index class.
    * - :doc:`TypeDefs<./TypeDefs>`
      - Definition of all compiletime parameters for Index.
    * - :doc:`InMemTypeDef<./InMemTypeDef>`
      - Type definitions for a RAM Index
    * - :doc:`CachedTypeDef<./CachedTypeDef>`
      - Type definitions for a Cached Index.
    * - :doc:`DiskTypeDef<./DiskTypeDef>`
      - Type definitions for a Disk Index.
    * - :doc:`DiskVecGenerator<./DiskVecGenerator>`
      - Generator class for a std::vector like.
    * - :doc:`RamVectorSorter<./RamVectorSorter>`
      - Sorter for a std::vector.
    * - :doc:`StdOutProgressStream<./StdOutProgressStream>`
      - Catches all print output.
    * - :doc:`SimpleVector<./SimpleVector>`
      - A naive, iteration-based index implementation for benchmarking.

.. toctree::
    :maxdepth: 2
    :hidden:

    AbstractIndex
    IntersectionType
    CppIndex
    TypeDefs
    InMemTypeDef
    CachedTypeDef
    DiskTypeDef
    DiskVecGenerator
    RamVectorSorter
    StdOutProgressStream
    SimpleVector
