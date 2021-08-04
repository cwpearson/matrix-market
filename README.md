# matrix-market

A C++11 single-header Matrix Market reader (raw file [here](https://raw.githubusercontent.com/cwpearson/matrix-market/master/include/mm/mm.hpp)).

Build examples and tests (tools are always built)
```bash
mkdir build && cd build
cmake .. -DMM_BUILD_EXAMPLES=ON -DMM_BUILD_TESTS=ON
make && make test
```

Example:
```c++
#include "mm/mm.hpp"

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <file.mtx>\n";
        exit(EXIT_FAILURE);
    }
    const std::string path = argv[1];

    {
        std::cout << "with Ordinal=int, Scalar=float, Offset=size_t" << std::endl;
        typedef int Ordinal;
        typedef float Scalar;
        typedef MtxReader<Ordinal, Scalar> reader_t;
        typedef typename reader_t::coo_type coo_t;
        typedef typename coo_t::entry_type entry_t;

        // read matrix as coo
        reader_t reader(path);
        coo_t coo = reader.read_coo();

        // non-zeros, rows, cols
        std::cout << coo.nnz() << std::endl;                               // size_t
        std::cout << coo.num_rows() << "," << coo.num_cols() << std::endl; // int

        // first entry
        entry_t e = coo.entries[0];
        std::cout << e.i << "," << e.j << std::endl; // int, int
        std::cout << e.e << std::endl;               // float
    }

    {
        std::cout << "with Ordinal=int64_t, Scalar=std::complex<float>, Offset=int" << std::endl;
        typedef MtxReader<int64_t, std::complex<float>, int> reader_t;
        typedef typename reader_t::coo_type coo_t;
        typedef typename coo_t::entry_type entry_t;

        // read matrix as coo
        reader_t reader(path);
        coo_t coo = reader.read_coo();

        // non-zeros, rows, cols
        std::cout << coo.nnz() << std::endl;                               // int
        std::cout << coo.num_rows() << "," << coo.num_cols() << std::endl; // complex<double>

        // first entry
        entry_t e = coo.entries[0];
        std::cout << e.i << "," << e.j << std::endl; // int64_t
        std::cout << e.e << std::endl;               // complex<double>
    }
    return 0;
}
```

## Principles

You decide statically what your `Ordinal`, `Scalar`, and `Offset` types will be:
* `Ordinal`: the type of i, j
* `Scalar`: the type of the matrix entries
* `Offset`: the type of the number of non-zeros

You provide that to your reader. Optionally, use typedefs for convenience:
```c++
typedef int64_t Ordinal;
typedef std::complex<float> Scalar;
typedef int Offset;
typedef MtxReader<Ordinal, Scalar, Offset> reader_t;
typedef typename reader_t::coo_type coo_t;
typedef typename coo_t::entry_type entry_t;
reader_t reader(/* path to mtx file */);
```

The `Scalar` conversions will behave as you expect:
* Each `Scalar` will be constructed from the read value: `Scalar(x)`.
* `pattern` matrices will yield you `Scalar(1)` as all values
* If the matrix is `complex` and the `Scalar` is not, the `Scalar` will be equal to the magnitude of the complex number.
* If `Scalar` is complex and the matrix is not, the imaginary part of the `Scalar` will be 0. 

You can then read the matrix as a COO file
```c++
coo_t coo = reader.read_coo();
```

The COO is also templated over `Ordinal`, `Scalar`, and `Offset`.
So COO::nnz() is an `Offset`, COO::num_rows() and COO::num_cols() are `Ordinal`.

You can then iterate over entries:
```c++
for (entry_t e : coo.entries) {
    e.i;
    e.j;
    e.e;
    // ...
}
```
Entries are also tempated over `Ordinal` and `Scalar`.
Entry::i and Entry::j are `Ordinal`, Entry::e is `Scalar`.

The "additional" non-zeros implied by `hermitian`, `symmetric`, `skew-symmetric` are added to the resulting COO.
This means that coo.nnz() may not be the number of entries listed in the `.mtx` file.

## Tools

Some tools are provided

### tools/mtx-stats
Print some statistics about provided mtx files

### tools/mtx-blocks
Print counts of dense blocks in the provided matrix.

### tools/mtx-to-ppm
Convert an mtx file to a PPM file visuzalizing non-zero entries.

* `tools/mtx-to-ppm in.mtx out.ppm 2048`
  * create `out.ppm` with a visualization of the non-zeros in `in.mtx`. THe longest image dimension will be 2048 pixels.
* `tools/mtx-to-ppm in.mtx out.ppm 537 223`
  * create `out.ppm` with a visualization of the non-zeros in `in.mtx`. The image will be 537 by 223 pixels.

PPM files are large, and probably should be converted to another lossless format before distribution.
A comment will be included in the PMM file with some information about the source matrix and the mapping of non-zeros to pixel values.
For example:
```
# created by github.com/cwpearson/matrix-market/tools/mtx-to-ppm
# source matrix: 208918 x 208918 w/ 5540838 nnz
# each pixel approx 102 x 102 entries
# approx pixel's nnz count vs value
# each row is pixel value, then nnz count for val ... val+9
#   0:   482 470 459 448 437 427 417 407 397 388
#  10:   378 369 360 352 343 335 327 319 312 304
#  20:   297 290 283 276 269 263 257 251 245 239
#  30:   233 227 222 217 211 206 201 197 192 187
#  40:   183 179 174 170 166 162 158 154 151 147
#  50:   144 140 137 133 130 127 124 121 118 115
#  60:   113 110 107 105 102 100  97  95  93  91
#  70:    88  86  84  82  80  78  76  75  73  71
#  80:    69  68  66  65  63  61  60  59  57  56
#  90:    54  53  52  51  49  48  47  46  45  44
# 100:    43  42  41  40  39  38  37  36  35  34
# 110:    34  33  32  31  30  30  29  28  28  27
# 120:    26  26  25  24  24  23  23  22  22  21
# 130:    21  20  20  19  19  18  18  17  17  17
# 140:    16  16  15  15  15  14  14  14  13  13
# 150:    13  12  12  12  12  11  11  11  10  10
# 160:    10  10  10   9   9   9   9   8   8   8
# 170:     8   8   7   7   7   7   7   7   6   6
# 180:     6   6   6   6   6   5   5   5   5   5
# 190:     5   5   5   4   4   4   4   4   4   4
# 200:     4   4   4   4   3   3   3   3   3   3
# 210:     3   3   3   3   3   3   3   3   2   2
# 220:     2   2   2   2   2   2   2   2   2   2
# 230:     2   2   2   2   2   2   2   2   2   1
# 240:     1   1   1   1   1   1   1   1   1   1
# 250:     1   1   1   1   1   0
```

## Test Data

* `data/abb313.mtx`: rectangular, coordinate pattern general
* `data/08blocks.mtx`: square, coordinate integer general 
* `data/Trefethen_20b`: coordinate integer symmetric
* `data/mhd1280b`: coordinate complex Hermitian

## License

This code is released under the GPLv3 license (please see the `LICENSE` file).
The test matrix data is provided by the SuiteSparse Matrix Collection and is not covered by the license.