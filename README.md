# matrix-market

A C++11 single-header Matrix Market reader.

Build examples and tests
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

## Test Data

* `data/abb313.mtx`: rectangular, coordinate pattern general
* `data/08blocks.mtx`: square, coordinate integer general 
* `data/Trefethen_20b`: coordinate integer symmetric
* `data/mhd1280b`: coordinate complex Hermitian

## License

This code is released under the GPLv3 license (please see the `LICENSE` file).
The test matrix data is provided by the SuiteSparse Matrix Collection and is not covered by the license.