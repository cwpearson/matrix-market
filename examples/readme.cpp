// Copyright (C) 2021 Carl Pearson
// This code is released under the GPLv3 license

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