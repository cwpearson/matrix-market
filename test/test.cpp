// Copyright (C) 2021 Carl Pearson
// This code is released under the GPLv3 license

#include "mm/mm.hpp"

template <typename Ordinal, typename Scalar, typename Offset = size_t>
int test_read(const std::string &path, Ordinal nrows, Ordinal ncols, Offset nnz)
{

    typedef COO<Ordinal, Scalar, Offset> coo_type;
    typedef typename coo_type::entry_type entry_type;

    MtxReader<Ordinal, Scalar, Offset> reader(path);
    if (!reader)
    {
        std::cerr << "ERR: bad reader for " << path << "\n";
        return 1;
    }

    coo_type coo = reader.read_coo();

    if (nnz != coo.nnz())
    {
        std::cerr << "ERR: expected " << nnz << " got " << coo.nnz() << " nnz in " << path << "\n";
        return 1;
    }
    if (nrows != coo.num_rows())
    {
        std::cerr << "ERR: expected " << nrows << " got " << coo.num_rows() << " rows\n";
        return 1;
    }
    if (ncols != coo.num_cols())
    {
        std::cerr << "ERR: expected " << ncols << " got " << coo.num_cols() << " cols\n";
        return 1;
    }

    if (std::string::npos != path.find("abb313.mtx"))
    {
        if (entry_type(9, 0, from_integer<Scalar>(1)) != coo.entries[2]) {
            std::cerr << "ERR: unexpected entry in " << path << "\n";
            return 1;
        }
    }
    else if (std::string::npos != path.find("08blocks.mtx"))
    {
        if (entry_type(36, 1, from_integer<Scalar>(33)) != coo.entries[2]) {
            std::cerr << "ERR: unexpected entry in " << path << "\n";
            return 1;
        }
    }
    else if (std::string::npos != path.find("Trefethen_20b.mtx"))
    {
        if (entry_type(16, 0, from_integer<Scalar>(1)) != coo.entries[9]) {
            std::cerr << "ERR: unexpected entry in " << path << "\n";
            auto entry = coo.entries[9];
            std::cerr << entry.i << " " << entry.j << " " << entry.e << "\n";
            return 1;
        }
    }
    else if (std::string::npos != path.find("mhd1280b.mtx"))
    {
        if (entry_type(34, 1, from_complex<Scalar>(std::complex<double>(7.21908598e-5, -6.04225745e-19))) != coo.entries[5]) {
            return 1;
        }
    }
    else
    {
        std::cerr << "ERR: unexpected test file " << path << "\n";
        return 1; // unexpected test file
    }

    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        return 1;
    }

    std::string dataDir = argv[1];

    if (test_read<int, float>(dataDir + "/abb313.mtx", 313, 176, 1557))
        return 1;
    if (test_read<int, std::complex<float>>(dataDir + "/abb313.mtx", 313, 176, 1557))
        return 1;
        
    if (test_read<int, float>(dataDir + "/08blocks.mtx", 300, 300, 592))
        return 1;
    if (test_read<int, std::complex<float>>(dataDir + "/08blocks.mtx", 300, 300, 592))
        return 1;

    if (test_read<int, float>(dataDir + "/Trefethen_20b.mtx", 19, 19, 147))
        return 1;
    if (test_read<int, std::complex<float>>(dataDir + "/Trefethen_20b.mtx", 19,19, 147))
        return 1;

    if (test_read<int, float>(dataDir + "/mhd1280b.mtx", 1280, 1280, 12029))
        return 1;
    if (test_read<int, std::complex<float>>(dataDir + "/mhd1280b.mtx", 1280, 1280, 12029))
        return 1;

    return 0;
}