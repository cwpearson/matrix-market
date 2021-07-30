#include "mm/mm.hpp"

template <typename Ordinal, typename Scalar, typename Offset = size_t>
int test_read(const std::string &path, Ordinal nrows, Ordinal ncols, Offset nnz)
{

    MtxReader<Ordinal, Scalar, Offset> reader(path);
    if (!reader)
    {
        std::cerr << "ERR: bad reader for " << path << "\n";
        return 1;
    }

    COO<Ordinal, Scalar, Offset> coo = reader.read_coo();

    if (nnz != coo.nnz())
    {
        std::cerr << "ERR: expected " << nnz << " got " << coo.nnz() << " nnz\n";
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
        if (1 != coo.entries[2].i)
        {
            std::cerr << "ERR1\n";
            return 1;
        }
        if (10 != coo.entries[2].j)
        {
            std::cerr << "ERR2\n";
            return 1;
        }
        if (Scalar(1) != coo.entries[2].e)
        {
            std::cerr << "ERR3\n";
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
    std::cerr << dataDir << "\n";

    if (test_read<int, float>(dataDir + "/abb313.mtx", 313, 176, 1557))
        return 1;
    if (test_read<int, std::complex<float>>(dataDir + "/abb313.mtx", 313, 176, 1557))
        return 1;

    return 0;
}