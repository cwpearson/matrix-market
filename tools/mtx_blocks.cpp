// Copyright (C) 2021 Carl Pearson
// This code is released under the GPLv3 license

#include "mm/mm.hpp"
#include "kd.hpp"

#include <algorithm>
#include <limits>
#include <map>

struct Block {
    int i;
    int j;
    int pop;

    Block(int _i, int _j) : i(_i), j(_j) {} // for lower_bound / upper_bound only
    Block(int _i, int _j, int _pop) : i(_i), j(_j), pop(_pop) {}

    bool contains(int _i, int _j) const {
        return _i >= i && _i < i+16 && _j >= j && _j < j+16;
    }

    static bool by_ij(const Block &a, const Block &b) {
        if (a.i < b.i) {
            return true;
        } else if (a.i > b.i) {
            return false;
        } else {
            return a.j < b.j;
        }
    }
};

#if 0
/* count nnz covered by blocks of density
*/
int nnz_in_blocks(const int nrows, const int ncols, const std::vector<Entry> &entries,
                   const KD &kd, const float density) {

    std::vector<Block> blocks;
    // find dense blocks centered on non-zeros
    // sort so that blocks are in order

    for (const Entry &e : entries) {
        int ci = e.i - 8;
        int cj = e.j - 8;
        ci = std::min(ci, nrows-16);
        ci = std::max(ci, 0);
        cj = std::min(cj, ncols-16);
        cj = std::max(cj, 0);

        // blocks is sorted in i,j order by construction,
        // so restrict the blocks that need to be searched through
        // our our block contains up to +16 in each direction,
        // so we need to check any block that could fit those points
        auto lb = std::lower_bound(blocks.begin(), blocks.end(), Block(ci-15,cj-15), Block::by_ij);
        auto ub = std::upper_bound(blocks.begin(), blocks.end(), Block(ci+31,cj+31), Block::by_ij);
        // lb = blocks.begin();
        // ub = blocks.end();
        
        // skip this block if any point is already in a block
        bool any = false;
        for (auto bi = lb; bi < ub; ++bi) {
            if (bi->contains(ci,cj) 
                || bi->contains(ci,cj+16) 
                || bi->contains(ci+16,cj) 
                || bi->contains(ci+16,cj+16) ) {
                any=true;
                break;
            }
        }
        if (any) {
            continue;
        }

        int pop = kd.range_count(ci, ci+16, cj, cj+16);

        if (float(pop) / (16*16) >= density) {
            blocks.push_back(Block(ci,cj, pop));
            continue;
        }

    }
    int total = 0;
    for (const Block &block : blocks) {
        total += block.pop;
    }
    return total;
}
#endif

#if 0
/* count non-zeros in blocks aligned to blockSize x blockSize 
*/
int nnz_aligned_blocks(const int nrows, const int ncols, const KD &kd, const int blockSize, const float density) {
    int tot = 0;
    for (int bi = 0; bi < nrows; bi += blockSize) {
        for (int bj = 0; bj < ncols; bj += blockSize) {
            int pop = kd.range_count(bi, bi+blockSize, bj, bj+blockSize);
            if ((float(pop) / (blockSize * blockSize)) >= density) {
                tot += pop;
            }
        }
    }
    return tot;
}
#endif

using Ordinal = int64_t;
using Scalar = float;
using Offset = size_t;
using reader_t = MtxReader<Ordinal, Scalar, Offset>;
using coo_t = reader_t::coo_type;
using entry_t = coo_t::entry_type;

#if 0
/* count non-zeros in blocks aligned to blockSize x blockSize 
*/
int nnz_aligned_blocks2(const coo_t &mat, const int blockSize, const float density) {
    int tot = 0;
    size_t i_lb = 0;
    for (int bi = 0; bi < mat.num_rows(); bi += blockSize) {
        // advance pointer to first entry in row bi
        for(; mat.entries[i_lb].i < bi; ++i_lb) {}

        for (int bj = 0; bj < mat.num_cols(); bj += blockSize) {

            // count population of block at bi,bj
            int pop = 0;
            for (size_t i = i_lb; mat.entries[i].i < bi + blockSize; ++i) {
                const entry_t &e = mat.entries[i];
                if (e.j >= bj && e.j < bj + blockSize) {
                    ++pop;
                }
            }

            if ((float(pop) / (blockSize * blockSize)) >= density) {
                tot += pop;
            }

        }
    }
    return tot;
}
#endif

struct Point {
    int i;
    int j;

    bool operator<(const Point &rhs) const {
        if (i < rhs.i) {
            return true;
        } else if (i > rhs.i) {
            return false;
        } else {
            return j < rhs.j;
        }
    }
};


struct Result {
    int blocks; // number of blocks
    int nnz; // number of non-zeros in blocks
    Result() : blocks(0), nnz(0) {}
};

/* count non-zeros in blocks aligned to blockSize x blockSize 
*/
std::vector<Result> nnz_aligned_blocks3(const coo_t &mat, const int blockSize, const std::vector<float> &densities) {

    // number of blocks for each provided density
    std::vector<Result> counts(densities.size());

    // count population of all blocks
    std::map<Point, int> blockPops;
    for (const entry_t &e : mat.entries) {
        Point b;
        b.i = e.i / blockSize;
        b.j = e.j / blockSize;

        // insert p,0 if it doesn't exist.
        // either way, return iterator to blockPops[p]
        blockPops.emplace(b, 0).first->second += 1;
    }


    // sum up population of dense blocks
    for (const auto &kv : blockPops) {
        for (size_t di = 0; di < densities.size(); ++di) {
            if (float(kv.second) / (blockSize * blockSize) >= densities[di]) {
                counts[di].blocks += 1;
                counts[di].nnz += kv.second;
            }
        }
    }

    return counts;
}

int main(int argc, char **argv) {

    if (argc <= 1 ) {
        std::cerr << "USAGE: " << argv[0] << " input.mtx...\n";
    }

    const std::vector<float> densities{0.1, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};


    std::cout << "file,rows,cols,nnz";
    for (float density : densities) {
        std::cout << "," << "blocks (" << density << "),nnz (" << density << ")";
    }
    std::cout << "\n";

    // read as coo data
    for (int arg = 1; arg < argc; ++arg) {

        std::cout << argv[arg] << std::flush;
        reader_t reader(argv[arg]);
        coo_t mat;
        try {
            mat = reader.read_coo();

        } catch (const std::exception &e) {
            // on error, blank, but print failure reason
            std::cout << ",,,,,,,,," << e.what() << "\n";
            continue;
        }

        std::cout << "," << mat.num_rows() << "," << mat.num_cols()
                  << "," << mat.entries.size() << std::flush;

#if 0
        // sort by i,j
        std::sort(mat.entries.begin(), mat.entries.end(), entry_t::by_ij);
        for (float density : densities) {
            std::cout << "," << nnz_aligned_blocks2(mat, 16, density) << std::flush;
            break;
        }
        std::cout << "\n";
#endif

#if 0
        // convert to KD::Point
        std::vector<KD::Point> ps;
        for (entry_t &e : mat.entries) {
            ps.push_back(KD::Point(e.i, e.j));
        }

        // build kd tree
        KD kd(ps);

        for (float density : densities) {
            std::cout << "," << nnz_aligned_blocks(mat.num_rows(), mat.num_cols(), kd, 16, density) << std::flush;
        }
        std::cout << "\n";
#endif

#if 1
        auto counts = nnz_aligned_blocks3(mat, 16, densities);
        for (const Result &count : counts) {
            std::cout << "," << count.blocks  << "," << count.nnz << std::flush;
        }
        std::cout << "\n";
#endif


    } 
}