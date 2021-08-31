// Copyright (C) 2021 Carl Pearson
// This code is released under the GPLv3 license

#include "mm/mm.hpp"

#include <algorithm>
#include <limits>
#include <random>
#include <iomanip>
#include <cmath>


using Ordinal = int64_t;
using Scalar = float;
using Offset = size_t;
using reader_t = MtxReader<Ordinal, Scalar, Offset>;
using coo_t = reader_t::coo_type;
using entry_t = coo_t::entry_type;

int main(int argc, char **argv) {

    if (argc <= 1 ) {
        std::cerr << "USAGE: " << argv[0] << " input.mtx...\n";
    }

    std::cout << "file,rows,cols,nnz,max abs,max nnz/row,avg nnz/row,diags,bandwidth,diagness,hopkins,err\n";

    // read as coo data
    for (int arg = 1; arg < argc; ++arg) {

        std::cout << argv[arg] << std::flush;
        reader_t reader(argv[arg]);
        coo_t res;
        try {
            res = reader.read_coo();

        } catch (const std::exception &e) {
            // on error, blank, but print failure reason
            std::cout << ",,,,,,,,,," << e.what() << "\n";
            continue;
        }

        std::cout.precision(10);
        std::cerr.precision(10);
        std::cout << "," << res.num_rows() << "," << res.num_cols() << "," << res.nnz() << std::flush;

        // find maximum absolute value of entries
        {
            Scalar sMax = -1;
            for (entry_t &e : res.entries) {
                sMax = std::max(sMax, std::abs(e.e));
            }
            std::cout << "," << sMax << std::flush;
        }

        // find max nnz per row
        {
            // histogram nnz per row
            std::vector<Ordinal> nnzs(res.num_rows(), 0);
            for (entry_t &e : res.entries) {
                ++nnzs[e.i];
            }

            // find max nnz per row
            {
                Ordinal nMax = 0;
                for (Ordinal u : nnzs) {
                    nMax = std::max(u, nMax);
                }
                std::cout << "," << nMax << std::flush;
            }
            // find avg nnz per row
            {
                uint64_t acc = 0;
                for (Ordinal u : nnzs) {
                    acc += u;
                }
                std::cout << "," << double(acc) / res.num_rows() << std::flush;
            }
        }

        // find avg
        {
            Scalar sMax = -1;
            for (entry_t &e : res.entries) {
                sMax = std::max(sMax, std::abs(e.e));
            }
            std::cout << "," << sMax << std::flush;
        }

        // count diagonal entries
        {
            Offset d = 0;
            for (size_t i = 0; i < res.entries.size(); ++i) {
                if (res.entries[i].i == res.entries[i].j) {
                    ++d;
                }
            }
            std::cout << "," << d << std::flush;
        }

        {
            // count bandwidth
            // smallest K such that A(i,j) = 0 for |i-j| > K
            Ordinal K = -1;
            for (size_t i = 0; i < res.entries.size(); ++i) {
                K = std::max(K, std::abs(res.entries[i].i - res.entries[i].j));
            }
            std::cout << "," << K;
        }


        // count diagonal-ness
        // entry e at i,j contributes e paired (i,j) samples.
        // THis probably can be extended to non-negative reals, but I'm not sure how
        // for the mean time, we just care about the position of the non-zero entries, not
        // their values, so
        // each entry contributes (i,j) and we just correlate all the x coordinates
        // with all the y coordinates
        {
            double xbar = 0;
            double ybar = 0;
            for (size_t i = 0; i < res.entries.size(); ++i) {
                xbar += res.entries[i].i;
                ybar += res.entries[i].j;
            }
            xbar /= res.entries.size();
            ybar /= res.entries.size();

            // pcc = A / (BC)
            double a=0, b=0, c=0;
            for (size_t i = 0; i < res.entries.size(); ++i) {
                a += (res.entries[i].i - xbar) * (res.entries[i].j - ybar);
                b += std::pow(res.entries[i].i - xbar, 2.0);
                c += std::pow(res.entries[i].j - ybar, 2.0);
            }
            b = std::sqrt(b);
            c = std::sqrt(c);
            const double r = a / (b*c);
            std::cout << "," << r << std::flush;
        }

        {
            // hopkins statistic
            const int m = 100;
            double su = 0;
            for (int mi = 0; mi < m; ++mi) {

                // random point in matrix
                int i = rand() % res.num_rows();
                int j = rand() % res.num_cols();

                // find closest non-zero
                double mind = std::numeric_limits<double>::infinity();
                for (size_t ei = 0; ei < res.entries.size(); ++ei) {
                    int xi = res.entries[ei].i;
                    int xj = res.entries[ei].j;
                    double d = std::sqrt(std::pow(i - xi, 2) + std::pow(j - xj, 2));
                    mind = std::min(mind, d);
                    
                }
                su += mind;
            }

            double sw = 0;
            const size_t min = 0;
            const size_t max = res.entries.size();
            std::default_random_engine generator;
            std::uniform_int_distribution<size_t> distribution(min,max);
            
            for (int mi = 0; mi < m; ++mi) {

                // random non-zero in matrix
                size_t ii =  distribution(generator);
                int i = res.entries[ii].i;
                int j = res.entries[ii].j;

                // find closest non-zero
                double mind = std::numeric_limits<double>::infinity();
                for (size_t ei = 0; ei < res.entries.size(); ++ei) {
                    if (ei == ii) {
                        continue; // skip self
                    }
                    int xi = res.entries[ei].i;
                    int xj = res.entries[ei].j;
                    double d = std::sqrt(std::pow(i - xi, 2.0) + std::pow(j - xj, 2.0));
                    mind = std::min(mind, d);
                }
                sw += mind;
            }

            std::cout << "," << su / (su + sw);
        }

        // no error
        std::cout << "," << std::endl;
    }
}