#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>

#include "mm/mm.hpp"

// assume maxval is 255
void ppm_banner(std::ofstream &fs, int width, int height, std::vector<std::string> comments = {}) {
    fs << "P6"; // magic number
    fs << "\n";
    fs << "# created by github.com/cwpearson/matrix-market/tools/mtx-to-ppm\n";
    for (const std::string &c : comments) {
        fs << "# " << c << "\n";
    }
    fs << std::to_string(width);
    fs << " ";
    fs << std::to_string(height);
    fs << " ";
    fs << "255";
    fs << "\n";
}

// assume maxval is 255
// data should be raster rows, of RGB, one byte each channel
void ppm_data(std::ofstream &fs, const char *data, int width, int height) {
    fs.write(data, width*height*3);
}

using Ordinal = int64_t;
using Scalar = float;
using Offset = size_t;
using reader_t = MtxReader<Ordinal, Scalar, Offset>;
using coo_t = reader_t::coo_type;
using entry_t = coo_t::entry_type;

int main(int argc, char **argv) {
    if (argc > 5 || argc < 4) {
        std::cerr << "USAGE:\n";
        std::cerr << " " << argv[0] << "input.mtx output.pbm width height\n";
        std::cerr << " " << argv[0] << "input.mtx output.pbm maxdim\n";
        exit(1);
    }

    std::cerr << "open " << argv[2] << std::endl;
    std::ofstream outf(argv[2]);
    if (!outf) {
        std::cerr << "some error opening " << argv[2] << std::endl;
        exit(1);
    }

    std::cerr << "read " << argv[1] << std::endl;
    reader_t reader(argv[1]);
    coo_t coo = reader.read_coo();

    int width = -1, height = -1;
    if (5 == argc) {
        width = std::atoi(argv[3]);
        height = std::atoi(argv[4]);
    } else if (4 == argc) {
        if (coo.num_rows() > coo.num_cols()) {
            height = std::atoi(argv[3]);
            width = double(coo.num_cols()) * height / coo.num_rows() + 0.5;
        } else {
            width = std::atoi(argv[3]);
            height = double(coo.num_rows()) * width / coo.num_cols() + 0.5;
        }
    }

    if (width <= 0 && height <= 0) {
        std::cerr << "need to specify width and/or height > 0\n";
        exit(1);
    } else {
        std::cerr << "output image will be " << width << " x " << height << "\n";
    }

    // histogram all matrix entries
    std::vector<double> hist(width * height, 0);

    // map to image pixel
    for (const entry_t &e : coo.entries) {
        int px = double(e.j) / coo.num_cols() * width;
        int py = double(e.i) / coo.num_rows() * height;
        px = std::max(0, std::min(px, width));
        py = std::max(0, std::min(py, height));
        hist[py * width + px] += 1.0;
    }

    // log of histogram
    for (size_t i = 0; i < hist.size(); ++i) {
        double h = hist[i];
        if (h != 0) {
            h = std::log2(hist[i]);
        }
        hist[i] = h;
    }

    // normalize to 0-255
    double hMax = 0;
    for (size_t i = 0; i < hist.size(); ++i) {
        hMax = std::max(hMax, hist[i]);
    }
    std::cerr << "max pixel val: " << hMax << "\n";
    std::vector<unsigned char> data(hist.size() * 3 /*channels*/);
    for (size_t i = 0; i < hist.size(); ++i) {
        double h = std::max(0.0, std::min(hist[i] / hMax * 255, 255.0));
        data[i*3+0] = (255 - h) + 0.5;
        data[i*3+1] = (255 - h) + 0.5;
        data[i*3+2] = (255 - h) + 0.5;
    }


    // generate comments
    std::vector<std::string> comments;
    {
        // source matrix info
        std::stringstream ss;
        ss << "source matrix: " << coo.num_rows() << " x " << coo.num_cols() << " w/ " << coo.nnz() << " nnz";
        comments.push_back(ss.str());
        ss.str(""); // clear

        // pixel size
        ss << "each pixel approx " 
           << uint64_t(double(coo.num_rows()) / height + 0.5) << " x "
           << uint64_t(double(coo.num_cols()) / width  + 0.5) << " entries";
        comments.push_back(ss.str()); 
        ss.str(""); // clear

        ss << "approx pixel's nnz count vs value\n";
        ss << "# each row is pixel value, then nnz count for val ... val+9";

        const int fieldWidth = std::ceil(std::log10(std::pow(2, hMax)));

        // pixel values
        for (int i = 0; i <= 255; ++i) {

            if (i % 10 == 0) {
                ss << "\n# ";
                ss << std::setfill(' ') << std::setw(3) << i << ":  ";
            }
            uint64_t u;
            if (i == 255) {
                u = 0;
            } else {
                double h = (255 - i) / 255.0 * hMax;
                h = std::pow(2.0, h);
                u = h + 0.5; // round
            }
            ss << " " << std::setfill(' ') << std::setw(fieldWidth) << u;
        }
        comments.push_back(ss.str()); 
        ss.str(""); // clear
    }

    std::cerr << "write " << argv[2] << std::endl;
    ppm_banner(outf, width, height, comments);
    ppm_data(outf, (char*)data.data(), width, height);
    outf.close();

    return 0;

}