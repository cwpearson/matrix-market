// Copyright (C)2021 Carl Pearson
// This code is released under the GPLv3 license

#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <complex>

struct Info
{
    enum class Format
    {
        unknown,
        COORDINATE,
        ARRAY
    };
    enum class Scalar
    {
        unknown,
        PATTERN,
        REAL,
        COMPLEX,
        INTEGER
    };
    enum class Symmetry
    {
        unknown,
        SYMMETRIC,
        SKEW,
        HERMITIAN,
        GENERAL
    };
    int nrows;
    int ncols;
    int nnz;
    Format format;
    Scalar scalar;
    Symmetry symmetry;

    Info() : nrows(-1), ncols(-1), nnz(-1), format(Format::unknown), scalar(Scalar::unknown), symmetry(Symmetry::unknown) {}

    bool operator==(const Info &rhs) const
    {
        return nrows == rhs.nrows && ncols == rhs.ncols && nnz == rhs.nnz && format == rhs.format && scalar == rhs.scalar && symmetry == rhs.symmetry;
    }
    bool operator!=(const Info &rhs) const
    {
        return !(*this == rhs);
    }

    operator bool() const
    {
        return *this != Info();
    }
};

template <typename Ordinal, typename Scalar, typename Offset = size_t>
class COO
{
private:
    Ordinal nrows_;
    Ordinal ncols_;

public:
    COO() : nrows_(0), ncols_(0) {}
    COO(Ordinal nrows, Ordinal ncols) : nrows_(nrows), ncols_(ncols) {}

    struct Entry
    {
        Ordinal i;
        Ordinal j;
        Scalar e;

        // for use with std::sort
        static bool by_ij(const Entry &a, const Entry &b)
        {
            if (a.i < b.i)
            {
                return true;
            }
            else if (a.i > a.i)
            {
                return false;
            }
            else
            {
                return a.j < b.j;
            }
        }

        Entry() = default;
        Entry(Ordinal _i, Ordinal _j, Scalar _e) : i(_i), j(_j), e(_e) {}
        bool operator==(const Entry &rhs) const
        {
            return i == rhs.i && j == rhs.j && e == rhs.e;
        }
        bool operator!=(const Entry &rhs) const
        {
            return !(*this == rhs);
        }
    };
    typedef Entry entry_type;

    std::vector<Entry> entries;

    Offset nnz() const { return Offset(entries.size()); }
    Ordinal num_rows() const { return nrows_; }
    Ordinal num_cols() const { return ncols_; }
};

/* convert `pattern` matrix to Scalar S*/
template <typename S>
S from_pattern() { return S(1); }

/* convert `real` matrix to Scalar S*/
template <typename S>
S from_real(double re) { return S(re); }

/* convert `integer` matrix to Scalar S*/
template <typename S>
S from_integer(int64_t i) { return S(i); }

/* convert `complex` matrix to Scalar S*/
template <typename S>
S from_complex(std::complex<double> c) { return S(std::abs(c)); }
template <>
std::complex<float> from_complex(std::complex<double> c) { return std::complex<float>(c.real(), c.imag()); }
template <>
std::complex<double> from_complex(std::complex<double> c) { return std::complex<double>(c.real(), c.imag()); }

/* define complex conjugate operation for non-complex types */
template <typename S>
S conj(S s) { return s; }
template <>
std::complex<float> conj(std::complex<float> s) { return std::conj(s); }
template <>
std::complex<double> conj(std::complex<double> s) { return std::conj(s); }

template <typename Ordinal, typename Scalar, typename Offset = size_t>
class MtxReader
{
private:
    Info info_;

    Info read_banner()
    {
        Info ret;
        std::ifstream inf(path_);
        if (!inf)
        {
            std::stringstream ss;
            ss << "couldn't open " << path_;
            throw std::runtime_error(ss.str());
        }
        ret = read_banner(inf);
        return ret;
    }

    static std::string to_lower(std::string s)
    {
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c)
                       { return std::tolower(c); });
        return s;
    }

    static Info read_banner(std::ifstream &inf)
    {
        Info ret;
        // read banner
        for (std::string line; std::getline(inf, line);)
        {

            if ('%' == line[0] && '%' == line[1])
            {

                std::string _;        // junk
                std::string format;   // (coordinate, array)
                std::string scalar;   // (pattern, real, complex, integer)
                std::string symmetry; // (general, symmetric, skew-symmetric, hermitian)

                // get matrix kind
                std::stringstream ss;
                ss << line;
                ss >> _; // skip '%%MatrixMarket
                ss >> _; // skip matrix
                ss >> format;
                ss >> scalar;
                ss >> symmetry;

                if ("coordinate" == format)
                {
                    ret.format = Info::Format::COORDINATE;
                }
                if ("pattern" == scalar)
                {
                    ret.scalar = Info::Scalar::PATTERN;
                }
                else if ("real" == scalar)
                {
                    ret.scalar = Info::Scalar::REAL;
                }
                else if ("complex" == scalar)
                {
                    ret.scalar = Info::Scalar::COMPLEX;
                }
                else if ("integer" == scalar)
                {
                    ret.scalar = Info::Scalar::INTEGER;
                }
                if ("symmetric" == symmetry)
                {
                    ret.symmetry = Info::Symmetry::SYMMETRIC;
                }
                else if ("general" == symmetry)
                {
                    ret.symmetry = Info::Symmetry::GENERAL;
                }
                else if ("hermitian" == to_lower(symmetry))
                {
                    ret.symmetry = Info::Symmetry::GENERAL;
                }
                else if ("skew-symmetric" == to_lower(symmetry))
                {
                    ret.symmetry = Info::Symmetry::SKEW;
                }
            }
            else if ('%' == line[0])
            {
                // skip comment
                continue;
            }
            else
            {
                //first line is matrix size, then done with banner
                std::stringstream ss;
                ss << line;
                ss >> ret.nrows;
                ss >> ret.ncols;
                ss >> ret.nnz;
                break;
            }
        }
        return ret;
    }

public:
    using coo_type = COO<Ordinal, Scalar, Offset>;
    using coo_entry_type = typename coo_type::Entry;

    MtxReader(const std::string &path) : path_(path)
    {
        info_ = read_banner();
    }

    operator bool() const
    {
        return bool(info_);
    }

    coo_type read_coo()
    {

        std::ifstream inf(path_);
        if (!inf)
        {
            throw std::logic_error("get_as_coo: couldn't open input file");
        }
        Info info = read_banner(inf);

        if (info.format == Info::Format::ARRAY)
        {
            throw std::logic_error("get_as_coo: array format");
        }

        coo_type coo(info.nrows, info.ncols);

        for (std::string line; std::getline(inf, line);)
        {
            if ('%' == line[0])
            {
                continue;
            }

            coo_entry_type entry;
            std::stringstream ss;
            ss << line;
            ss >> entry.i;
            ss >> entry.j;
            if (!ss)
            {
                throw std::logic_error("get_as_coo: unexpected format");
            }

            --entry.i;
            --entry.j;
            if (entry.i < 0 || entry.j < 0)
            {
                std::logic_error("row/col is too small (not 1-indexed?)");
            }

            // pattern has all non-zeros are 1
            switch (info.scalar)
            {
            case Info::Scalar::unknown:
                throw std::logic_error("221");
            case Info::Scalar::PATTERN:
                entry.e = from_pattern<Scalar>();
                break; // no more read
            case Info::Scalar::REAL:
            {
                double re;
                ss >> re;
                if (0.0 == re)
                    continue; // skip explicit 0
                entry.e = from_real<Scalar>(re);
                break;
            }
            case Info::Scalar::INTEGER:
            {
                int64_t i;
                ss >> i;
                if (0 == i)
                    continue; // skip explicit 0
                entry.e = from_integer<Scalar>(i);
                break;
            }
            case Info::Scalar::COMPLEX:
            {
                double real, imag;
                ss >> real;
                ss >> imag;
                if (real == 0 && imag == 0)
                    continue; // skip 0
                entry.e = from_complex<Scalar>(std::complex<double>(real, imag));
                break;
            }
            default:
                throw std::logic_error("get_as_coo: unsupported scalar type");
            }

            coo.entries.push_back(entry);

            // add any symmetric entries
            switch (info.symmetry)
            {
            case Info::Symmetry::unknown:
                throw std::logic_error("251");
            case Info::Symmetry::GENERAL:
                break;                      // no-op
            case Info::Symmetry::SYMMETRIC: // fall
            case Info::Symmetry::SKEW:
            {
                if (entry.i != entry.j)
                {
                    std::swap(entry.i, entry.j);
                    coo.entries.push_back(entry);
                }
                break;
            }
            case Info::Symmetry::HERMITIAN:
            {
                if (entry.i != entry.j)
                {
                    std::swap(entry.i, entry.j);
                    entry.e = conj(entry.e);
                    coo.entries.push_back(entry);
                }
            }
            default:
                throw std::logic_error("must be general, skew-symmetric or symmetric");
            }
        }

        return coo;
    }

private:
    std::string path_;
};