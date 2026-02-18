// File       : Common.h
// Created    : Fri Mar 14 2025 12:23:29 (+0100)
// Author     : Fabian Wermelinger
// Description: Common solver utilities
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef COMMON_H
#define COMMON_H

#include <ostream>
#include <streambuf>

namespace linearSolver
{

// Null ostream:
// https://stackoverflow.com/a/760353
template <class cT, class traits = std::char_traits<cT>>
class basic_nullbuf : public std::basic_streambuf<cT, traits>
{
    typename traits::int_type overflow(typename traits::int_type c)
    {
        return traits::not_eof(c); // indicate success
    }
};

template <class cT, class traits = std::char_traits<cT>>
class basic_onullstream : public std::basic_ostream<cT, traits>
{
public:
    basic_onullstream()
        : std::basic_ios<cT, traits>(&sbuf_),
          std::basic_ostream<cT, traits>(&sbuf_)
    {
        this->init(&sbuf_);
    }

private:
    basic_nullbuf<cT, traits> sbuf_;
};

using onullstream = basic_onullstream<char>;

} /* namespace linearSolver */

#endif /* COMMON_H */
