#pragma once

#include <iostream>
#include <stdexcept>

#define Rcout std::cout
#define Rcerr std::cerr
#define eigen_assert(condition) do { if (!(condition)) throw std::out_of_range(#condition); } while (false)
#define TMBAD_ASSERT2(condition, message) do { if (!(condition)) throw std::runtime_error(message); } while (false)
#define TMBAD_ASSERT(condition) TMBAD_ASSERT2(condition, #condition)
#define TMBAD_ASSERT_DEFINED

#include <TMBad/TMBad.hpp>
#include <TMBad/tmbad_allow_comparison.hpp>
#include <TMBad/eigen_numtraits.hpp>
