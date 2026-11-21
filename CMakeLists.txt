cmake_minimum_required(VERSION 3.18)
project(karstpy LANGUAGES C CXX)

find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

set(rmath_dir third_party/rmath)
set(rmath_sources
    bd0 chebyshev cospi d1mach dbinom dnorm dpois fmax2 fmin2 fsign gamma gammalims
    i1mach imax2 imin2 lgamma lgammacor mlutils pbeta pbinom pgamma pnorm qbinom qnorm rbinom rgamma
    rmultinom rpois sexp snorm stirlerr toms708
)
list(TRANSFORM rmath_sources PREPEND "${rmath_dir}/")
list(TRANSFORM rmath_sources APPEND ".c")

pybind11_add_module(_core NO_EXTRAS
    cpp/bindings.cpp cpp/rng.cpp cpp/tmb_ad.cpp cpp/models/karst_models.cpp
    ${rmath_sources}
)
target_include_directories(_core SYSTEM PRIVATE ${rmath_dir} third_party/eigen third_party/tmb)
target_compile_definitions(_core PRIVATE MATHLIB_STANDALONE EIGEN_MPL2_ONLY)
set_target_properties(_core PROPERTIES C_STANDARD 11 C_STANDARD_REQUIRED ON)
target_compile_features(_core PRIVATE cxx_std_17)
set_target_properties(_core PROPERTIES CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(_core PRIVATE -ffp-contract=off)
endif()

install(TARGETS _core LIBRARY DESTINATION karstpy)
