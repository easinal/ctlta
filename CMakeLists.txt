cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(RoutingFramework CXX)

# Directories containing .cc source files.
set(SOURCE_DIRECTORIES Launchers RawData)

# Determine the language standard.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Flags when building for the Devel configuration.
if(CMAKE_CXX_COMPILER_ID MATCHES GNU|Clang)
  set(CMAKE_CXX_FLAGS_DEVEL -O3)
endif()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Devel)
endif()

# Enable the compiler to use extended instructions in generated code.
if(CMAKE_CXX_COMPILER_ID MATCHES GNU|Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
option(DISABLE_AVX "Disable use of instructions in the AVX extended instruction set.")
if(DISABLE_AVX)
  if(CMAKE_CXX_COMPILER_ID MATCHES GNU|Clang)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-avx")
  endif()
endif()

# Include frame pointer information and avoid inlining for tagged functions in debug builds for better stack traces.
set(NOINLINE_TOGGLE "false")
if(CMAKE_CXX_COMPILER_ID MATCHES GNU|Clang)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
      set(NOINLINE_TOGGLE "true")
  endif()
endif()


# Request warnings.
if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
  set(FULL_WARNINGS "-Werror" "-Wfatal-errors" "-Wpedantic" "-pedantic-errors" "-Wall" "-Wextra" "-ftemplate-backtrace-limit=1" "-Wno-unknown-pragmas" "-Wno-sign-compare")
elseif (CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set(FULL_WARNINGS "-Werror" "-Wall" "-Wextra" "-pedantic-errors" "-ferror-limit=1" "-ftemplate-backtrace-limit=1" "-Wno-sign-compare")
endif ()

include_directories(${CMAKE_SOURCE_DIR})

find_package(Boost REQUIRED)

find_package(OpenMP)
find_package(CGAL REQUIRED)
find_library(cairo_LIBRARY cairo)
#message("OpenMP library path: ${OpenMP_CXX_LIBRARIES}")

include(FetchContentDependencies.cmake)
add_subdirectory(External)

# Clang does not support the option -frounding-math, thus remove it from the CGAL compile options.
if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set_target_properties(CGAL::CGAL PROPERTIES INTERFACE_COMPILE_OPTIONS "")
endif()

foreach(SRC_DIR ${SOURCE_DIRECTORIES})
  add_subdirectory(${SRC_DIR})
endforeach()
