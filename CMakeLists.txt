cmake_minimum_required(VERSION 3.16)

project(c-coursework)
set(CMAKE_BUILD_TYPE Debug)

include_directories("/usr/include/")
# -D_USEMINGW_ANSI_STDIO=1
set(DCMAKE_CC_COMPILER=gcc)

set(GCC_COVERAGE_COMPILE_FLAGS "-pg -g -L=/usr/include/ -static -Wall -Wextra -pedantic -O3 -ggdb3 -Werror -Wdouble-promotion -Wformat=2 -Wconversion")
set(GCC_COVERAGE_LINK_FLAGS "-libm -pg")

add_executable(a.out 23762.c)
target_link_libraries(a.out m)

