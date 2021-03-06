cmake_minimum_required(VERSION 3.16)

project(c-coursework)
set(CMAKE_BUILD_TYPE Release)

include_directories("/usr/include/")
# -D_USEMINGW_ANSI_STDIO=1
set(DCMAKE_CC_COMPILER=gcc)

set(GCC_COVERAGE_COMPILE_FLAGS "-L=/usr/include/ -static -Wall -Wextra -pedantic -O3")
set(GCC_COVERAGE_LINK_FLAGS "-libm")

add_executable(coursework 23762.c)
target_link_libraries(coursework m)

