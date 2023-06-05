cmake_minimum_required(VERSION 3.10)
project(ROSIA_c)


include_directories(include /usr/include/aarch64-linux-gnu/atlas/)

#'/usr/include/x86_64-linux-gnu/atlas/'

set(CMAKE_CXX_STANDARD 14)


set(SRC_FILES
        src/binarytree.cpp
        src/geometry.cpp
        src/imp_rot3.cpp
        src/reg_rtree.cpp
        src/util_sort.cpp
        main.cpp)

add_executable(ROSIA_c ${SRC_FILES})
# Add the libatlas library to the linker flags

target_link_libraries(ROSIA_c -L/usr/lib/aarch64-linux-gnu -lblas -llapack -lf77blas -latlas)


