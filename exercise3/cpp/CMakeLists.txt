cmake_minimum_required(VERSION 3.24)
project(Huckel)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

set(OpenBLAS_LIBRARIES /usr/lib/libopenblas.so)
include_directories(/usr/include/openblas)

# Add executable
add_executable(${PROJECT_NAME} main.cpp)
# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE ${OpenBLAS_LIBRARIES} OpenMP::OpenMP_CXX)