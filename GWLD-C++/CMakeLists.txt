# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.9)

# set the project name
project (GWLD)

#specify the C++ standard 
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

include_directories(include)
aux_source_directory(src SRCS)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

find_package(Armadillo REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(${ARMADILLO_INCLUDE_DIRS})
add_executable(GWLD ${SRCS})
target_link_libraries(GWLD ${ARMADILLO_LIBRARIES} OpenMP::OpenMP_CXX)