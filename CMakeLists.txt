cmake_minimum_required(VERSION 3.20)
project(kmer_distribution)

set(CMAKE_CXX_STANDARD 14)
find_package(ZLIB)

add_executable(kmer_distribution main.cpp)
target_link_libraries(kmer_distribution ZLIB::ZLIB)