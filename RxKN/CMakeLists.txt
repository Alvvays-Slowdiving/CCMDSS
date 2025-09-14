cmake_minimum_required(VERSION 3.10)
project(RxKN)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(Eigen)

file(GLOB SOURCES "src/*.cpp")

add_executable(RxKN_app ${SOURCES})

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
