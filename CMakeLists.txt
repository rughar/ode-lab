cmake_minimum_required(VERSION 3.20)

project(ode_lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(libs)
add_subdirectory(examples)

enable_testing()
add_subdirectory(tests)