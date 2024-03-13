cmake_minimum_required(VERSION 3.14)

set(PROJECT_MAIN "${PROJECT_NAME}_main")
project(${PROJECT_MAIN} LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable("Soda-solver" soda-solver-lib_main.cpp)
add_executable("Couette-solver" couette-solver-lib_main.cpp)
target_link_libraries("Soda-solver" FVM-RG-Solver)
target_link_libraries("Couette-solver" FVM-RG-Solver)

