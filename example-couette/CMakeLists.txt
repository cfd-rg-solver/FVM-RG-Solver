cmake_minimum_required(VERSION 3.14)

set(PROJECT_MAIN "${PROJECT_NAME}_main")
project(${PROJECT_MAIN} LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable("${PROJECT_MAIN}_1" soda-solver-lib_main.cpp)
add_executable("${PROJECT_MAIN}_2" couette-solver-lib_main.cpp)
target_link_libraries("${PROJECT_MAIN}_1" FVM-RG-Solver)
target_link_libraries("${PROJECT_MAIN}_2" FVM-RG-Solver)

