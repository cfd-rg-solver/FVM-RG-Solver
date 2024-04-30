cmake_minimum_required(VERSION 3.14)

project(FVM-RG-Solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SRC_FOLDER "${PROJECT_SOURCE_DIR}/src")
set(HDR_FOLDER "${PROJECT_SOURCE_DIR}/hdr")


set(${PROJECT_NAME}_SRC  ${SRC_FOLDER}/abstractsolver.cpp
                         ${SRC_FOLDER}/macroparam.cpp
                         ${SRC_FOLDER}/godunovsolver.cpp
			             ${SRC_FOLDER}/coeffsolver.cpp
                         ${SRC_FOLDER}/mixture.cpp
                         ${SRC_FOLDER}/bordercondition.cpp
                         ${SRC_FOLDER}/startcondition.cpp
                         ${SRC_FOLDER}/datawriter.cpp
                         ${SRC_FOLDER}/lawchecker.cpp
                         ${SRC_FOLDER}/observer.cpp
                         ${SRC_FOLDER}/systemofequation.cpp
                         ${SRC_FOLDER}/riemannsolver.cpp
                         ${SRC_FOLDER}/energycalc.cpp
                         ${SRC_FOLDER}/numeric.cpp
)

set(${PROJECT_NAME}_HDR  ${HDR_FOLDER}/abstractsolver.h
                         ${HDR_FOLDER}/macroparam.h
                         ${HDR_FOLDER}/godunovsolver.h
                         ${HDR_FOLDER}/coeffsolver.h
                         ${HDR_FOLDER}/global.h    
                         ${HDR_FOLDER}/mixture.h
                         ${HDR_FOLDER}/bordercondition.h
                         ${HDR_FOLDER}/startcondition.h
                         ${HDR_FOLDER}/datawriter.h
                         ${HDR_FOLDER}/nn.h  
                         ${HDR_FOLDER}/lawchecker.h  
                         ${HDR_FOLDER}/observer.h
                         ${HDR_FOLDER}/systemofequation.h
                         ${HDR_FOLDER}/riemannsolver.h
                         ${HDR_FOLDER}/energycalc.h
                         ${HDR_FOLDER}/numeric.h
)

#add_subdirectory(example-couette)
add_subdirectory(example-shockwave)

find_package(OpenMP)
IF(OPENMP_FOUND)
    message("Found OpenMP! ^_^")
    # add flags for OpenMP
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ELSE()
    message("Missed OpenMP! x_x")
ENDIF()

add_library(${PROJECT_NAME} ${${PROJECT_NAME}_SRC} ${${PROJECT_NAME}_HDR})

TARGET_LINK_LIBRARIES(${PROJECT_NAME}  ${OpenMP_CXX_LIBRARIES} ${PROJECT_LINK_LIBS})

target_compile_definitions(${PROJECT_NAME} PRIVATE -D_USE_MATH_DEFINES)

set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")


target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/hdr>"
                                                  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>")

