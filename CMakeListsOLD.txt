cmake_minimum_required(VERSION 3.16)

project(process_PC_to_vegetation_metrics VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

if (MSVC)
    add_definitions( "/D_CRT_SECURE_NO_WARNINGS /MP /openmp" )
endif(MSVC)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/gdal_include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/GEOGRAPHICLIB/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/QUAZIP/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/LAStools_include)

file(GLOB CPP_FILES src/*.cpp)
file(GLOB H_FILES src/*.h)

set(SOURCE_LIST ${CPP_FILES} ${H_FILES})

add_executable(${PROJECT_NAME} ${SOURCE_LIST})
target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX
${CMAKE_CURRENT_SOURCE_DIR}/libs/gdal.lib
${CMAKE_CURRENT_SOURCE_DIR}/libs/GEOGRAPHICLIB/GeographicLib-i.lib
${CMAKE_CURRENT_SOURCE_DIR}/libs/LASlib.lib
)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_LIST})

# Copy necessary .dll files
add_custom_command(TARGET process_PC_to_vegetation_metrics POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${PROJECT_SOURCE_DIR}/libs/dlls_to_copy"
        ${PROJECT_SOURCE_DIR}/build/release)

        add_custom_command(TARGET process_PC_to_vegetation_metrics POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${PROJECT_SOURCE_DIR}/libs/dlls_to_copy"
            ${PROJECT_SOURCE_DIR}/build/RelWithDebInfo)