cmake_minimum_required(VERSION 3.16)

project(process_PC_to_vegetation_metrics VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- GDAL & OSGeo4W Configuration ---
# Setting the path to your new OSGeo4W installation
set(OSGEO4W_ROOT "C:/Users/Maria/AppData/Local/Programs/OSGeo4W")

include_directories("${OSGEO4W_ROOT}/include")
link_directories("${OSGEO4W_ROOT}/lib")
# ------------------------------------

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

if (MSVC)
    add_definitions( "/D_CRT_SECURE_NO_WARNINGS /MP /openmp" )
endif(MSVC)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# REMOVED: include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/gdal_include) <- This was likely the old version
include_directories(BEFORE "${OSGEO4W_ROOT}/include")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/GEOGRAPHICLIB/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/QUAZIP/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/LAStools_include)

file(GLOB CPP_FILES src/*.cpp)
file(GLOB H_FILES src/*.h)

set(SOURCE_LIST ${CPP_FILES} ${H_FILES})

add_executable(${PROJECT_NAME} ${SOURCE_LIST})

target_link_libraries(${PROJECT_NAME} PUBLIC 
    MPI::MPI_CXX
    gdal_i                              # LINKING TO OSGeo4W GDAL
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/GEOGRAPHICLIB/GeographicLib-i.lib
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/LASlib.lib
)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_LIST})

# UPDATED: Copying DLLs from OSGeo4W bin so ZSTD and other dependencies are found
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        "${OSGEO4W_ROOT}/bin/gdal312.dll" # Note: Check your bin folder; name might be gdal309.dll or similar
        "${OSGEO4W_ROOT}/bin/zstd.dll"
	"${OSGEO4W_ROOT}/share/proj/proj.db" $<TARGET_FILE_DIR:${PROJECT_NAME}>
        $<TARGET_FILE_DIR:${PROJECT_NAME}>
    COMMENT "Copying runtime DLLs to output directory"
)