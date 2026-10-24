cmake_minimum_required(VERSION 3.16)
project(process_PC_to_vegetation_metrics VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT WIN32)
    message(FATAL_ERROR "The bundled GDAL and LASlib libraries require Windows.")
endif()
find_package(OpenMP REQUIRED COMPONENTS CXX)

add_library(vegetation_core STATIC
    src/PointCloudMetrics.cpp src/PointCloudMetrics.h
    src/stat_functions.cpp src/stat_functions.h
    src/PolygonGeometry.cpp src/PolygonGeometry.h
    src/Forest.cpp src/Forest.h
    src/TextUtils.cpp src/TextUtils.h
    src/ForestManager.cpp src/ForestManager.h
    src/RasterIO.cpp src/RasterIO.h
    src/FeatureVector.cpp src/FeatureVector.h
)
target_include_directories(vegetation_core PUBLIC
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/libs"
    "${PROJECT_SOURCE_DIR}/libs/gdal_include"
    "${PROJECT_SOURCE_DIR}/libs/LAStools_include"
)
target_compile_definitions(vegetation_core PUBLIC NOMINMAX _USE_MATH_DEFINES)
if(MSVC)
    target_compile_definitions(vegetation_core PUBLIC _CRT_SECURE_NO_WARNINGS)
    target_compile_options(vegetation_core PRIVATE /MP)
endif()

target_link_libraries(vegetation_core PUBLIC
    OpenMP::OpenMP_CXX
    "${PROJECT_SOURCE_DIR}/libs/gdal_i.lib"
    "${PROJECT_SOURCE_DIR}/libs/LASlib.lib"
)

target_link_libraries(vegetation_core PRIVATE Normaliz)

add_library(app_support STATIC src/apps/AppSupport.cpp src/apps/AppSupport.h)
target_include_directories(app_support PUBLIC "${PROJECT_SOURCE_DIR}/src/apps")
target_link_libraries(app_support PUBLIC vegetation_core)

function(add_vegetation_app name workflow)
    add_executable(${name} src/apps/${name}_main.cpp src/apps/${workflow}.cpp)
    target_link_libraries(${name} PRIVATE app_support)
    set_target_properties(${name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/$<CONFIG>")
    add_dependencies(${name} vegetation_runtime)
endfunction()

# One shared copy step avoids parallel writes by the four executable targets.
add_custom_target(vegetation_runtime ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${PROJECT_SOURCE_DIR}/libs/dlls_to_copy" "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
    VERBATIM
)

add_vegetation_app(lidar_metrics LidarMetrics)
add_vegetation_app(segment_features SegmentFeatures)
add_vegetation_app(crop_laz CropLaz)
add_vegetation_app(inspect_forests InspectForests)

include(CTest)
if(BUILD_TESTING)
    foreach(app IN ITEMS lidar_metrics segment_features crop_laz inspect_forests)
        add_test(NAME ${app}_help COMMAND ${app} --help)
        set_tests_properties(${app}_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage:")
        add_test(NAME ${app}_missing_arguments COMMAND ${app})
        set_tests_properties(${app}_missing_arguments PROPERTIES WILL_FAIL TRUE)
    endforeach()
endif()
