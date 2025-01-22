cmake_minimum_required(VERSION 3.10)

# Project settings
project(Caustic_Design)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Qt5 components
find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)

# Include directories for Qt5
include_directories(${Qt5Widgets_INCLUDE_DIRS})
include_directories(${Qt5OpenGL_INCLUDE_DIRS})

# Find OpenGL
find_package(OpenGL REQUIRED)
if(NOT OPENGL_FOUND)
    message(FATAL_ERROR "OpenGL not found!")
endif()
include_directories(${OPENGL_INCLUDE_DIRS})
link_directories(${OPENGL_LIBRARY_DIRS})
add_definitions(${OPENGL_DEFINITIONS})

# Find LBFGS
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")
find_package(LBFGS REQUIRED)

# Find SuiteSparse
find_package(SUITESPARSE REQUIRED)

# Find TinyXML
find_package(TinyXML REQUIRED)

# Find CGAL
find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})

# Compiler flags
set(CMAKE_CXX_FLAGS "-frounding-math -pthread -O3 ${CMAKE_CXX_FLAGS}")

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Qt5Widgets_INCLUDE_DIRS}
    ${Qt5OpenGL_INCLUDE_DIRS}
    ${OPENGL_INCLUDE_DIR}
    ${LBFGS_INCLUDE_DIRS}
    ${SUITESPARSE_INCLUDE_DIRS}
    ${TinyXML_INCLUDE_DIRS}
)

# Header files
set(Caustic_Design_HDRS
    matrix/sparse_array.h
    matrix/sparse_matrix.h
    matrix/suite_sparse_qr.h
    console_color.h
    convex_polygon.h
    domain.h
    enriched_segment.h
    grid.h
    line_search.h
    pgm.h
    pixel.h
    primitives.h
    pw_line_search.h
    ramp.h
    random.h
    rt2.h
    scene.h
    timer.h
    types.h
    util.h
    interpolation.h
    optimal_transport.h
    gradientdescent.h
)

# Source files
set(Caustic_Design_SRCS
    matrix/sparse_array.cpp
    matrix/sparse_matrix.cpp
    assign.cpp
    energy.cpp
    glviewer.cpp
    init.cpp
    io.cpp
    main.cpp
    render.cpp
    sites.cpp
    window.cpp
    interpolation.cpp
    optimal_transport.cpp
    gradientdescent.cpp
)

# UI files
set(Caustic_Design_UIS
    caustic.ui
    dialog.ui
)

# MOC files
set(Caustic_Design_MOCS
    dialog.h
    glviewer.h
    window.h
)

# Wrap UI and MOC files
qt5_wrap_ui(UIS_WRAPPED ${Caustic_Design_UIS})
qt5_wrap_cpp(MOCS_WRAPPED ${Caustic_Design_MOCS})

# Add executable
add_executable(Caustic_Design
    ${Caustic_Design_SRCS}
    ${UIS_WRAPPED}
    ${MOCS_WRAPPED}
)

# Link libraries
target_link_libraries(Caustic_Design
    Qt5::Widgets
    Qt5::OpenGL
    ${OPENGL_LIBRARIES}
    ${LBFGS_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${TinyXML_LIBRARIES}
    ${CGAL_LIBRARIES}
    tbb
    gmp
    mpfr
    blas
    lapack
	CGAL::CGAL        # Link CGAL
)
