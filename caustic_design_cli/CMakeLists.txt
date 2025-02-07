project(caustic_design)

cmake_minimum_required(VERSION 3.14)
include(FetchContent)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(USE_CIMAGE true)

if(USE_CIMAGE)
    add_definitions(-D_cimg_)
    add_definitions(-Dcimg_display=0)

    find_package(PNG)
    if(PNG_FOUND)
        add_definitions(-Dcimg_use_png)
        set(ALLLIBS ${ALLLIBS} PNG::PNG)
    endif()

    find_package(JPEG)
    if(JPEG_FOUND)
        add_definitions(-Dcimg_use_jpeg)
        set(ALLLIBS ${ALLLIBS} JPEG::JPEG)
    endif()

   if(WIN32)
      add_definitions(-Dcimg_OS=2)
   else()
      add_definitions(-Dcimg_OS=1)
   endif()
endif()

include_directories(extern/)
include_directories(src/)

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})

FetchContent_Declare(
    libLBFGS
    GIT_REPOSITORY https://github.com/chokkan/liblbfgs.git
    GIT_TAG        master  # Or use a specific commit/tag
)

FetchContent_MakeAvailable(libLBFGS)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")
#find_package(LBFGS REQUIRED)
#include_directories(${LBFGS_INCLUDE_DIRS})

## FILES ##########################################################################

set(CD_SRC_FILES
    src/normal_integration/mesh.cpp
    src/normal_integration/normal_integration.cpp
    src/optimal_transport/optimal_transport.cpp
    src/optimal_transport/gradientdescent.cpp
    src/optimal_transport/init.cpp
    src/optimal_transport/interpolation.cpp
    src/optimal_transport/assign.cpp
)

###################################################################################

add_definitions(-DDATA_DIR="${PROJECT_SOURCE_DIR}/data")

if(NOT WIN32)
   add_definitions(-DENABLE_SSE_MODE)
else()
    add_definitions(-D_USE_MATH_DEFINES)
endif()

add_library(cdesign ${CD_SRC_FILES})



#set(APPS_SRC_FILES
    
#)

#add_library(otapputils ${APPS_SRC_FILES})

add_executable(caustic_design apps/caustic_design.cpp)
target_link_libraries(caustic_design cdesign ${ALLLIBS} ${CERES_LIBRARIES} ${CGAL_LIBRARIES} lbfgs CGAL::CGAL)
