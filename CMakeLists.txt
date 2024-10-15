cmake_minimum_required(VERSION 2.8)

project(BA_compare)

#set(MATRIX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/matrix)

#find_package(Ceres REQUIRED)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")
set(DEFAULT_CXX_STANDARD ${CMAKE_CXX_STANDARD})
if (NOT DEFAULT_CXX_STANDARD)
  set(DEFAULT_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_STANDARD ${DEFAULT_CXX_STANDARD} CACHE STRING
  "C++ standard (minimum 14)" FORCE)
# Restrict CMAKE_CXX_STANDARD to the valid versions permitted and ensure that
# if one was forced via -D that it is in the valid set.
set(ALLOWED_CXX_STANDARDS 14 17 20)
set_property(CACHE CMAKE_CXX_STANDARD PROPERTY STRINGS ${ALLOWED_CXX_STANDARDS})
list(FIND ALLOWED_CXX_STANDARDS ${CMAKE_CXX_STANDARD} POSITION)
if (POSITION LESS 0)
  message(FATAL_ERROR "Invalid CMAKE_CXX_STANDARD: ${CMAKE_CXX_STANDARD}. "
    "Must be one of: ${ALLOWED_CXX_STANDARDS}")
endif()
# Specify the standard as a hard requirement, otherwise CMAKE_CXX_STANDARD is
# interpreted as a suggestion that can decay *back* to lower versions.
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "")
mark_as_advanced(CMAKE_CXX_STANDARD_REQUIRED)


set(ANN_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/ann_1.2_char/include)
set(CLAPACK_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/clapack_3.1/INCLUDE)
set(MRC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/mrc)
set(BLASWRAP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/clapack_3.1/BLAS/WRAP)
set(CBLAS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/cblas/INCLUDE)
set(SBA_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/sba_1.6)
set(MATRIX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/matrix)
set(LEVMAR_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/levmar_2.5)

#add_executable(${PROJECT_NAME} ba_demo.cpp baproblem.cpp parametersse3.cpp main.cpp)
#target_link_libraries(${PROJECT_NAME} ${CERES_LIBRARIES} )

#add_executable(${PROJECT_NAME} main.cpp parameter.h parameter.cpp BA.h BA.cpp noise.cpp noise.h)
#target_link_libraries(${PROJECT_NAME} ${CERES_LIBRARIES})


include_directories(${PROJECT_SOURCE_DIR}/lib ${ANN_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/src
        ${CLAPACK_INCLUDE_DIR} ${BLASWRAP_INCLUDE_DIR} ${CBLAS_INCLUDE_DIR} ${MRC_INCLUDE_DIR} ${SBA_INCLUDE_DIR}
        ${MATRIX_INCLUDE_DIR} ${OPENCV_INCLUDE_DIR} ${LEVMAR_INCLUDE_DIR})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_subdirectory(lib)
add_subdirectory(src)
