CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(LoopClosureDetector)


# Set the output directory for the build executables and libraries
ADD_DEFINITIONS(-DUNIX)
MESSAGE(STATUS "Compiling on Unix")

SET(CMAKE_CXX_FLAGS "-Wall -std=c++11") 
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG") 
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -DNDEBUG") 
SET(CMAKE_BUILD_TYPE Release)

# Find Eigen3
# FIND_PACKAGE(Eigen3 REQUIRED)
include_directories(include $ENV{EIGEN3_INCLUDE_DIR})

# Find Opencv
FIND_PACKAGE(OpenCV  REQUIRED)
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
#include_directories(include $ENV{OpenCV_INCLUDE_DIRS})
MESSAGE( STATUS "variable is: " ${OpenCV_INCLUDE_DIRS} )



add_executable(
  display_test
  src/utils.cpp
  src/points_utils.cpp
  src/display.cpp
  src/display_test.cpp
  )

target_link_libraries(
  display_test
  ${OpenCV_LIBS} 
)


add_executable(
  correct_alignment_test
  src/utils.cpp
  src/points_utils.cpp
  src/display.cpp
  src/p3p_solver.cpp
src/distance_map_correspondence_finder.cpp
src/distance_map.cpp
src/linear_relaxation_solver.cpp
  src/correct_alignment_test.cpp
)

target_link_libraries(
  correct_alignment_test
  ${OpenCV_LIBS} 
)

add_executable(
  closure_finder_test
  src/utils.cpp
  src/points_utils.cpp
  src/display.cpp
  src/p3p_solver.cpp
src/distance_map_correspondence_finder.cpp
src/distance_map.cpp
  src/closure_finder_test.cpp
)

target_link_libraries(
  closure_finder_test
  ${OpenCV_LIBS} 
)



