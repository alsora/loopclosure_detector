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
FIND_PACKAGE(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

# Find Opencv
FIND_PACKAGE(OpenCV  REQUIRED)
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})




add_executable(
  display_test
  src/utils.cpp
  src/points_utils.cpp
  src/database.cpp
  src/display_test.cpp
  )

target_link_libraries(
  display_test
  ${OpenCV_LIBS} 
)


add_executable(
  alignment_test
  src/utils.cpp
  src/points_utils.cpp
  src/database.cpp
  src/solver.cpp
src/distance_map_correspondence_finder.cpp
src/distance_map.cpp
  src/alignment_test.cpp
)

target_link_libraries(
  alignment_test
  ${OpenCV_LIBS} 
)

add_executable(
  closure_finder
  src/utils.cpp
  src/points_utils.cpp
  src/database.cpp
  src/solver.cpp
src/distance_map_correspondence_finder.cpp
src/distance_map.cpp
  src/closure_finder.cpp
)

target_link_libraries(
  closure_finder
  ${OpenCV_LIBS} 
)



