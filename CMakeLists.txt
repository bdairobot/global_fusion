cmake_minimum_required(VERSION 2.8.3)
project(global_optimization)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")
#-DEIGEN_USE_MKL_ALL")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  nav_msgs
  sensor_msgs
)

find_package(Ceres REQUIRED)
find_package(Eigen3)

add_subdirectory(./ThirdParty/GeographicLib/)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ./ThirdParty/GeographicLib/include/
)


catkin_package()

add_executable(global_fusion
  src/globalNode.cpp
  src/globalOpt.cpp
  src/utility/utility.cpp)

target_link_libraries(global_fusion ${catkin_LIBRARIES} ${CERES_LIBRARIES} libGeographiccc) 