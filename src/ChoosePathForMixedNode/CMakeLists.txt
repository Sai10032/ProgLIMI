cmake_minimum_required(VERSION 3.6)
project(ChoosePathForMixedNode)

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES Add_ECMP_Paths.cpp Choose_Path_For_Mixed_Node.cpp Path_Selection.cpp  Read_Graph.cpp)
add_executable(ChoosePathForMixedNode ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(ChoosePathForMixedNode ${Boost_LIBRARIES})


