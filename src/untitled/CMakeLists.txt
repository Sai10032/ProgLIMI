cmake_minimum_required(VERSION 3.6)
project(untitled)

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES Add_ECMP_Paths.cpp  main.cpp  Path_Selection.cpp  Read_Graph.cpp)
add_executable(untitled ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(untitled ${Boost_LIBRARIES})


