cmake_minimum_required(VERSION 3.2)

project(a3-finite-elements-3d)

#Important CMAKE stuff
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# for profiling
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

#setup libigl properly
# libigl
option(LIBIGL_WITH_OPENGL               "Use OpenGL"            ON)
option(LIBIGL_WITH_OPENGL_GLFW          "Use GLFW"              ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI    "Use ImGui"             ON)

find_package(LIBIGL REQUIRED QUIET)
find_package(TBB REQUIRED)
# set(LIBIGL_DIR ${CMAKE_SOURCE_DIR}/extern/libigl)
# set(LIBIGL_INCLUDE_DIR ${LIBIGL_DIR}/include)
# list(APPEND CMAKE_MODULE_PATH "${LIBIGL_INCLUDE_DIR}/../cmake")
# add_subdirectory(${LIBIGL_DIR})
# include(oneTBB)

#include directories
include_directories(${PROJECT_SOURCE_DIR}/shared/include)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${LIBIGL_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extern/oneTBB)
# include_directories(${CMAKE_SOURCE_DIR}/extern/libigl/include)

#include source code
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB H_SHARED_FILES ${PROJECT_SOURCE_DIR}/shared/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB CPP_SHARED_FILES ${PROJECT_SOURCE_DIR}/shared/src/*.cpp)

#build executable 
add_executable(${PROJECT_NAME} main.cpp assignment_setup.h ${H_FILES} ${H_SHARED_FILES} ${CPP_FILES} ${CPP_SHARED_FILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)
target_link_libraries(${PROJECT_NAME} TBB::tbb)
# target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::glfw igl::imgui)


