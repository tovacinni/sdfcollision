cmake_minimum_required(VERSION 2.6)
project(sdfcontact LANGUAGES CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#######################################
# Variable Setup
#######################################

# libigl
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
find_package(LIBIGL REQUIRED)

#######################################
# Core Library
#######################################

include_directories("include/")
if(USE_SOLUTION)
  file(GLOB SRCFILES solution/*.cpp)
else()
  file(GLOB SRCFILES src/*.cpp)
endif()

add_library(
    core 
    ${SRCFILES}
)

target_include_directories(
    core PUBLIC
)

target_link_libraries(
    core PUBLIC 
    igl::core 
    igl::opengl 
    igl::opengl_glfw
)

#######################################
# Executable
#######################################

add_executable(
    sdfcontact
    "main.cpp"
)

target_link_libraries(
    sdfcontact 
    core
    igl::core 
    igl::opengl 
    igl::opengl_glfw 
)

