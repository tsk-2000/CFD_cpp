cmake_minimum_required(VERSION 3.20)
project(my_solver LANGUAGES CXX)

add_executable(solver
  src/main.cpp
  src/flux.cpp
  src/timeint.cpp
)

target_include_directories(solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)