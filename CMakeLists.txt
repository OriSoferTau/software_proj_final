cmake_minimum_required(VERSION 3.21)
project( software_proj_final)

set(CMAKE_C_STANDARD 90)
add_compile_options(-ansi -Wall -Wextra -Werror -pedantic-errors -lm)
add_executable(kmeans spkmeans.c)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(kmeans ${PYTHON_LIBRARIES})
