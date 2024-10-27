cmake_minimum_required(VERSION 3.11)
project(matrix)

set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GTest REQUIRED)
enable_testing()

set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/inc)

aux_source_directory(src SRC_FILES)

add_executable(matrix ${SRC_FILES})

target_include_directories(matrix PRIVATE ${INCLUDE_DIR} ${GTEST_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wpedantic -g -O0 -DDEBUG \
                          -Wmissing-declarations -Wcast-align      \
                          -Wunused -Wshadow -Wsign-conversion      \
                          -fsanitize=address -fsanitize=undefined")

set_target_properties(
    matrix PROPERTIES
    CXX_STANDARD          20
    CXX_STANDARD_REQUIRED ON
)

add_subdirectory(tests)

target_link_libraries(matrix ${GTEST_BOTH_LIBRARIES})