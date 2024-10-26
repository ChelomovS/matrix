cmake_minimum_required(VERSION 3.11)

project(gtest)
cmake_minimum_required(VERSION 3.11)

if(NOT DEFINED INCLUDE_DIR)
    message(WARNING "INCLUDE_DIR is not defined.")
endif()

find_package(GTest REQUIRED)
enable_testing()

set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES_FOR_TESTING
    ../src/double_compare.cpp
)

aux_source_directory(unit UNIT_TESTS)

add_executable(unit_tests ${SOURCES_FOR_TESTING} ${UNIT_TESTS})

target_include_directories(unit_tests 
                            PRIVATE 
                            ${INCLUDE_DIR}
                            ${GTEST_INCLUDE_DIRS}
)

target_link_libraries(unit_tests gtest gtest_main)