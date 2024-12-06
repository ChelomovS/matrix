cmake_minimum_required(VERSION 3.11)

project(gtest)

if(NOT DEFINED INCLUDE_DIR)
    message(WARNING "INCLUDE_DIR is not defined.")
endif()

find_package(GTest REQUIRED)

set(SOURCES_FOR_TESTING
    ../src/double_compare.cpp
)

set(HEADERS 
    unit/controllanle.hpp)
aux_source_directory(unit UNIT_TESTS)

add_executable(unit_tests ${SOURCES_FOR_TESTING} ${UNIT_TESTS})

target_include_directories(unit_tests PRIVATE ${INCLUDE_DIR} ${GTEST_INCLUDE_DIRS} ${HEADERS})

enable_testing()

target_compile_features(unit_tests PUBLIC cxx_std_20)
target_link_libraries(unit_tests ${GTEST_BOTH_LIBRARIES})

add_test(NAME unit_tests COMMAND unit_tests)