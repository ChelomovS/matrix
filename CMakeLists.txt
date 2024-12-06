cmake_minimum_required(VERSION 3.15)
project(matrix)

find_package(GTest REQUIRED)
enable_testing()

set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/inc)

aux_source_directory(src SRC_FILES)

add_executable(matrix ${SRC_FILES})

target_include_directories(matrix PUBLIC ${INCLUDE_DIR} ${GTEST_INCLUDE_DIRS})

target_compile_options(matrix PRIVATE 
    $<$<CONFIG:DEBUG>:
        -Wall
        -Wextra
        -Wpedantic
        -g
        -O0
        -DDEBUG
        -Wmissing-declarations
        -Wcast-align
        -Wunused
        -Wshadow
        -Wsign-conversion
        -fsanitize=address
        -fsanitize=undefined>
)

add_subdirectory(tests)

target_compile_features(matrix PUBLIC cxx_std_20)
target_link_libraries(matrix ${GTEST_BOTH_LIBRARIES} -fsanitize=address)

