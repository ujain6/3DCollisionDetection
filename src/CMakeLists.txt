cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(collide LANGUAGES CXX)

include(CTest)
enable_testing()

# set(CMAKE_CXX_FLAGS "--std c++14 -lm")
set(CMAKE_CXX_FLAGS "--std c++11 -lm")

add_executable(
    collide collide.cpp
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(collide PUBLIC OpenMP::OpenMP_CXX)
endif()

##############
## TESTS
## Uncomment these if you want to test the ability of each component
##############
# add_executable( test_geometry tests/test_geometry.cpp)
# add_test( test_geometry test_geometry )
# add_executable( test_collision tests/test_collision.cpp)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(test_collision PUBLIC OpenMP::OpenMP_CXX)
# endif()
# add_test( test_collision test_collision )
# add_executable( test_read_csv tests/test_read_csv.cpp)
# add_test( test_read_csv test_read_csv )


include_directories(.)
