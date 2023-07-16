cmake_minimum_required (VERSION 2.8)
project (mytest)
include_directories(include)

add_executable (mytest main.cpp voronoi.cpp)
find_package(CGAL REQUIRED COMPONENTS Qt5)

target_link_libraries(mytest PUBLIC gmp CGAL::CGAL_Qt5)

target_include_directories(mytest PUBLIC include)