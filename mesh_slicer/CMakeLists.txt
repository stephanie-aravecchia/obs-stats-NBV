cmake_minimum_required(VERSION 2.4.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_BUILD_TYPE DEBUG)
#set(CMAKE_BUILD_TYPE RELWITHDEBINFO)


set(PROJECT_NAME crystals)
set($ENV{CPATH} "")

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)


find_path(BOOST_INCLUDE_DIR boost/thread.hpp
    /usr/local/include
    /usr/include
)
set(BOOST_LIBS boost_system boost_thread boost_filesystem boost_date_time)

find_package(OpenCV REQUIRED)

find_package(CGAL)
# If using  CGAL
SET(USE_CGAL YES)
IF(USE_CGAL)
    find_path(CGAL_INCLUDE_DIR CGAL/Polygon_2.h
        /usr/local/include
        /usr/include
    )
    add_definitions(-DUSE_CGAL  -frounding-math)
    set(CGAL_LIBS CGAL::CGAL)
ENDIF(USE_CGAL)
# Otherwise, PCL must be linked against.

# eigen 2 or 3
find_path(EIGEN_INCLUDE_DIR Eigen/Core
	/usr/local/include/eigen3
	/usr/local/include/eigen2
	/usr/local/include/eigen
	/usr/include/eigen3
	/usr/include/eigen2
	/usr/include/eigen
	/opt/ros/diamondback/stacks/geometry/eigen/include
)

FIND_PACKAGE(Ceres REQUIRED PATHS /usr/local /usr)

add_compile_options(-std=c++14)
add_definitions(-g -ggdb -Wall -O3 -DUSE_GETOPT)


# include all libs so far
include_directories(
    include
    ${PROJECT_SOURCE_DIR}/include 
    ${EIGEN_INCLUDE_DIR} 
    ${CERES_INCLUDES}
    ${CGAL_INCLUDE_DIR}
    ${BOOST_INCLUDE_DIR}
    ${OpenCV_INCLUDE_DIR}
)

link_directories(${CGAL_LIB_DIR})

add_executable(aabb_test src/AABB_polyhedron_facet_distance_example.cpp)
target_link_libraries(aabb_test ${CGAL_LIBS})

add_executable(shortestpath_test src/shortest_paths.cpp)
target_link_libraries(shortestpath_test ${CGAL_LIBS})

add_executable(mesh_segment_intersection_test src/mesh_segment_intersection.cpp)
target_link_libraries(mesh_segment_intersection_test ${CGAL_LIBS})

add_executable(mesh_segment_slicer src/mesh_segment_slicer.cpp)
target_link_libraries(mesh_segment_slicer ${CGAL_LIBS})

add_executable(mesh_segment_vslicer src/mesh_segment_vslicer.cpp)
target_link_libraries(mesh_segment_vslicer ${CGAL_LIBS})

add_executable(off2off src/off2off.cpp)
target_link_libraries(off2off ${CGAL_LIBS})

add_executable(polyhedron_translate src/polyhedron_translate.cpp)
target_link_libraries(polyhedron_translate ${CGAL_LIBS})

add_executable(clean_polyhedron src/clean_polyhedron.cpp)
target_link_libraries(clean_polyhedron ${CGAL_LIBS})

add_executable(optimize_path src/optimize_path.cpp)
target_link_libraries(optimize_path ${CERES_LIBRARIES} gflags)

add_executable(mesh_voxelizer src/mesh_voxelizer.cpp)
target_link_libraries(mesh_voxelizer ${CGAL_LIBS} ${OpenCV_LIBRARIES})
