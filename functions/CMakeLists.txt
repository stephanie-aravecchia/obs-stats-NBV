cmake_minimum_required(VERSION 2.4.6)
project(functions)

# Activate GSL. For now this is used only to compute spline approximation 
# and could be replaced by Eigen in the future
set(USE_GSL OFF)


if(USE_GSL)
    find_path(GSL_INCLUDE_DIR gsl/gsl_version.h
        /usr/local/include
        /usr/include
        )
    find_library(GSL_LIBRARIES NAMES gsl PATHS /usr/lib /usr/local/lib)
    find_library(GSLCBLAS_LIBRARIES NAMES gslcblas PATHS /usr/lib /usr/local/lib)
    set(GSL_LIBRARIES ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
    add_definitions(-DUSE_GSL)
endif()


add_definitions(-Wall -g -ggdb -O3)
find_package(catkin REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME} 
  #CATKIN_DEPENDS
  #DEPENDS
)
include_directories(include 
 ${catkin_INCLUDE_DIRS}
 ${GSL_INCLUDE_DIR})
#common commands for building c++ executables and libraries
add_library(${PROJECT_NAME} 
    src/Function.cpp
    )
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES} ${GSL_LIBRARIES})

add_executable(test_function test/function.cpp)
target_link_libraries(test_function ${PROJECT_NAME} ${catkin_LIBRARIES} ${GSL_LIBRARIES})
   


if(USE_GSL)
    add_executable(test_approx test/function_approx.cpp)
    target_link_libraries(test_approx ${PROJECT_NAME} ${catkin_LIBRARIES}  ${GSL_LIBRARIES})
endif()

