project(Jacobi_MGPCG)
set(PROJECT_NAME "Jacobi_MGPCG")

cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)#Default is 14
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")

file(GLOB_RECURSE PROJECT_SRC "src/*.cpp")
include_directories(include)

#include Eigen
set(EIGEN3_INCLUDE_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})
message("${EIGEN3_INCLUDE_DIR}")

#using Apples' built in blas lapack in Accelerate framework
set(BLAS_LAPACK_LIBS "-framework Accelerate")

#link Suite-sparse
set(SUITESPARSE_INCLUDES "/usr/local/opt/suite-sparse/include")
include_directories(${SUITESPARSE_INCLUDES})
set(SUITESPARSE_LIBRARY "/usr/local/opt/suite-sparse/lib")
#set(SUITESPARSE_LIB_IND umfpack amd cholmod suitesparseconfig csparse 
	#cxsparse ldl)
list(APPEND SHARED_LIBS ${SUITESPARSE_LIBRARY})
#If no variable named <list> exists in the current scope its value is
#treated as empty and the elements are appended to that empty list. Also
#only apply to current scope.

add_executable(main ${PROJECT_SRC})
target_link_libraries(main ${SHARED_LIBS} ${BLAS_LAPACK_LIBS})



