# CMakeLists.txt example for STXXL

project(VCF_RLZ)
cmake_minimum_required(VERSION 2.8)

# disallow in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# enable warnings (always good)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")

# INCLUDE LIBRARIES
#   include the STXXL library in my-project
add_subdirectory(stxxl)
add_subdirectory(RLZ)

# ADD FLAGS REQUIRED
# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# ADD INCLUDES PATH OF OUR LIBRARIES
# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})
include_directories(${RLZ_INCLUDE_DIRS})

# build a program and link it with STXXL.
add_executable(project main.cpp)
target_link_libraries(project ${STXXL_LIBRARIES})