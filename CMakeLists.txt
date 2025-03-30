project(VCF_RLZ)
cmake_minimum_required(VERSION 3.10)

# disallow in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# enable warnings (always good)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")
set(CMAKE_BUILD_TYPE Release)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -O5")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -O5")
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
endif()

# TODO: al final hardcodear todo
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

# INCLUDE LIBRARIES
add_subdirectory(stxxl)
add_subdirectory(RLZ)

# ADD FLAGS REQUIRED
# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# ADD INCLUDES PATH OF OUR LIBRARIES
# add STXXL includes path
set(VCF_RLZ_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include/)
include_directories(${STXXL_INCLUDE_DIRS} ${RLZ_INCLUDE_DIRS} ${VCF_RLZ_INCLUDE_DIRS})

link_directories(${VCF_RLZ_INCLUDE_DIRS})
add_library(vcf_rlz STATIC ${SOURCES} ${VCF_RLZ_INCLUDE_DIRS})
# target_include_directories(vcf_rlz PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(fullbuild GenerateFromVCFFiles.cpp)
add_executable(prebuild GenerateFromPreloadedFile.cpp)
add_executable(timeC "${PROJECT_SOURCE_DIR}/experiments/ConstructionTime.cpp")
add_executable(timeF "${PROJECT_SOURCE_DIR}/experiments/FindSnippetTime.cpp")

#target_include_directories(process PUBLIC ${STXXL_INCLUDE_DIRS} ${VCF_RLZ_INCLUDE_DIRS})
target_link_libraries(vcf_rlz stxxl sdsl divsufsort divsufsort64 RLZ)
target_link_libraries(fullbuild vcf_rlz sdsl RLZ stxxl sdsl divsufsort divsufsort64)
target_link_libraries(prebuild vcf_rlz sdsl RLZ sdsl divsufsort divsufsort64)
target_link_libraries(timeC vcf_rlz sdsl RLZ stxxl sdsl divsufsort divsufsort64)
target_link_libraries(timeF vcf_rlz sdsl RLZ sdsl divsufsort divsufsort64)
