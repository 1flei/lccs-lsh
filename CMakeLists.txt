PROJECT(LCSB)

CMAKE_MINIMUM_REQUIRED(VERSION 3.8)
SET(CMAKE_CXX_COMPILER /usr/bin/g++-8)
set (CMAKE_CXX_STANDARD 17)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-variable -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

FIND_PACKAGE(GSL) 
include_directories(${GSL_INCLUDE_DIR} )
include_directories("lshkit-0.2.1/include")

AUX_SOURCE_DIRECTORY(. DIR_SRCS)
ADD_EXECUTABLE(lcsb ${DIR_SRCS}
                "pri_queue.cc"
                "util.cpp"
                "register.cpp"
                "myndarray.cc"
                "hashAlg/srp.cpp"
                "hashAlg/pivots.cpp"
                "bucketAlg/lcs.cpp"
                "bucketAlg/lcs_int_reorder.cpp"
                "app/composible_index.cpp"
                "app/ground_truth.cc"
                "app/mplsh.cpp"

                "bucketAlg/lcs_test.cpp"
              )

TARGET_LINK_LIBRARIES( lcsb LINK_PUBLIC ${Boost_LIBRARIES} 
  ${CMAKE_SOURCE_DIR}/mih/build/libmihlib.a ${CMAKE_SOURCE_DIR}/lshkit-0.2.1/build/lib/liblshkit.a
  GSL::gsl GSL::gslcblas)