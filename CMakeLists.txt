PROJECT(LCSB)

CMAKE_MINIMUM_REQUIRED(VERSION 3.4)
set (CMAKE_CXX_STANDARD 14)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

AUX_SOURCE_DIRECTORY(. DIR_SRCS)
ADD_EXECUTABLE(lcsb ${DIR_SRCS}
                "pri_queue.cc"
                "util.cpp"
                "register.cpp"
                "myndarray.cc"
                "hashAlg/srp.cpp"
                "hashAlg/pivots.cpp"
                "bucketAlg/lcs.cpp"
                "app/composible_index.cpp"
                "app/ground_truth.cc"
                "app/srp_lcs.cpp"
                "app/rh_lcs.cpp"
              )

TARGET_LINK_LIBRARIES( lcsb LINK_PUBLIC ${Boost_LIBRARIES} )