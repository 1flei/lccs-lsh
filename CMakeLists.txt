# SET(CMAKE_CXX_COMPILER /usr/bin/g++-8)
# SET(CMAKE_CXX_COMPILER clang++)
PROJECT(LCCS_LSH)

CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

#since the code uses structure binding which is the standard of C++17, 
#there is the minimum version requisite for the compiler
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0)
        message(FATAL_ERROR "GCC version must be at least 7.0!")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8)
        message(FATAL_ERROR "Clang version must be at least 8!")
    endif()
else()
    message(WARNING "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
endif()

set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wno-unused-variable -Wno-int-in-bool-context -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -fsanitize=undefined")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

FIND_PACKAGE(GSL) 
#GSL for lshkit
include_directories(${GSL_INCLUDE_DIR} )
include_directories("lshkit-0.2.1/include")
include_directories("FALCONN/src/include")
include_directories("FALCONN/external/simple-serializer")
include_directories("FALCONN/external/eigen")

aux_source_directory(. DIR_SRCS)
add_executable(lccs ${DIR_SRCS}
        "pri_queue.cc"
        "util.cpp"
        "register.cpp"
        "hashAlg/srp.cpp"
        "hashAlg/pivots.cpp"
        "app/composible_index.cpp"
        "app/ground_truth.cc"
        "app/mplsh.cpp"
        "app/falconn.cpp"
        "app/mp_lccs.cpp"
        "app/srs.cpp"
        "app/qalsh.cpp"
      )
      
add_custom_command(
  OUTPUT ${CMAKE_SOURCE_DIR}/lshkit-0.2.1/build/lib/liblshkit.a
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/lshkit-0.2.1/build
  COMMAND cmake ..
  COMMAND make -j
)

add_custom_target(lshkit
    DEPENDS ${CMAKE_SOURCE_DIR}/lshkit-0.2.1/build/lib/liblshkit.a
)

add_custom_command(
  OUTPUT ${CMAKE_SOURCE_DIR}/SRS/libsrs.a 
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/SRS
  COMMAND make -j
)

add_custom_target(srs
    DEPENDS ${CMAKE_SOURCE_DIR}/SRS/libsrs.a 
)


add_custom_command(
  OUTPUT ${CMAKE_SOURCE_DIR}/QALSH_Mem/libqalsh.a
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/QALSH_Mem
  COMMAND make -j
)

add_custom_target(qalsh
    DEPENDS ${CMAKE_SOURCE_DIR}/QALSH_Mem/libqalsh.a
)


target_link_libraries( lccs LINK_PUBLIC ${Boost_LIBRARIES}
  ${CMAKE_SOURCE_DIR}/lshkit-0.2.1/build/lib/liblshkit.a 
  ${CMAKE_SOURCE_DIR}/SRS/libsrs.a 
  ${CMAKE_SOURCE_DIR}/QALSH_Mem/libqalsh.a
  GSL::gsl GSL::gslcblas pthread)

add_dependencies(lccs lshkit srs qalsh)