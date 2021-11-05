cmake_minimum_required(VERSION 3.11)

project(Stag)

include(GenerateExportHeader)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall -w  ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -w ")

# Check C++11 or C++0x support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenCV 4.5.2 REQUIRED COMPONENTS ximgproc)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

file(GLOB HEADER_FILES "src/*.h")
file(GLOB SRC_FILE1 "src/*.c*")
#file(GLOB SRC_FILE2 "src/ED/*.c*")


include_directories(
    ${OpenCV_INCLUDE_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR}
    src/
    #src/ED/
)

add_library(
    LibStag SHARED
    ${HEADER_FILES}
    ${SRC_FILE1}
    #${SRC_FILE2}
)

set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)
generate_export_header(LibStag)

target_link_libraries(
    LibStag
    ${OpenCV_LIBS}
)

install(TARGETS LibStag DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES
 ${HEADER_FILES}
 ${PROJECT_BINARY_DIR}/libstag_export.h
 DESTINATION ${CMAKE_INSTALL_PREFIX}/include/libstag
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/app)

add_executable(
    testrun
    main.cpp
)

target_link_libraries(
    testrun
    LibStag
)
