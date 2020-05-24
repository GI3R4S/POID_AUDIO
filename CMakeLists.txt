cmake_minimum_required(VERSION 3.10.0)

project(POID_AUDIO_PROJECT VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

set(PYTHON_LIBRARY "/usr/lib/python2.7/config-x86_64-linux-gnu/libpython2.7.so")
set(PYTHON_INCLUDE_DIR "/usr/include/python2.7")

add_executable(POID_AUDIO main.cpp ${SOURCES})

target_include_directories(POID_AUDIO PRIVATE ${PROJECT_SOURCE_DIR}/include
                                              ${PYTHON_INCLUDE_DIR})

set_target_properties(
  POID_AUDIO
  PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
             LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
             RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_BUILD_TYPE Debug)
target_compile_definitions(
  POID_AUDIO PUBLIC RESOURCES_DIR="${PROJECT_SOURCE_DIR}/resources/")

target_link_libraries(POID_AUDIO ${PYTHON_LIBRARY} pthread X11)

