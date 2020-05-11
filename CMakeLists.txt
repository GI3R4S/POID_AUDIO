cmake_minimum_required(VERSION 3.10.0)

project(POID_AUDIO_PROJECT VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB_RECURSE SOURCES
  "${PROJECT_SOURCE_DIR}/src/*.cpp"
)

add_executable(POID_AUDIO main.cpp ${SOURCES})

target_include_directories(POID_AUDIO PRIVATE ${PROJECT_SOURCE_DIR}/include)

set_target_properties(
  POID_AUDIO
  PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
             LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
             RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_BUILD_TYPE Debug)
target_compile_definitions(POID_AUDIO PUBLIC RESOURCES_DIR="${PROJECT_SOURCE_DIR}/resources/")

# find_package(Python2 COMPONENTS Development NumPy)
# target_include_directories(POID PRIVATE ${Python2_INCLUDE_DIRS} ${Python2_NumPy_INCLUDE_DIRS})
# target_link_libraries(POID Python2::Python Python2::NumPy)

if (UNIX)
    target_link_libraries(POID_AUDIO pthread  X11)
endif (UNIX)