cmake_minimum_required(VERSION 3.5)
project(geometric_algebra_classes)

#POLICIES
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 OLD)
endif(COMMAND cmake_policy)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -U__STRICT_ANSI__ -O3 -ftemplate-depth-1200 -Wno-switch -Wno-deprecated-declarations -Wno-int-to-pointer-cast")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -U__STRICT_ANSI__ -nostdinc++ -O3 -ftemplate-depth-1200 -Wno-switch -Wno-int-to-pointer-cast")

#PROJECT INCLUDE AND LINK
include_directories(${PROJECT_SOURCE_DIR}/versor)
include_directories(${PROJECT_SOURCE_DIR}/versor/include)
include_directories(${PROJECT_SOURCE_DIR}/versor/include/vsr)
include_directories(${PROJECT_SOURCE_DIR}/versor/include/vsr/detail)
include_directories(${PROJECT_SOURCE_DIR}/versor/include/vsr/space)
include_directories(${PROJECT_SOURCE_DIR}/versor/include/vsr/form)
#include_directories(${PROJECT_SOURCE_DIR}/versor/include/vsr/draw)
include_directories(${PROJECT_SOURCE_DIR}/versor/include/vsr/util)

#file(GLOB core /versor/src/space/*.cpp)


link_directories(${PROJECT_BINARY_DIR})

#CORE TEST TARGETS
#set( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )


set(vsr_libraries vsr)

#OPTIONAL ADD GFX AS SUBDIRECTORY TO BUILD GRAPHICS DEMOS
#add_subdirectory(versor/ext/gfx)

#GET INCLUDE DIRECTORIES, LINK DIRECTORIES, AND LINKED LIBS
#get_property(GFX_LIBRARIES GLOBAL PROPERTY GFX_LIBRARIES)
#get_property(GFX_INCLUDE_DIR GLOBAL PROPERTY GFX_INCLUDE_DIR)
#get_property(GFX_LINK_DIR GLOBAL PROPERTY GFX_LINK_DIR)

#include_directories(${GFX_INCLUDE_DIR})
#link_directories(${GFX_LINK_DIR})

#file(GLOB graphicsGL versor/src/draw/GL/*.cpp)
#file(GLOB graphicsGLES versor/src/draw/GLES/*.cpp)

#set(vsr_libraries  ${vsr_libraries} vsrGfx)
#set(libraries ${vsr_libraries} ${GFX_LIBRARIES} )

#DISABLE GFX
set(vsr_libraries  ${vsr_libraries})
set(libraries ${vsr_libraries})



add_executable(geometric_algebra_classes main.cpp data.cpp data.h EQ.cpp EQ.h Afun.cpp Afun.h numAlgo.cpp numAlgo.h)
target_link_libraries(geometric_algebra_classes ${libraries})


