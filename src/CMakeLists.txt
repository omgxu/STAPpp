add_compile_options(-fexec-charset=UTF-8)

cmake_minimum_required(VERSION 4.0)
PROJECT(stap++)

set (CMAKE_CXX_STANDARD 11)

OPTION(STAP++_DEBUG "Output detailed results for debugging." OFF)
IF(STAP++_DEBUG)
   ADD_DEFINITIONS(-D_DEBUG_)
ENDIF()

INCLUDE_DIRECTORIES(h)

AUX_SOURCE_DIRECTORY(cpp SRC)
AUX_SOURCE_DIRECTORY(cpp/Elements SRC_ELEMENTS)
list(APPEND SRC ${SRC_ELEMENTS})
FILE(GLOB_RECURSE HEAD h/*.h)

source_group(SOURCE\ FILES FILES ${SRC})
source_group(HEADER\ FILES FILES ${HEAD})

ADD_EXECUTABLE(stap++ ${SRC} ${HEAD})
