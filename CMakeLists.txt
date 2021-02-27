cmake_minimum_required(VERSION 3.0)
project(GRMT)

SET(CMAKE_BUILD_TYPE Release)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

add_subdirectory(lib)

include_directories (
	${GRMT_SOURCE_DIR}/lib/bayesopt/include
	${GRMT_SOURCE_DIR}/lib/config
	${GRMT_SOURCE_DIR}/lib/inputparser
	${GRMT_SOURCE_DIR}/lib/matrix
	${GRMT_SOURCE_DIR}/lib/mydefine
	${GRMT_SOURCE_DIR}/lib/paramodel
	${GRMT_SOURCE_DIR}/lib/split
	${GRMT_SOURCE_DIR}/lib/threadpool
	${GRMT_SOURCE_DIR}/lib/treecaller
)

add_executable(grmt src/grmt.cpp)

target_link_libraries(grmt mydefine)

SET_TARGET_PROPERTIES(grmt PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${GRMT_BINARY_DIR}/bin)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
	set_property(TARGET grmt PROPERTY COMPILE_OPTIONS "-pthread")
	set_property(TARGET grmt PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
	target_link_libraries(grmt "${CMAKE_THREAD_LIBS_INIT}")
endif()

find_package(Boost REQUIRED)
if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
   find_library(Boost boost PATHS /opt/local/lib)
   include_directories(${Boost_LIBRARY_PATH})
endif()

INSTALL(TARGETS grmt
	RUNTIME DESTINATION bin
)
