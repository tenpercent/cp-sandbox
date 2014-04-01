# Ani3D/src/aniAFT/src/aniFRT/CAD

include(${cubit_util_INCLUDES}/UseUtil.cmake)
file (GLOB C_SRCS *.c)
file (GLOB CPP_SRCS *.cpp)
file (GLOB HDRS *.h)

foreach(src_file ${C_SRCS})
	set_source_files_properties(${src_file} PROPERTIES LANGUAGE C)
endforeach(src_file ${C_SRCS})

if (WIN32)
	add_definitions(-DWNT)
endif (WIN32)
if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
	SET( OCC64 "-D_OCC64")
endif(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
set_source_files_properties( ${CPP_SRCS} COMPILE_FLAGS "-DHAVE_IOSTREAM -DHAVE_IOMANIP -DHAVE_FSTREAM -DHAVE_LIMITS_H ${OCC64}" )

include_directories(${cubit_util_INCLUDES})
include_directories(${cubit_geom_INCLUDES})
include_directories(${cubit_occ_INCLUDES})
include_directories(${OCC_INCLUDES})
include_directories(${aniFRT3D_PRM_INCLUDES})
link_directories(${OCC_LIBRARY_DIR})

add_library(aniFRT3D_CAD ${C_SRCS} ${CPP_SRCS} ${HDRS})
target_link_libraries(aniFRT3D_CAD aniFRT3D_PRM cubit_OCC cubit_geom cubit_util ${MATH_LIB})
