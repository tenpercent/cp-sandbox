# Ani3D/src/aniAFT/src/aniFRT/PRM

set(aniFRT3D_PRM_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "aniFrtPrm include dir")

file (GLOB C_SRCS *.c)
file (GLOB HDRS *.h)

foreach(src_file ${C_SRCS})
	set_source_files_properties(${src_file} PROPERTIES LANGUAGE C)
endforeach(src_file ${C_SRCS})

include_directories(${aniAFT3D_INCLUDES})

add_library(aniFRT3D_PRM ${C_SRCS} ${HDRS})
target_link_libraries(aniFRT3D_PRM aniAFT3D ${MATH_LIB})
