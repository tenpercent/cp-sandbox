# Ani3D/src/aniAFT/src/aniTET

set(aniAFT3D_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "aniAFT3D include dir")

set(aniAFT3D_C_SRCS
	aft3d.c
	ugly.c
	refine.c
	qual.c
	libaft.c
	helper.c
	cdt.c
	det.c
	check.c
	delaunay.c
	aft2d.c
)

foreach(src_file ${aniAFT3D_C_SRCS})
	set_source_files_properties(${src_file} PROPERTIES LANGUAGE C)
endforeach(src_file ${aniAFT3D_C_SRCS})

file(GLOB aniAFT3D_C_HDRS *.h)

add_library(aniAFT3D ${aniAFT3D_C_SRCS} ${aniAFT3D_C_HDRS})
target_link_libraries(aniAFT3D ${MATH_LIB})
