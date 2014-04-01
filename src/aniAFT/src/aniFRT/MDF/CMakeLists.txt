# Ani3D/src/aniAFT/src/aniFRT/MDF


file (GLOB C_SRCS *.c)
file (GLOB CPP_SRCS *.cpp)
file (GLOB HDRS *.h)

foreach(src_file ${C_SRCS})
	set_source_files_properties(${src_file} PROPERTIES LANGUAGE C)
endforeach(src_file ${C_SRCS})

include_directories(${aniAFT3D_INCLUDES})

add_library(aniFRT3D_MDF ${C_SRCS} ${CPP_SRCS} ${HDRS})
target_link_libraries(aniFRT3D_MDF aniAFT3D ${MATH_LIB})
