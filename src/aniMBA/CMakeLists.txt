# Ani3D/src/aniMBA

set(aniMBA3D_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "aniMBA include dir")

FortranCInterface_HEADER(Ani3D_FC.h MACRO_NAMESPACE Ani3DInterface SYMBOLS seconds fpcheck)
include_directories(${CMAKE_CURRENT_BINARY_DIR})


file(GLOB aniMBA3D_FORTRAN_SRCS *.f)
file(GLOB aniMBA3D_FORTRAN_DCLS *.fd)
file(GLOB aniMBA3D_C_SRCS *.c)
foreach(src_file ${aniMBA3D_C_SRCS})
	set_source_files_properties(${src_file} PROPERTIES LANGUAGE C)
endforeach(src_file ${aniMBA3D_C_SRCS})
file(GLOB aniMBA3D_C_HDRS *.h)

add_library(aniMBA3D_C ${aniMBA3D_C_SRCS} ${aniMBA3D_C_HDRS})
add_library(aniMBA3D ${aniMBA3D_FORTRAN_SRCS} ${aniMBA3D_FORTRAN_DCLS})

target_link_libraries(aniMBA3D aniMBA3D_C ${LAPACK_LIBRARIES})
