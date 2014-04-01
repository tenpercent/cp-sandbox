# Ani3D/src/aniFEM

set(aniFEM3D_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "aniFEM include dir")

file(GLOB aniFEM3D_FORTRAN_SRCS *.f)
file(GLOB aniFEM3D_FORTRAN_DCLS *.fd)

add_library(aniFEM3D ${aniFEM3D_FORTRAN_SRCS} ${aniFEM3D_FORTRAN_DCLS})
target_link_libraries(aniFEM3D aniMBA3D)
