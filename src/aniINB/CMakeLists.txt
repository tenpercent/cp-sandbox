# Ani3D/src/aniINB

file(GLOB aniINB3D_FORTRAN_SRCS *.f)

add_library(aniINB3D ${aniINB3D_FORTRAN_SRCS})
target_link_libraries(aniINB3D ${LAPACK_LIBRARIES})
