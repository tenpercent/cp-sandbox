# Ani3D/src/aniILU

file(GLOB aniILU3D_FORTRAN_SRCS *.f)

add_library(aniILU3D ${aniILU3D_FORTRAN_SRCS})
target_link_libraries(aniILU3D ${BLAS_LIBRARIES})
