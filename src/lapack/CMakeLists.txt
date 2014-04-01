# Ani3D/src/lapack

file(GLOB lapack_FORTRAN_SRCS *.f)

add_library(lapack_ani3d ${lapack_FORTRAN_SRCS})

target_link_libraries(lapack_ani3d ${BLAS_LIBRARIES})

