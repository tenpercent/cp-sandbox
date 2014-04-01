# Ani3D/src/blas

file(GLOB blas_FORTRAN_SRCS *.f)

add_library(blas_ani3d ${blas_FORTRAN_SRCS})
