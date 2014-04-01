# Ani3D/src/iluk

file(GLOB iluk_FORTRAN_SRCS *.f)

add_library(iluk ${iluk_FORTRAN_SRCS})
#target_link_libraries(iluk ${BLAS_LIBRARIES})
