# Ani3D/src/aniLMR


file(GLOB aniLMR3D_FORTRAN_SRCS *.f)

include_directories(${aniMBA3D_INCLUDES})
include_directories(${aniFEM3D_INCLUDES})

add_library(aniLMR3D ${aniLMR3D_FORTRAN_SRCS})

target_link_libraries(aniLMR3D aniMBA3D)
