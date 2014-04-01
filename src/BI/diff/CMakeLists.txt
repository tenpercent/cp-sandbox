# Ani3D/src/BI/diff

set(Diff_SRCS
	main.f
	forlibfem.f
	solver.f
	mmax.fd
	iobin.c
	sens.c
)

set(DiffSimple_SRCS
	simple.f
	forlibfem.f
	solver.f
	mmax.fd
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../aniFEM)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../aniMBA)

add_executable(Diff ${Diff_SRCS})
add_executable(DiffSimple ${DiffSimple_SRCS})


target_link_libraries(Diff iluk aniFEM3D aniMBA3D aniVIEW3D aniILU3D ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(DiffSimple iluk aniFEM3D aniMBA3D aniVIEW3D aniILU3D ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(WIN32)
   set_target_properties(Diff PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(DiffSimple PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
endif(WIN32)
