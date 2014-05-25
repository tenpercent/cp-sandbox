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

set(DiffThermal_SRCS
        thermoconductivity.f
        forlibfem.f
        solver.f
        mmax.fd
)

set(DiffThermal_Stationary_SRCS
        thermoconductivity_stationary.f
        lininterp.f
        pointintet.f
        det.f
        dirac.f
        solver.f
        mmax.fd
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../aniFEM)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../aniMBA)

add_executable(Diff ${Diff_SRCS})
add_executable(DiffSimple ${DiffSimple_SRCS})
add_executable(DiffThermal ${DiffThermal_SRCS})
# add_executable(TestDirac dirac.f test_dirac.f)
add_executable(DiffThermalStat ${DiffThermal_Stationary_SRCS})

target_link_libraries(Diff iluk aniFEM3D aniMBA3D aniVIEW3D aniILU3D ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(DiffSimple iluk aniFEM3D aniMBA3D aniVIEW3D aniILU3D ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(DiffThermal iluk aniFEM3D aniMBA3D aniVIEW3D aniILU3D ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
# target_link_libraries(TestDirac iluk aniFEM3D aniMBA3D aniVIEW3D aniILU3D ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(DiffThermalStat iluk aniFEM3D aniMBA3D aniVIEW3D aniILU3D ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

if(WIN32)
# not on win32 anyway; fix later
   set_target_properties(Diff PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(DiffSimple PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(DiffThermal PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(TestDirac PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
endif(WIN32)
