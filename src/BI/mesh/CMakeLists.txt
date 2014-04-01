# Ani3D/src/BI/mesh

FortranCInterface_HEADER(Ani3D_FC.h MACRO_NAMESPACE Ani3DInterface SYMBOLS mbaFixShape ANI_Metric_Eucl saveMani)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

file(GLOB Mesh_HDRS *.h)
file(GLOB Mesh_SRCS *.c)
foreach(src_file ${Mesh_SRCS})
	set_source_files_properties(${src_file} PROPERTIES LANGUAGE C)
endforeach(src_file ${Mesh_SRCS})

set(Misc_SRCS
	fix.c
)


add_executable(Vena vena.c ${Misc_SRCS} ${Mesh_HDRS})
target_link_libraries(Vena aniAFT3D aniFRT3D_PRM aniMBA3D)

add_executable(Mesh mesh.c ${Misc_SRCS} ${Mesh_HDRS})
target_link_libraries(Mesh aniAFT3D aniFRT3D_PRM aniMBA3D)

add_executable(cMesh cmesh.c ${Misc_SRCS} ${Mesh_HDRS})
target_link_libraries(cMesh aniAFT3D aniFRT3D_PRM aniMBA3D)

add_executable(Cyl cyl.c ${Misc_SRCS} ${Mesh_HDRS})
target_link_libraries(Cyl aniAFT3D aniFRT3D_PRM aniMBA3D)

add_executable(cCyl ccyl.c ${Misc_SRCS} ${Mesh_HDRS})
target_link_libraries(cCyl aniAFT3D aniFRT3D_PRM aniMBA3D)

add_executable(Hspace hspace.c ${Misc_SRCS} ${Mesh_HDRS})
target_link_libraries(Hspace aniAFT3D aniFRT3D_PRM aniMBA3D)

