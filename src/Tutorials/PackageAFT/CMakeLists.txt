# Ani3D/src/Tutorials/PackageAFT

file(GLOB C_SRCS *.c)
foreach(c_file ${C_SRCS})
	set_source_files_properties(${c_file} PROPERTIES LANGUAGE C)
endforeach(c_file ${C_SRCS})

add_executable       (FrontToTet_aft main_aft.c)
target_link_libraries(FrontToTet_aft aniAFT3D)

if (ANI3D_USE_FRT_PRM)
    add_executable       (ParamToTet_aft main_prm.c)
    target_link_libraries(ParamToTet_aft aniFRT3D_PRM aniAFT3D)
endif (ANI3D_USE_FRT_PRM)

if (ANI3D_USE_FRT_MDF)
    add_executable       (FrtModifToTet_aft main_mdf.c)
    target_link_libraries(FrtModifToTet_aft aniFRT3D_MDF aniAFT3D)
endif (ANI3D_USE_FRT_MDF)

if (ANI3D_USE_FRT_SCG)
    add_executable       (ScgToTet_aft main_scg.c)
    target_link_libraries(ScgToTet_aft aniFRT3D_SCG aniAFT3D)
endif (ANI3D_USE_FRT_SCG)

if (ANI3D_USE_CGM_OCC)
    link_directories     (${OCC_LIBRARY_DIR})
    add_executable       (CadToTet_aft main_cad.c)
    target_link_libraries(CadToTet_aft aniFRT3D_CAD aniAFT3D)
endif (ANI3D_USE_CGM_OCC)
