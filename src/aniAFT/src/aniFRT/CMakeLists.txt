# Ani3D/src/aniAFT/src/aniFRT

if (ANI3D_USE_FRT_MDF)
  add_subdirectory(MDF)
endif (ANI3D_USE_FRT_MDF)
if (ANI3D_USE_FRT_PRM)
  add_subdirectory(PRM)
  if (ANI3D_USE_CGM_OCC)
    add_subdirectory(CAD)
  endif (ANI3D_USE_CGM_OCC)
endif (ANI3D_USE_FRT_PRM)
if (ANI3D_USE_FRT_SCG)
    add_subdirectory(SCG)
endif (ANI3D_USE_FRT_SCG)
