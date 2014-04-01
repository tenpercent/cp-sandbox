# Ani3D/src/Tutorials/PackageMBA

add_executable       (mbaAnalytic main_analytic.f)
target_link_libraries(mbaAnalytic aniMBA3D aniLMR3D)

add_executable       (mbaNodal    main_nodal.f)
target_link_libraries(mbaNodal    aniMBA3D aniLMR3D)

add_executable       (mbaReadAft  mainReadAftMesh.f)
target_link_libraries(mbaReadAft  aniMBA3D aniLMR3D)

add_executable       (mbaFixAft   mainFixAftMesh.f)
target_link_libraries(mbaFixAft   aniMBA3D aniLMR3D)

add_executable       (mbaRefine   mainRefineMesh.f)
target_link_libraries(mbaRefine   aniMBA3D aniLMR3D)


if(WIN32)
   set_target_properties(mbaAnalytic PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(mbaNodal    PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(mbaReadAft  PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(mbaFixAft   PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(mbaRefine   PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
endif(WIN32)
