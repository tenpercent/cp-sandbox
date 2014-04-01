# Ani3D/src/Tutorials/PackageILU

add_executable       (aniILU0 main_ilu0.f)
target_link_libraries(aniILU0 aniILU3D)

add_executable       (aniILU2 main_ilu2.f)
target_link_libraries(aniILU2 aniILU3D)

if(WIN32)
   set_target_properties(aniILU0 PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(aniILU2 PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
endif(WIN32)
