# Ani3D/src/Tutorials/PackageINB

add_executable       (aniINBsimple main_simple.f)
target_link_libraries(aniINBsimple aniINB3D aniILU3D)

add_executable       (aniINBbratu  main_bratu.f)
target_link_libraries(aniINBbratu  aniINB3D aniILU3D)

if(WIN32)
   set_target_properties(aniINBsimple PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(aniINBbratu  PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
endif(WIN32)
