# Ani3D/src/Tutorials/PackageRCB

add_executable       (rcbBisect main.f)
target_link_libraries(rcbBisect aniRCB3D aniMBA3D aniVIEW3D)

if(WIN32)
   set_target_properties(rcbBisect PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
endif(WIN32)
