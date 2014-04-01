if (MSVC)
	add_definitions(-DWNT)
endif (MSVC)

# do util and geom and test subdirs
SUBDIRS( util geom )


