# Recurse into the "Hello" and "Demo" subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.

#this is a comment line

IF(BUILD_BIOMESH3D)
  SUBDIRS(util image mat vol)
ELSE(BUILD_BIOMESH3D)
  SUBDIRS(util image mat vol scan voxmodel diffuse param)
ENDIF(BUILD_BIOMESH3D)

#ADD_DEFINITIONS( -g )
