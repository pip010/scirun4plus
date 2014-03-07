# Recurse into the "Hello" and "Demo" subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.

IF(BUILD_VISPACK_LAPACK)
  SUBDIRS(netlib)
ELSE(BUILD_VISPACK_LAPACK)
  SUBDIRS(netlib fits tiff jpeg)
ENDIF(BUILD_VISPACK_LAPACK)
