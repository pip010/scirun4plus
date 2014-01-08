#Then for each subdirectory listed in the SUBDIRS command,
# CMakeLists.txt files are created. In the ./Hello directory, the
# following CMakeLists.txt file is created:

# Create a library called "Hello" which includes the source file "hello.cxx".
# Any number of sources could be listed here.

SET(vol_SRCS
      vol.cxx
      volume.cxx
      volumefile.cxx
      volindexlist.cxx
      march.cxx
   )

IF(BUILD_BIOMESH3D)
  IF(WIN32)
    # Heavily templated source code does not play well with
    # exporting DLL symbols.
    ADD_LIBRARY(vol
                STATIC
                ${vol_SRCS}
               )
  ELSE(WIN32)
    ADD_LIBRARY(vol ${vol_SRCS})
  ENDIF(WIN32)
ELSE(BUILD_BIOMESH3D)

  ADD_LIBRARY(vol ${vol_SRCS})

ENDIF(BUILD_BIOMESH3D)

TARGET_LINK_LIBRARIES(vol util image)

INCLUDE_DIRECTORIES(${VISPACK_PACKAGES_DIR}/tiff
                    ${VISPACK_PACKAGES_DIR}/fits
                    ${VISPACK_PACKAGES_DIR}/jpeg)
