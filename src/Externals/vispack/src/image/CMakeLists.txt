SET(image_SRCS
      image.cxx
      imagefile.cxx
      rgba.cxx
      imageRGBA.cxx
      im.cxx
      indexlist.cxx
      vispacktiffio.c
   )

IF(BUILD_BIOMESH3D)
  #Then for each subdirectory listed in the SUBDIRS command,
  # CMakeLists.txt files are created. In the ./Hello directory, the
  # following CMakeLists.txt file is created:

IF(WIN32)
  # Heavily templated source code does not play well with
  # exporting DLL symbols.
  ADD_LIBRARY(image
              STATIC
              ${image_SRCS}
             )
ELSE(WIN32)
  ADD_LIBRARY(image
              ${image_SRCS}
             )
ENDIF(WIN32)

ELSE(BUILD_BIOMESH3D)
  #Then for each subdirectory listed in the SUBDIRS command,
  # CMakeLists.txt files are created. In the ./Hello directory, the
  # following CMakeLists.txt file is created:

  # Create a library called "Hello" which includes the source file "hello.cxx".
  # Any number of sources could be listed here.
  ADD_LIBRARY(image ${image_SRCS})

ENDIF(BUILD_BIOMESH3D)

TARGET_LINK_LIBRARIES(image jpeg tiff fitsio util)

INCLUDE_DIRECTORIES(${VISPACK_PACKAGES_DIR}/tiff
                    ${VISPACK_PACKAGES_DIR}/fits
                    ${VISPACK_PACKAGES_DIR}/jpeg)
