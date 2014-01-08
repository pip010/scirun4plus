IF(BUILD_BIOMESH3D)

  # Make sure the compiler can find include files from our Hellolibrary.
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${VISPACK_INCLUDE_DIRECTORIES} ${TEEM_INCLUDE_DIRS})

  ADD_EXECUTABLE(morphsmooth morphsmooth.cxx volutils.cxx morphology.cxx)
  TARGET_LINK_LIBRARIES(morphsmooth
    ${TEEM_LIBRARIES}
    matrix
    vol
    image
    jpeg
    tiff
    fitsio
    util
   ${M_LIBRARY})

ELSE(BUILD_BIOMESH3D)

  #SET(CMAKE_CXX_FLAGS -O3)
  #SET(CMAKE_CXX_FLAGS -g)

  # Make sure the compiler can find include files from our Hellolibrary.
  INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES} ${TEEM_INCLUDE_DIRS})

  # Make sure the linker can find the Hello library once it is built.
  LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib ${TEEM_LIBRARY_DIRS})

  ADD_EXECUTABLE(morphsmooth morphsmooth.cxx volutils.cxx morphology.cxx)
  TARGET_LINK_LIBRARIES(morphsmooth ${TEEM_LIBRARIES} matrix vol image jpeg tiff fitsio util lapack2 lapack1 blas F77 m)

ENDIF(BUILD_BIOMESH3D)
