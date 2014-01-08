IF(BUILD_BIOMESH3D)

  # Make sure the compiler can find include files from our Hellolibrary.
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${VISPACK_INCLUDE_DIRECTORIES} ${TEEM_INCLUDE_DIRS})

  # Make sure the linker can find the Hello library once it is built.  
  LINK_DIRECTORIES(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${TEEM_LIBRARY_DIRS})

  ADD_EXECUTABLE(vol2nrrd
    vol2nrrd.cxx
    volutils.cxx
    AnalyticScalarFields.cxx
    ScalarField.cxx
    mtxlib.cxx
    ReconstructionField.cxx
    IOScalarField.cxx)
  TARGET_LINK_LIBRARIES(vol2nrrd
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


  ADD_EXECUTABLE(distance3 distance3.cxx volutils.cxx AnalyticScalarFields.cxx
  ScalarField.cxx mtxlib.cxx ReconstructionField.cxx IOScalarField.cxx)
  TARGET_LINK_LIBRARIES(distance3
  ${TEEM_LIBRARIES}
  matrix 
  vol 
  image 
  jpeg 
  tiff 
  fitsio 
  util 
  lapack2 
  lapack1 
  blas 
  F77 
  m)

  ADD_EXECUTABLE(distance_parallel distance_parallel.cxx volutils.cxx AnalyticScalarFields.cxx
  ScalarField.cxx mtxlib.cxx ReconstructionField.cxx IOScalarField.cxx)
  TARGET_LINK_LIBRARIES(distance_parallel
  ${TEEM_LIBRARIES}
  matrix 
  vol 
  image 
  jpeg 
  tiff 
  fitsio 
  util 
  lapack2 
  lapack1 
  blas 
  F77 
  m)


  ADD_EXECUTABLE(vol2nrrd vol2nrrd.cxx volutils.cxx AnalyticScalarFields.cxx
  ScalarField.cxx mtxlib.cxx ReconstructionField.cxx IOScalarField.cxx)
  TARGET_LINK_LIBRARIES(vol2nrrd
  ${TEEM_LIBRARIES}
  matrix 
  vol 
  image 
  jpeg 
  tiff 
  fitsio 
  util 
  lapack2 
  lapack1 
  blas 
  F77 
  m)

ENDIF(BUILD_BIOMESH3D)



