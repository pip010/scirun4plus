IF(BUILD_BIOMESH3D)

  # Make sure the compiler can find include files from our Hellolibrary.
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${VISPACK_INCLUDE_DIRECTORIES} ${NRRDIO_INCLUDE})

  # Make sure the linker can find the Hello library once it is built.
  LINK_DIRECTORIES(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${NRRDIO_LIB})

  ADD_EXECUTABLE(linear_test linear_test.cxx)
  TARGET_LINK_LIBRARIES(linear_test matrix
    util
    lapack2
    lapack1 
    blas
    F77
    ${M_LIBRARY})

  ADD_EXECUTABLE(linear_exp linear_exp.cxx)
  TARGET_LINK_LIBRARIES(linear_exp
    matrix
    util
    lapack2
    lapack1 
    blas
    F77
    ${M_LIBRARY})

  ADD_EXECUTABLE(linear_run linear_run.cxx)
  TARGET_LINK_LIBRARIES(linear_run
    matrix
    util
    lapack2
    lapack1 
    blas
    F77
    ${M_LIBRARY})

  ADD_EXECUTABLE(linear_labels linear_labels.cxx)
  TARGET_LINK_LIBRARIES(linear_labels
    matrix
    util
    lapack2
    lapack1
    blas
    F77
    ${M_LIBRARY})

  ADD_EXECUTABLE(linear_labels2 linear_labels2.cxx)
  TARGET_LINK_LIBRARIES(linear_labels2
    matrix
    util
    lapack2
    lapack1
    blas
    F77
    ${M_LIBRARY})

  ADD_EXECUTABLE(linear_labels3 linear_labels3.cxx)
  TARGET_LINK_LIBRARIES(linear_labels3
    matrix
    util
    lapack2
    lapack1
    blas
    F77
    ${M_LIBRARY})

  ADD_EXECUTABLE(linear_labels4 linear_labels4.cxx)
  TARGET_LINK_LIBRARIES(linear_labels4
    matrix
    util
    lapack2
    lapack1
    blas
    F77
    ${M_LIBRARY})

  ADD_EXECUTABLE(linear_labels_restrict linear_labels_restrict.cxx)
  TARGET_LINK_LIBRARIES(linear_labels_restrict
    matrix
    util
    lapack2
    lapack1
    blas
    F77
    ${M_LIBRARY})

ELSE(BUILD_BIOMESH3D)

  #SET(CMAKE_CXX_FLAGS -O3)
  #SET(CMAKE_CXX_FLAGS -g)

  # Make sure the compiler can find include files from our Hellolibrary.
  INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES} ${NRRDIO_INCLUDE})

  # Make sure the linker can find the Hello library once it is built.
  LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib ${NRRDIO_LIB})

  ADD_EXECUTABLE(linear_test linear_test.cxx)
  TARGET_LINK_LIBRARIES(linear_test matrix util lapack2 lapack1 
  #lapack 
  blas F77 m)

  ADD_EXECUTABLE(linear_exp linear_exp.cxx)
  TARGET_LINK_LIBRARIES(linear_exp matrix util lapack2 lapack1 
  #lapack 
  blas F77 m)

  ADD_EXECUTABLE(linear_run linear_run.cxx)
  TARGET_LINK_LIBRARIES(linear_run matrix util lapack2 lapack1 
  #lapack 
  blas F77 m)

  ADD_EXECUTABLE(linear_labels linear_labels.cxx)
  TARGET_LINK_LIBRARIES(linear_labels matrix util lapack2 lapack1 
  #lapack 
  blas F77 m)

  ADD_EXECUTABLE(linear_labels2 linear_labels2.cxx)
  TARGET_LINK_LIBRARIES(linear_labels2 matrix util lapack2 lapack1 
  #lapack 
  blas F77 m)

  ADD_EXECUTABLE(linear_labels3 linear_labels3.cxx)
  TARGET_LINK_LIBRARIES(linear_labels3 matrix util lapack2 lapack1 
  #lapack 
  blas F77 m)

  ADD_EXECUTABLE(linear_labels4 linear_labels4.cxx)
  TARGET_LINK_LIBRARIES(linear_labels4 matrix util lapack2 lapack1 
  #lapack 
  blas F77 m)

  ADD_EXECUTABLE(linear_labels_restrict linear_labels_restrict.cxx)
  TARGET_LINK_LIBRARIES(linear_labels_restrict matrix util lapack2 lapack1 
  #lapack 
  blas F77 m)

ENDIF(BUILD_BIOMESH3D)
