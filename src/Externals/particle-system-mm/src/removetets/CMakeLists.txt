
ADD_EXECUTABLE( removetets 
		main.cpp 
		ScalarField.cxx
		IOScalarField.cxx
		Surface.cxx
		SurfaceParameters.cxx
		mtxlib.cxx
		)

ADD_EXECUTABLE( removetetsall 
		main_all.cpp 
		ScalarField.cxx
		IOScalarField.cxx
		Surface.cxx
		SurfaceParameters.cxx
		mtxlib.cxx
		)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}  ${PARTICLEUTIL_INCLUDES})
TARGET_LINK_LIBRARIES(removetetsall particleutil)

TARGET_LINK_LIBRARIES(removetets particleutil)