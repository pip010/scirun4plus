ADD_LIBRARY(jpeg
jcapimin.c
jcapistd.c
jccoefct.c
jccolor.c
jcdctmgr.c
jchuff.c
jcinit.c
jcmainct.c
jcmarker.c
jcmaster.c
jcomapi.c
jcparam.c
jcphuff.c
jcprepct.c
jcsample.c
jctrans.c
jdapimin.c
jdapistd.c
jdatadst.c
jdatasrc.c
jdcoefct.c
jdcolor.c
jddctmgr.c
jdhuff.c
jdinput.c
jdmainct.c
jdmarker.c
jdmaster.c
jdmerge.c
jdphuff.c
jdpostct.c
jdsample.c
jdtrans.c
jerror.c
jfdctflt.c
jfdctfst.c
jfdctint.c
jidctflt.c
jidctfst.c
jidctint.c
jidctred.c
jquant1.c
jquant2.c
jutils.c
jmemmgr.c
jmemansi.c
)

#jmemnobs.c
#jmemname.c
#jmemdos.c

IF(WIN32)
  SET_TARGET_PROPERTIES(jpeg
                        PROPERTIES
                          COMPILE_FLAGS "-DBUILD_jpeg"
                          DEFINE_SYMBOL _USRDLL
                       )
ENDIF(WIN32)

SET(JPEG_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "JPEG Include." FORCE)
SET(JPEG_LIBRARY jpeg CACHE INTERNAL "JPEG Library name." FORCE)


#add_custom_target(jpeg_test ALL)
#add_custom_command(TARGET jpeg_test COMMAND echo srcdir: ${JPEG_INCLUDE})
