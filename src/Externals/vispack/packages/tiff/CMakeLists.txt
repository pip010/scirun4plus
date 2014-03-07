ADD_DEFINITIONS(
  -DHAVE_IEEEFP=1
  -DPACKBITS_SUPPORT
  -DLZW_SUPPORT 
  -DCCITT_SUPPORT 
  -DTHUNDER_SUPPORT  
  -DNEXT_SUPPORT  
  -DUSE_VARARGS=0 
  -DUSE_PROTOTYPES=1 
  -DCOLORIMETRY_SUPPORT 
  -DJPEG_SUPPORT 
  -DYCBCR_SUPPORT 
  -DCMYK_SUPPORT 
  -DBSDTYPES 
)

IF(UNIX)
  ADD_DEFINITIONS(
    -Dunix
    -DMMAP_SUPPORT
  ) 
ENDIF(UNIX)

ADD_LIBRARY(tiff
  tif_fax3.c
  tif_fax4.c
  tif_aux.c
  tif_ccittrle.c
  tif_close.c
  tif_compat.c 
  tif_compress.c
  tif_dir.c
  tif_dirinfo.c
  tif_dirread.c
  tif_dirwrite.c
  tif_dumpmode.c
  tif_error.c
  tif_getimage.c
  tif_jpeg.c
  tif_flush.c
  tif_lzw.c
  tif_next.c
  tif_open.c
  tif_packbits.c
  tif_print.c
  tif_read.c
  tif_swab.c
  tif_strip.c
  tif_thunder.c
  tif_tile.c
  #tif_unix.c
  tif_version.c
  tif_warning.c
  tif_write.c
)


IF(WIN32)
  SET_TARGET_PROPERTIES(tiff
                        PROPERTIES
                          COMPILE_FLAGS "-DBUILD_tiff"
                          DEFINE_SYMBOL _USRDLL
                       )
ENDIF(WIN32)
