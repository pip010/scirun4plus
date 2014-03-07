SET(__i386__ TRUE)
ADD_LIBRARY(fitsio
fitscore.c
buffers.c
cfileio.c
checksum.c
getcol.c
getcolb.c
getcold.c
getcole.c
getcoli.c
getcolj.c
getcolk.c
getcoll.c
getcols.c
getcolui.c
getcoluj.c
editcol.c
edithdu.c
getkey.c
modkey.c
putcol.c
putcolb.c
putcold.c
putcole.c
putcoli.c
putcolj.c
putcolk.c
putcoll.c
putcols.c
putcolui.c
putcolu.c
putcoluj.c
putkey.c
scalnull.c
swapproc.c
wcsutil.c
)

IF(WIN32)
  SET_TARGET_PROPERTIES(fitsio
                        PROPERTIES
                          COMPILE_FLAGS "-DBUILD_fitsio"
                          DEFINE_SYMBOL _USRDLL
                       )
ENDIF(WIN32)
