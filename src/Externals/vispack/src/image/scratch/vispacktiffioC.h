/* contains items that are common to both library and interface */

#ifndef tiffio_c_h
#define tiffio_c_h

#define CONV24_8BIT  0
#define CONV24_24BIT 1
#define PIC8  CONV24_8BIT
#define PIC24 CONV24_24BIT
#define PIC16  3
#define F_TIFF      9

typedef struct 
{ 
    byte *pic;                  /* image data */
    int   w, h;                 /* size */
    int   type;                 /* PIC8 or PIC24 */

    byte  r[256],g[256],b[256];
				/* colormap, if PIC8 */

    int   frmType;              /* def. Format type to save in */
    int   colType;              /* def. Color type to save in */
    char  fullInfo[128];        /* Format: field in info box */
    char  shrtInfo[128];        /* short format info */
    char *comment;              /* comment text */
		 
    int   numpages;             /* # of page files, if >1 */
    char  pagebname[64];        /* basename of page files */
} PICINFO;


void destroyPICINFO(PICINFO *struct_ptr)
{
    if (struct_ptr->comment != NULL)
	{
	    free(struct_ptr->comment);
	    //printf("got destroy pic info\n");
	}
    else
	{
	    // printf("width %d and height %d\n", struct_ptr->w, struct_ptr->h);
	}
    struct_ptr->comment = NULL;
    if (struct_ptr->pic != NULL)
	free(struct_ptr->pic);
    struct_ptr->pic = NULL;
}

/* these are color styles */

#define F_FULLCOLOR 0
#define F_GREYSCALE 1
#define F_BWDITHER  2
#define F_REDUCED   3
#define F_COLORINDEXED 4
#define F_RGB 5

#endif

