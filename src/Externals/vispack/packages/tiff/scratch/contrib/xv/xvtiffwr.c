/*
 * xvtiffwr.c - write routine for TIFF pictures
 *
 * WriteTIFF(fp,pic,w,h,r,g,b,numcols,style,raw,fname)
 */
#include "xv.h"
#include "tiffio.h"

static void setupColormap(tif, r, g, b)
  TIFF *tif;
  byte *r, *g, *b;
{
  short red[256], green[256], blue[256];
  int i;

  /* convert 8-bit colormap to 16-bit */
  for (i=0; i<256; i++) {
#define	SCALE(x)	(((x)*((1L<<16)-1))/255)
      red[i] = SCALE(r[i]);
      green[i] = SCALE(g[i]);
      blue[i] = SCALE(b[i]);
  }
  TIFFSetField(tif, TIFFTAG_COLORMAP, red, green, blue);
}

/*******************************************/
static int WriteTIFF(fp,pic,w,h,rmap,gmap,bmap,numcols,colorstyle,fname,comp)
FILE *fp;
byte *pic;
int   w,h;
byte *rmap, *gmap, *bmap;
int   numcols, colorstyle;
char *fname;
{
  TIFF *tif;
  byte *pix;
  int   i,j;

  tif = TIFFOpen(fname, "w");
  if (!tif)
	return 0;
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, comp);
  if (comp == COMPRESSION_CCITTFAX3)
      TIFFSetField(tif, TIFFTAG_GROUP3OPTIONS,
	  GROUP3OPT_2DENCODING+GROUP3OPT_FILLBITS);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, h);

  /* write the image data */

  if (colorstyle==0) {                  /* 8bit Palette RGB */
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_PALETTE);
    setupColormap(tif, rmap, gmap, bmap);
    TIFFWriteEncodedStrip(tif, 0, pic, w*h);
  } else if (colorstyle==1) {             /* 8-bit greyscale */
    byte rgb[256];
    byte *tpic = (byte *)malloc(w*h);
    byte *tp = tpic;
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    for (i=0; i<numcols; i++) rgb[i] = MONO(rmap[i],gmap[i],bmap[i]);
    for (i=0, pix=pic; i<w*h; i++,pix++) {
      if ((i&0x7fff)==0) WaitCursor();
      *tp++ = rgb[*pix];
    }
    TIFFWriteEncodedStrip(tif, 0, tpic, w*h);
    free(tpic);
  } else if (colorstyle==2) {             /* 1-bit B/W stipple */
    int bit,k;
    byte *tpic, *tp;
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 1);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    tpic = (byte *)malloc(TIFFStripSize(tif));
    tp = tpic;
    for (i=0, pix=pic; i<h; i++) {
      if ((i&15)==0) WaitCursor();
      for (j=0, bit=0, k=0; j<w; j++, pix++) {
	k = (k << 1) | *pix;
	bit++;
	if (bit==8) {
	  *tp++ = ~k;
	  bit = k = 0;
	}
      } /* j */
      if (bit) {
	k = k << (8-bit);
	*tp++ = ~k;
      }
    }
    TIFFWriteEncodedStrip(tif, 0, tpic, TIFFStripSize(tif));
    free(tpic);
  }
  TIFFClose(tif);

  return 0;
}

/**** Stuff for TIFFDialog box ****/

#define TWIDE 280
#define THIGH 160
#define T_NBUTTS 2
#define T_BOK    0
#define T_BCANC  1
#define BUTTH    24

#ifdef __STDC__
static void drawTD(int, int, int, int);
static void clickTD(int, int);
static void doCmd(int);
static void writeTIFF(void);
#else
static void drawTD(), doCmd(), clickTD(), writeTIFF();
#endif


/* local variables */
static char *filename;
static int   colorType;
static BUTT  tbut[T_NBUTTS];
static RBUTT *compRB;


/***************************************************/
void CreateTIFFW()
{
  CARD32     data[2];
  Atom       prop;
  int	     y;

  tiffW = CreateWindow("xv tiff", NULL, TWIDE, THIGH, infofg, infobg);
  if (!tiffW) FatalError("can't create tiff window!");

  data[0] = (CARD32) XInternAtom(theDisp, "WM_DELETE_WINDOW", FALSE);
  data[1] = (CARD32) time((long *)0);
  prop = XInternAtom(theDisp, "WM_PROTOCOLS", FALSE),
  XChangeProperty(theDisp, tiffW, prop, prop,
		  32, PropModeReplace, (unsigned char *) data, 2);

  XSelectInput(theDisp, tiffW, ExposureMask | ButtonPressMask | KeyPressMask);

  BTCreate(&tbut[T_BOK], tiffW, TWIDE-140-1, THIGH-10-BUTTH-1, 60, BUTTH, 
	   "Ok", infofg, infobg);

  BTCreate(&tbut[T_BCANC], tiffW, TWIDE-70-1, THIGH-10-BUTTH-1, 60, BUTTH, 
	   "Cancel", infofg, infobg);

  y = 55;
  compRB = RBCreate(NULL, tiffW, 36, y,   "None", infofg, infobg);
  RBCreate(compRB, tiffW, 36, y+18,       "LZW", infofg, infobg);
  RBCreate(compRB, tiffW, 36, y+36,       "PackBits", infofg, infobg);
  RBCreate(compRB, tiffW, TWIDE/2, y,     "CCITT Group3", infofg, infobg);
  RBCreate(compRB, tiffW, TWIDE/2, y+18,  "CCITT Group4", infofg, infobg);
  RBCreate(compRB, tiffW, TWIDE/2, y+36,  "JPEG", infofg, infobg);

  XMapSubwindows(theDisp, tiffW);
}
  

/***************************************************/
void TIFFDialog(vis)
int vis;
{
  if (vis) {
    CenterMapWindow(tiffW, tbut[T_BOK].x + tbut[T_BOK].w/2,
		    tbut[T_BOK].y + tbut[T_BOK].h/2, TWIDE, THIGH);
  }
  else     XUnmapWindow(theDisp, tiffW);
  tiffUp = vis;
}


/***************************************************/
int TIFFCheckEvent(xev)
XEvent *xev;
{
  /* check event to see if it's for one of our subwindows.  If it is,
     deal accordingly, and return '1'.  Otherwise, return '0' */

  int rv;
  rv = 1;

  if (!tiffUp) return 0;

  if (xev->type == Expose) {
    int x,y,w,h;
    XExposeEvent *e = (XExposeEvent *) xev;
    x = e->x;  y = e->y;  w = e->width;  h = e->height;

    if (e->window == tiffW)       drawTD(x, y, w, h);
    else rv = 0;
  }

  else if (xev->type == ButtonPress) {
    XButtonEvent *e = (XButtonEvent *) xev;
    int x,y;
    x = e->x;  y = e->y;

    if (e->button == Button1) {
      if      (e->window == tiffW)     clickTD(x,y);
      else rv = 0;
    }  /* button1 */
    else rv = 0;
  }  /* button press */


  else if (xev->type == KeyPress) {
    XKeyEvent *e = (XKeyEvent *) xev;
    char buf[128];  KeySym ks;  XComposeStatus status;  
    int stlen;
	
    stlen = XLookupString(e,buf,128,&ks,&status);
    buf[stlen] = '\0';

    if (e->window == tiffW) {
      if (stlen) {
	if (buf[0] == '\r' || buf[0] == '\n') { /* enter */
	  FakeButtonPress(&tbut[T_BOK]);
	}
	else if (buf[0] == '\033') {            /* ESC */
	  FakeButtonPress(&tbut[T_BCANC]);
	}
      }
    }
    else rv = 0;
  }
  else rv = 0;

  if (rv==0 && (xev->type == ButtonPress || xev->type == KeyPress)) {
    XBell(theDisp, 50);
    rv = 1;   /* eat it */
  }

  return rv;
}


/***************************************************/
void TIFFSaveParams(fname, col)
char *fname;
int col;
{
  filename = fname;
  colorType = col;
  if (colorType == F_BWDITHER) {
      RBSetActive(compRB,3,1);
      RBSetActive(compRB,4,1);
      RBSetActive(compRB,5,0);
  } else {
      RBSetActive(compRB,3,0);
      RBSetActive(compRB,4,0);
      RBSetActive(compRB,5,1);
  }
}


/***************************************************/
static void drawTD(x,y,w,h)
int x,y,w,h;
{
  char *title  = "Save TIFF file...";
  int  i;
  XRectangle xr;

  xr.x = x;  xr.y = y;  xr.width = w;  xr.height = h;
  XSetClipRectangles(theDisp, theGC, 0,0, &xr, 1, Unsorted);

  XSetForeground(theDisp, theGC, infofg);
  XSetBackground(theDisp, theGC, infobg);

  for (i=0; i<T_NBUTTS; i++) BTRedraw(&tbut[i]);

  ULineString(tiffW, "Compression", compRB->x-16, compRB->y-3-DESCENT);
  RBRedraw(compRB, -1);

  XDrawString(theDisp, tiffW, theGC, 20, 29, title, strlen(title));

  XSetClipMask(theDisp, theGC, None);
}


/***************************************************/
static void clickTD(x,y)
int x,y;
{
  int i;
  BUTT *bp;

  /* check BUTTs */

  /* check the RBUTTS first, since they don't DO anything */
  if ( (i=RBClick(compRB, x,y)) >= 0) { 
    (void) RBTrack(compRB, i);
    return;
  }


  for (i=0; i<T_NBUTTS; i++) {
    bp = &tbut[i];
    if (PTINRECT(x, y, bp->x, bp->y, bp->w, bp->h)) break;
  }

  if (i<T_NBUTTS) {  /* found one */
    if (BTTrack(bp)) doCmd(i);
  }
}



/***************************************************/
static void doCmd(cmd)
int cmd;
{
  switch (cmd) {
  case T_BOK:
	writeTIFF();
  case T_BCANC:
	TIFFDialog(0);
	break;
  default:
	break;
  }
}


/*******************************************/
static void writeTIFF()
{
  FILE *fp;
  int   w, h, nc, rv, comp;
  byte *rmap, *gmap, *bmap;
  byte *bwpic;

  nc = numcols;

  /* see if we can open the output file before proceeding */
  fp = OpenOutFile(filename);
  if (!fp) return;

  WaitCursor();
  
  bwpic = HandleBWandReduced(colorType, &nc, &rmap, &gmap, &bmap);

  if (colorType == F_REDUCED) colorType = F_FULLCOLOR;

  switch (RBWhich(compRB)) {
  case 0: comp = COMPRESSION_NONE; break;
  case 1: comp = COMPRESSION_LZW; break;
  case 2: comp = COMPRESSION_PACKBITS; break;
  case 3: comp = COMPRESSION_CCITTFAX3; break;
  case 4: comp = COMPRESSION_CCITTFAX4; break;
  case 5: comp = COMPRESSION_JPEG; break;
  default: comp = COMPRESSION_NONE; break;
  }
  rv = WriteTIFF(fp, (bwpic ? bwpic : epic), eWIDE, eHIGH,
      rmap, gmap, bmap, nc, colorType, filename, comp);

  SetCursors(-1);

  if (CloseOutFile(fp, filename, rv) == 0) {
    /* everything's cool! */
    DirBox(0);
    LoadCurrentDirectory();
  }
  if (bwpic) free(bwpic);

}
