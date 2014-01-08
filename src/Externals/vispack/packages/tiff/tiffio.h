/* $Header: /uusoc/res/image/CVS/vispack/packages/tiff/tiffio.h,v 1.1.1.1 2003/02/12 16:51:50 whitaker Exp $ */

/*
 * Copyright (c) 1988, 1989, 1990, 1991, 1992 Sam Leffler
 * Copyright (c) 1991, 1992 Silicon Graphics, Inc.
 *
 * Permission to use, copy, modify, distribute, and sell this software and 
 * its documentation for any purpose is hereby granted without fee, provided
 * that (i) the above copyright notices and this permission notice appear in
 * all copies of the software and related documentation, and (ii) the names of
 * Sam Leffler and Silicon Graphics may not be used in any advertising or
 * publicity relating to the software without the specific, prior written
 * permission of Sam Leffler and Silicon Graphics.
 * 
 * THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY 
 * WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 * IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR
 * ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
 * OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
 * WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF 
 * LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE 
 * OF THIS SOFTWARE.
 */

#ifndef _TIFFIO_
#define	_TIFFIO_

#include "tiffExports.h"

/*
 * TIFF I/O Library Definitions.
 */
#include "tiff.h"

/*
 * TIFF is defined as an incomplete type to hide the
 * library's internal data structures from clients.
 */
typedef	struct tiff TIFF;

#ifndef NULL
#define	NULL	0
#endif

/*
 * Flags to pass to TIFFPrintDirectory to control
 * printing of data structures that are potentially
 * very large.   Bit-or these flags to enable printing
 * multiple items.
 */
#define	TIFFPRINT_NONE		0x0		/* no extra info */
#define	TIFFPRINT_STRIPS	0x1		/* strips/tiles info */
#define	TIFFPRINT_CURVES	0x2		/* color/gray response curves */
#define	TIFFPRINT_COLORMAP	0x4		/* colormap */
#define	TIFFPRINT_JPEGQTABLES	0x100		/* JPEG Q matrices */
#define	TIFFPRINT_JPEGACTABLES	0x200		/* JPEG AC tables */
#define	TIFFPRINT_JPEGDCTABLES	0x200		/* JPEG DC tables */

#if defined(__STDC__) || defined(__EXTENDED__) || USE_CONST
extern const char TIFFVersion[];
extern const unsigned char TIFFBitRevTable[256];
extern const unsigned char TIFFNoBitRevTable[256];
#else
extern char TIFFVersion[];
extern unsigned char TIFFBitRevTable[256];
extern unsigned char TIFFNoBitRevTable[256];
#endif

/*
 * Macros for extracting components from the
 * packed ABGR form returned by TIFFReadRGBAImage.
 */
#define	TIFFGetR(abgr)	((abgr) & 0xff)
#define	TIFFGetG(abgr)	(((abgr) >> 8) & 0xff)
#define	TIFFGetB(abgr)	(((abgr) >> 16) & 0xff)
#define	TIFFGetA(abgr)	(((abgr) >> 24) & 0xff)

#if defined(c_plusplus) || defined(__cplusplus) || defined(__STDC__) || defined(__EXTENDED__) || USE_PROTOTYPES
#include <stdio.h>
#include <stdarg.h>

typedef	void (*TIFFErrorHandler)(char* module, char* fmt, va_list ap);

#if defined(__cplusplus)
extern "C" {
#endif
tiff_SHARE extern	void TIFFClose(TIFF*);
tiff_SHARE extern	int TIFFFlush(TIFF*);
tiff_SHARE extern	int TIFFFlushData(TIFF*);
tiff_SHARE extern	int TIFFGetField(TIFF*, int, ...);
tiff_SHARE extern	int TIFFVGetField(TIFF*, int, va_list);
tiff_SHARE extern	int TIFFGetFieldDefaulted(TIFF*, int, ...);
tiff_SHARE extern	int TIFFVGetFieldDefaulted(TIFF*, int, va_list);
tiff_SHARE extern	int TIFFReadDirectory(TIFF*);
tiff_SHARE extern	int TIFFScanlineSize(TIFF*);
tiff_SHARE extern	unsigned long TIFFStripSize(TIFF*);
tiff_SHARE extern	unsigned long TIFFVStripSize(TIFF*, unsigned long);
tiff_SHARE extern	unsigned long TIFFTileRowSize(TIFF*);
tiff_SHARE extern	unsigned long TIFFTileSize(TIFF*);
tiff_SHARE extern	unsigned long TIFFVTileSize(TIFF*, unsigned long);
tiff_SHARE extern	int TIFFFileno(TIFF*);
tiff_SHARE extern	int TIFFGetMode(TIFF*);
tiff_SHARE extern	int TIFFIsTiled(TIFF*);
tiff_SHARE extern	long TIFFCurrentRow(TIFF*);
tiff_SHARE extern	int TIFFCurrentDirectory(TIFF*);
tiff_SHARE extern	int TIFFCurrentStrip(TIFF*);
tiff_SHARE extern	int TIFFCurrentTile(TIFF*);
tiff_SHARE extern	int TIFFReadBufferSetup(TIFF*, char*, unsigned);
tiff_SHARE extern	int TIFFSetDirectory(TIFF*, int);
tiff_SHARE extern	int TIFFSetField(TIFF*, int, ...);
tiff_SHARE extern	int TIFFVSetField(TIFF*, int, va_list);
tiff_SHARE extern	int TIFFWriteDirectory(TIFF *);
#if defined(c_plusplus) || defined(__cplusplus)
tiff_SHARE extern	TIFF* TIFFOpen(const char*, const char*);
tiff_SHARE extern	TIFF* TIFFFdOpen(const int, const char*, const char*);
tiff_SHARE extern	const char* TIFFFileName(TIFF*);
tiff_SHARE extern	void TIFFError(const char*, const char*, ...);
tiff_SHARE extern	void TIFFWarning(const char*, const char*, ...);
tiff_SHARE extern	void TIFFPrintDirectory(TIFF*, FILE*, long = 0);
tiff_SHARE extern	int TIFFReadScanline(TIFF*, unsigned char*, unsigned, unsigned = 0);
tiff_SHARE extern	int TIFFWriteScanline(TIFF*, unsigned char*, unsigned, unsigned = 0);
tiff_SHARE extern	int TIFFReadRGBAImage(TIFF*, unsigned long, unsigned long, unsigned long*, int stop = 0);
#else
tiff_SHARE extern	TIFF* TIFFOpen(char*, char*);
tiff_SHARE extern	TIFF* TIFFFdOpen(int, char*, char*);
tiff_SHARE extern	char* TIFFFileName(TIFF*);
tiff_SHARE extern	void TIFFError(char*, char*, ...);
tiff_SHARE extern	TIFFErrorHandler TIFFSetErrorHandler(TIFFErrorHandler handler);
tiff_SHARE extern	void TIFFWarning(char*, char*, ...);
tiff_SHARE extern	TIFFErrorHandler TIFFSetWarningHandler(TIFFErrorHandler handler);
tiff_SHARE extern	void TIFFPrintDirectory(TIFF*, FILE*, long);
tiff_SHARE extern	int TIFFReadScanline(TIFF*, unsigned char*, unsigned, unsigned);
tiff_SHARE extern	int TIFFWriteScanline(TIFF*, unsigned char*, unsigned, unsigned);
tiff_SHARE extern	int TIFFReadRGBAImage(TIFF*, unsigned long, unsigned long, unsigned long*, int stop);
#endif
tiff_SHARE extern	unsigned int TIFFComputeTile(TIFF*, unsigned long, unsigned long, unsigned int, unsigned long);
tiff_SHARE extern	int TIFFCheckTile(TIFF*, unsigned long, unsigned long, unsigned long, unsigned);
tiff_SHARE extern	unsigned int TIFFNumberOfTiles(TIFF*);
tiff_SHARE extern	int TIFFReadTile(TIFF*, unsigned char*, unsigned long, unsigned long, unsigned long, unsigned);
tiff_SHARE extern	unsigned int TIFFComputeStrip(TIFF*, unsigned long, unsigned int);
tiff_SHARE extern	unsigned int TIFFNumberOfStrips(TIFF*);
tiff_SHARE extern	int TIFFReadEncodedStrip(TIFF*, unsigned, unsigned char*, unsigned);
tiff_SHARE extern	int TIFFReadRawStrip(TIFF*, unsigned, unsigned char*, unsigned);
tiff_SHARE extern	int TIFFReadEncodedTile(TIFF*, unsigned, unsigned char*, unsigned);
tiff_SHARE extern	int TIFFReadRawTile(TIFF*, unsigned, unsigned char*, unsigned);
tiff_SHARE extern	int TIFFWriteEncodedStrip(TIFF*, unsigned, unsigned char*, unsigned);
tiff_SHARE extern	int TIFFWriteRawStrip(TIFF*, unsigned, unsigned char*, unsigned);
tiff_SHARE extern	int TIFFWriteEncodedTile(TIFF*, unsigned, unsigned char*, unsigned);
tiff_SHARE extern	int TIFFWriteRawTile(TIFF*, unsigned, unsigned char*, unsigned);
tiff_SHARE extern	int TIFFSwabShort(unsigned short *);
tiff_SHARE extern	int TIFFSwabLong(unsigned long *);
tiff_SHARE extern	int TIFFSwabVISArrayOfShort(unsigned short *, int);
tiff_SHARE extern	int TIFFSwabVISArrayOfLong(unsigned long *, int);
tiff_SHARE extern	int TIFFReverseBits(unsigned char *, int);
#if defined(__cplusplus)
}
#endif
#else
typedef	void (*TIFFErrorHandler)();

tiff_SHARE extern	void TIFFClose();
tiff_SHARE extern	TIFF *TIFFOpen();
tiff_SHARE extern	TIFF *TIFFFdOpen();
tiff_SHARE extern	char* TIFFFileName();
tiff_SHARE extern	int TIFFFileno();
tiff_SHARE extern	int TIFFGetMode();
tiff_SHARE extern	int TIFFIsTiled();
tiff_SHARE extern	unsigned int TIFFComputeTile();
tiff_SHARE extern	long TIFFCurrentRow();
tiff_SHARE extern	int TIFFCurrentDirectory();
tiff_SHARE extern	int TIFFCurrentStrip();
tiff_SHARE extern	int TIFFCurrentTile();
tiff_SHARE extern	void TIFFError();
tiff_SHARE extern	TIFFErrorHandler TIFFSetErrorHandler();
tiff_SHARE extern	int TIFFFlush();
tiff_SHARE extern	int TIFFFlushData();
tiff_SHARE extern	int TIFFGetField();
tiff_SHARE extern	int TIFFVGetField();
tiff_SHARE extern	int TIFFGetFieldDefaulted();
tiff_SHARE extern	int TIFFVGetFieldDefaulted();
tiff_SHARE extern	unsigned int TIFFNumberOfTiles();
tiff_SHARE extern	void TIFFPrintDirectory();
tiff_SHARE extern	int TIFFReadDirectory();
tiff_SHARE extern	int TIFFReadBufferSetup();
tiff_SHARE extern	int TIFFReadScanline();
tiff_SHARE extern	int TIFFReadTile();
tiff_SHARE extern	unsigned int TIFFComputeStrip();
tiff_SHARE extern	unsigned int TIFFNumberOfStrips();
tiff_SHARE extern	int TIFFReadEncodedStrip();
tiff_SHARE extern	int TIFFReadRawStrip();
tiff_SHARE extern	int TIFFReadEncodedTile();
tiff_SHARE extern	int TIFFReadRGBAImage();
tiff_SHARE extern	int TIFFReadRawTile();
tiff_SHARE extern	int TIFFScanlineSize();
tiff_SHARE extern	unsigned long TIFFStripSize();
tiff_SHARE extern	unsigned long TIFFVStripSize();
tiff_SHARE extern	unsigned long TIFFTileRowSize();
tiff_SHARE extern	unsigned long TIFFTileSize();
tiff_SHARE extern	unsigned long TIFFVTileSize();
tiff_SHARE extern	int TIFFSetDirectory();
tiff_SHARE extern	int TIFFSetField();
tiff_SHARE extern	int TIFFVSetField();
tiff_SHARE extern	void TIFFWarning();
tiff_SHARE extern	TIFFErrorHandler TIFFSetWarningHandler();
tiff_SHARE extern	int TIFFWriteDirectory();
tiff_SHARE extern	int TIFFWriteScanline();
tiff_SHARE extern	int TIFFWriteEncodedStrip();
tiff_SHARE extern	int TIFFWriteRawStrip();
tiff_SHARE extern	int TIFFWriteEncodedTile();
tiff_SHARE extern	int TIFFWriteRawTile();
tiff_SHARE extern	int TIFFSwabShort();
tiff_SHARE extern	int TIFFSwabLong();
tiff_SHARE extern	int TIFFSwabVISArrayOfShort();
tiff_SHARE extern	int TIFFSwabVISArrayOfLong();
tiff_SHARE extern	int TIFFReverseBits();
tiff_SHARE extern	int TIFFCheckTile();
#endif
#endif /* _TIFFIO_ */
