// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: defs.h,v 1.1.1.1 2003/02/12 16:51:54 whitaker Exp $

#ifndef	UTIL_DEFS_H
#define	UTIL_DEFS_H

#include <stddef.h>
#include <stdlib.h>
//#include <sys/iostream>
#include <iostream>
using namespace std;
#include <memory.h>


#if (defined(_SVR4) ||\
     defined(_SVR4_SOURCE) ||\
     defined(SYSTYPE_SVR4)) && !defined(SVR4)
#define SVR4
#endif

#if defined(SVR4) && !defined(__sgi)


#	include	<sys/types.h>
#	define	boolean	bool
typedef	unsigned	bool;
#	undef	PI
#	define	PI	M_PI
#else
// Standard boolean declaration.  Shouldn't this be an unsigned char
// to save space?
typedef unsigned 	boolean;
#endif

// ---------------------------------------------------------------------------
// Comon Grasp types

typedef	float		VISReal;	// Coordinates are 32bit "C" floats
typedef unsigned char	byte;		// bytes...

// ---------------------------------------------------------------------------

#ifndef True
#define True 1
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef False
#define False 0
#endif
#ifndef FALSE
#define FALSE 0
#endif

// ---------------------------------------------------------------------------

#endif
