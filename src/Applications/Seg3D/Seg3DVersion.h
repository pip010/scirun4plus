//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : Seg3DVersion.h
//    Author : Allen Sanderson
//    Date   : Dec 2008

#ifndef SEG3D_SEG3DVERSION_H
#define SEG3D_SEG3DVERSION_H

#include <string>

namespace SCIRun {

/* Version numbers */
#define SEG3D_VERS_MAJOR    1	/* For major interface/format changes  	  */
#define SEG3D_VERS_MINOR    13	/* For minor interface/format changes  	  */
#define SEG3D_VERS_RELEASE  0	/* For tweaks, bug-fixes, or development  */
#define SEG3D_VERS_STRING "1.13.0"

unsigned int GET_MAJOR_VERSION( string version );

unsigned int GET_MINOR_VERSION( string version );

unsigned int GET_RELEASE_VERSION( string version );

int check_version( unsigned int major,
		   unsigned int minor,
		   unsigned int release );

} // end namespace SCIRun

#endif // SEG3D_SEG3DVERSION_H
