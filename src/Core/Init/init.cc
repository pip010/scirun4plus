/*
  For more information, please see: http://software.sci.utah.edu

  The MIT License

  Copyright (c) 2009 Scientific Computing and Imaging Institute,
  University of Utah.

  
  Permission is hereby granted, free of charge, to any person obtaining a
  copy of this software and associated documentation files (the "Software"),
  to deal in the Software without restriction, including without limitation
  the rights to use, copy, modify, merge, publish, distribute, sublicense,
  and/or sell copies of the Software, and to permit persons to whom the
  Software is furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included
  in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
  DEALINGS IN THE SOFTWARE.
*/


/*
 *  init.cc: 
 *
 *  Written by:
 *   McKay Davis
 *   Scientific Computing and Imaging Institute
 *   University of Utah
 *   Mar 2005
 *
 */

#include <Core/Init/init.h>
#include <Core/Util/soloader.h>
#include <Core/Util/Environment.h>
#include <Core/Util/StringUtil.h>
#include <Core/ImportExport/Matrix/MatrixIEPlugin.h>
#include <Core/ImportExport/Nrrd/NrrdIEPlugin.h>
#include <Core/ImportExport/Field/FieldIEPlugin.h>
#include <Core/Matlab/matlabfile.h>
#include <Core/Matlab/matlabarray.h>
#include <Core/Matlab/matlabconverter.h>
#include <Core/Datatypes/MatrixTypeConverter.h>



using namespace SCIRun;
using std::string;

#ifdef __APPLE__
  static string lib_ext = ".dylib";
  static string lib_base = "lib";
#elif defined(_WIN32)
  static string lib_ext = ".dll";
  static string lib_base = "";
#else
  static string lib_ext = ".so";
  static string lib_base = "lib";
#endif


// SCIRunInit is called from from main() (and from stantalone converters).
//
// This method calls ${PACKAGE}Init(0) for every package in the comma-seperated
// 'packages' string.  Called for every package in SCIRUN_LOAD_PACKAGE
// if packages is empty (which is the default)
//
// Note, the void * ${PACKAGE}Init(void *) function must be declared extern C
// in the SCIRun/src/Packages/${PACKAGE}/Core/Datatypes directory
//
// For example, for the BioPSE package:
//
// extern "C" void * BioPSEInit(void *param) 
// {
//   std::cerr << "BioPSEInit called.\n";
//   return 0;
// }

void
SCIRunInit(string packages) 
{
  string package_list = packages;

  // Force package to load
  if( package_list.empty() ) 
  {
    const char * env_value = sci_getenv("SCIRUN_LOAD_PACKAGE");
    if( env_value ) {
      package_list = env_value;
    }
  }

  std::vector<string> package = split_string(package_list, ',');
  typedef void *(*PackageInitFunc)(void *);

  for (unsigned int i = 0; i < package.size(); ++i) 
  {
    std::string libName = lib_base + "Packages_" + package[i] +
                          "_Core_Datatypes" + lib_ext;
    LIBRARY_HANDLE lib = findLib(libName);
    PackageInitFunc init = 0;
    if (lib)
    {
      std::string error;
      init = reinterpret_cast<PackageInitFunc>
                             (GetHandleSymbolAddress(lib, (package[i]+"Init"),
                                                     error));
      if(init == NULL) {
        std::cerr << "Loading package initialization function '"
                  << package[i] + "Init" << "' for package '" << libName << "'"
                  << " failed.  Error: " << error << std::endl;
      }
    }
    if (init) 
    {
      (*init)(0);
    }
  }
  
  // Force library to load
  LIBRARY_HANDLE plugin = findLib(lib_base+"Core_IEPlugin"+lib_ext);
  if (plugin == 0)
  {
    std::cerr << "Could not load plugins for reading and writing datafiles"<<std::endl;
  }
  
}


