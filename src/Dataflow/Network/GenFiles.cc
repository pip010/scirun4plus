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


/* GenFiles.cc */

#include <stdio.h>
#include <sys/stat.h>
#include <string>
#include <string.h>
#include <stdio.h>

#include <Dataflow/Network/ComponentNode.h>
#include <Dataflow/Network/SkeletonFiles.h>
#include <Dataflow/Network/GenFiles.h>
#include <Core/Util/FileUtils.h>
#include <Core/Util/Environment.h>
#include <Core/Util/Dir.h> // for MKDIR
#include <sci_debug.h>

#define PERM S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH

namespace SCIRun {

using std::string;

/*! these functions all assume that the files and/or directories
    about to be generated, do not already exist */ 

int
GenPackage(char* package, char* psepath)
{
#if DEBUG
  printf ("Begin GenPack\n");
#endif
  int check=0,checkall=0;
  char* strbuf = 0;
  FILE* file = 0;

  bool pse_core_dir = true;
  std::string packstring("");
  if (strcmp(package,"SCIRun"))
    {
      pse_core_dir = false;
      packstring = std::string("Packages/") + package + std::string("/");
    }
  const char *packdir = packstring.c_str();


  /* allocate enough space to hold the largest path */
  strbuf = new char[strlen(packdir)+strlen(psepath)+50];

  /* create all directories associated with a package */
  sprintf(strbuf,"%s/%s",psepath,packdir);
  // TODO: change to use Dir.cc
  checkall |= check = MKDIR(strbuf,PERM);
  if (check)
    printf("could not create directory \"%s\"\n",strbuf);

  sprintf(strbuf,"%s/%sDataflow",psepath,packdir);
  // TODO: change to use Dir.cc
  checkall |= check = MKDIR(strbuf,PERM);
  if (check)
    printf("could not create directory \"%s\"\n",strbuf);

  sprintf(strbuf,"%s/%sCore",psepath,packdir);
  // TODO: change to use Dir.cc
  checkall |= check = MKDIR(strbuf,PERM);
  if (check)
    printf("could not create directory \"%s\"\n",strbuf);

  sprintf(strbuf,"%s/%sDataflow/Modules",psepath,packdir);
  checkall |= check = MKDIR(strbuf,PERM);
  if (check)
    printf("could not create directory \"%s\"\n",strbuf);
  
  sprintf(strbuf,"%s/%sDataflow/GUI",psepath,packdir);
  checkall |= check = MKDIR(strbuf,PERM);
  if (check)
    printf("could not create directory \"%s\"\n",strbuf);

  sprintf(strbuf,"%s/%sDataflow/XML",psepath,packdir);
  checkall |= check = MKDIR(strbuf,PERM);
  if (check)
    printf("could not create directory \"%s\"\n",strbuf);

#if 0
  sprintf(strbuf,"%s/%sshare",psepath,packdir);
  checkall |= check = MKDIR(strbuf,PERM);
  if (check)
    printf("could not create directory \"%s\"\n",strbuf);
#endif

  sprintf(strbuf,"%s/%sCore/Datatypes",psepath,packdir);
  checkall |= check = MKDIR(strbuf,PERM);
  if (check)
    printf("could not create directory \"%s\"\n",strbuf);

  if (checkall) {
    printf("Could not create one or more directories.  Giving up.");
    return 0;
  }

#if 0
  /* create all the non-directory files associated with a package */
  sprintf(strbuf,"%s/%sshare/share.h",psepath,package);
  file = fopen(strbuf,"w");
  fprintf(file,share_skeleton,package,package,package,
	  package,package);
  fclose(file);

  sprintf(strbuf,"%s/%sshare/DllEntry.cc",psepath,packdir);
  file = fopen(strbuf,"w");
  fprintf(file,dllentry_skeleton,package);
  fclose(file);
#endif

  if (!pse_core_dir)
    {
      sprintf( strbuf, "%s/%sCMakeLists.txt", psepath, packdir );
      file = fopen(strbuf,"w");
      fprintf( file, package_cmakelists_skeleton, package );
      fclose(file);
    }

  sprintf(strbuf,"%s/%sDataflow/CMakeLists.txt", psepath,packdir);
  file = fopen(strbuf,"w");
  fprintf( file, dataflow_cmakelists_skeleton, packdir );
  fclose(file);

  sprintf( strbuf, "%s/%sCore/CMakeLists.txt", psepath, packdir );
  file = fopen(strbuf,"w");
  fprintf( file, core_cmakelists_skeleton, packdir );
  fclose(file);

  sprintf( strbuf, "%s/%sDataflow/Modules/CMakeLists.txt", psepath, packdir );
  file = fopen(strbuf,"w");
  fprintf( file, modules_cmakelists_skeleton, packdir );
  fclose(file);

  sprintf( strbuf, "%s/%sCore/Datatypes/CMakeLists.txt", psepath, packdir );
  file = fopen(strbuf,"w");
  fprintf( file, datatypes_cmakelists_skeleton, package, package, package );
  fclose(file);

#if 0 
  // GUI dirs don't have CMakeLists.txt file... at least not yet.
  sprintf(strbuf,"%s/%sDataflow/GUI/CMakeLists.txt",psepath,packdir);
  file = fopen(strbuf,"w");
  fprintf(file,gui_cmakelists_skeleton,packdir,packdir);
  fclose(file);
#endif

  delete[] strbuf;

#if DEBUG
  printf ("End GenPack\n");
#endif

  return 1;
}

int
GenCategory(char* catname, char* package, char* psepath)
{
#if DEBUG
  printf ("begin GenCat\n");
#endif
  int check;
  char* strbuf = 0;
  FILE* file = 0;

  std::string packstring("");
  if (strcmp(package,"SCIRun"))
    packstring = std::string("Packages/") + package + std::string("/");
  const char *packdir = packstring.c_str();


  strbuf = new char[strlen(packdir)+strlen(psepath)+
		    strlen(catname)+50];

  /* create category directory */
  sprintf(strbuf,"%s/%sDataflow/Modules/%s",
	  psepath,packdir,catname);
  check = MKDIR(strbuf,PERM);
  if (check) {
    printf("could not create directory \"%s\".  Giving up.\n",strbuf);
    return 0;
  }

  // Create category CMakeLists.txt file...
  sprintf( strbuf, "%s/%sDataflow/Modules/%s/CMakeLists.txt", psepath, packdir, catname );
  file = fopen(strbuf,"w");
  fprintf( file, category_cmakelists_skeleton, packdir, catname);
  fclose(file);

  // Edit the modules CMakeLists.txt file - add the new category.

  char* modname = new char[strlen(catname)+50];
  sprintf( modname,"  %s\n", catname );
  strbuf = new char[strlen(psepath)+strlen(packdir)+50];
  sprintf( strbuf, "%s/%sDataflow/Modules/CMakeLists.txt", psepath, packdir );
  InsertStringInFile(strbuf,"#[INSERT NEW CATEGORY DIR HERE]",modname);
  delete[] strbuf;

#if DEBUG
  printf ("end GenCat\n");
#endif
  return 1;
}

int
GenComponent(const ModuleInfo &mi, char* package, char* psepath)
{
#if DEBUG
  printf ("Begin GenComp\n");
#endif
  
  char* filename = 0;
  char* strbuf = 0;
  FILE* file = 0;
  int length;

  std::string packstring("");
  if (strcmp(package,"SCIRun"))
    packstring = std::string("Packages/") + package + std::string("/");
  const char *packdir = packstring.c_str();

  /* generate a skeleton .cc file */
  length = mi.module_name_.size() + strlen(psepath)+
    strlen(packdir) + mi.category_name_.size() + 50;
  filename = new char[length];
  sprintf(filename,"%s/%sDataflow/Modules/%s/%s.cc",
	  psepath, packdir, 
	  mi.category_name_.c_str(), mi.module_name_.c_str());
  file = fopen(filename,"w");

  if( file == NULL ) {
    printf("Error, fopen failed for file %s!\n", filename);
    return 0;
  }

  const char* mname = mi.module_name_.c_str();
  const char* cname = mi.category_name_.c_str();
  fprintf(file,component_skeleton,mname,
	  sci_getenv("USER"),"TODAY'S DATE HERE",
	  package,
	  mname,mname,mname,
	  mname,mname,mname,mname,
	  cname,package,mname,mname,
	  mname,mname,package);
  fclose(file);
  delete[] filename;

  /* generate a full component .xml file */
  length = strlen(mname)+strlen(psepath)+
    strlen(packdir)+50;
  filename = new char[length];
  sprintf(filename,"%s/%sDataflow/XML/%s.xml",psepath,
	  packdir,mname);
  write_component_file(mi, filename);
  delete[] filename;

  if (mi.has_gui_node_) {
    /* generate a skeleton .tcl file */
    length = strlen(mname)+strlen(psepath)+
      strlen(packdir)+50;
    filename = new char[length];
    sprintf(filename,"%s/%sDataflow/GUI/%s.tcl",psepath,
	    packdir,mname);
    file = fopen(filename,"w");
    fprintf(file,gui_skeleton,package,cname,mname,mname,filename);
    fclose(file);
    delete[] filename;
  }

  // Edit the category CMakeLists.txt file - add the new component:
  char* modname = new char[strlen(mname)+50];
  sprintf(modname,"\n  %s.cc",mname);
  strbuf = new char[strlen(psepath)+strlen(packdir)+strlen(cname)+50];
  sprintf( strbuf, "%s/%sDataflow/Modules/%s/CMakeLists.txt", psepath, packdir, cname);

  // Setup the MATCH string to find in the CMakeLists.txt file:
  //
  std::string package_directory = "";
  if( strlen(packdir) > 0 ) {
    std::string temp = packdir;
    std::string::size_type loc1 = temp.find( "/" );
    std::string::size_type loc2 = temp.rfind( "/" );
    if( loc1 == std::string::npos || loc2 == std::string::npos ) {
      printf( "Error: GenFiles.cc: GenComponent():  packdir (%s) is formatted badly.\n", packdir );
      return 0;
    }
    package_directory = temp.substr( 0, loc1 ) + "_" + temp.substr( loc1+1, loc2-loc1-1 ) + "_";
  }
  std::string match_string = "SET(" + package_directory + "Dataflow_Modules_" + cname + "_SRCS";

  InsertStringInFile( strbuf, match_string.c_str(), modname );
  delete[] strbuf;
  delete[] modname;

#if 0
  // There is no GUI/sub.mk now... now sure what, if anything, we should do here...
  if (mi.has_gui_node_) {
    /* edit the GUI sub.mk file - add the new component */
    modname = new char[strlen(mname)+50];
    sprintf(modname,"\t$(SRCDIR)/%s.tcl\\\n",mname);
    strbuf = new char[strlen(psepath)+strlen(packdir)+strlen(cname)+50];
    sprintf(strbuf,"%s/%sDataflow/GUI/sub.mk",psepath,packdir);
    InsertStringInFile(strbuf,"#[INSERT NEW TCL FILE HERE]",modname);
    delete[] strbuf;
  }
#endif
#if DEBUG
  printf ("End GenComp\n");
#endif
  return 1;
} // end GenComponent

} // End namespace SCIRun



