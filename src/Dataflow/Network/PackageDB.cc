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


// PackageDB.cc - Interface to module-finding and loading mechanisms

#ifdef ASSERT
#  undef ASSERT
#endif

#include <Dataflow/Network/PackageDB.h>
#include <Dataflow/Network/ComponentNode.h>
#include <Dataflow/Network/NetworkEditor.h>
#include <Core/XMLUtil/XMLUtil.h>
#include <Core/Util/StringUtil.h>
#include <Dataflow/GuiInterface/TCLInterface.h>
#include <Core/Util/Environment.h>
#include <Core/Util/soloader.h>
#include <Core/Util/FileUtils.h>
#include <Core/Util/Dir.h> // for LSTAT

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ctype.h>
#include <string>
#include <vector>
#include <map>

#include <libxml/xmlreader.h>
#include <libxml/catalog.h>

#include <sys/stat.h>

#ifdef __APPLE__
  static std::string lib_ext = ".dylib";
#elif defined(_WIN32)
  const std::string lib_ext = ".dll";
#else
  static std::string lib_ext = ".so";
#endif

namespace SCIRun {

  PackageDB* packageDB = 0;

  struct category {
    std::string name;
    std::map<std::string,ModuleInfo*> modules;
  };
  
  struct package {
    std::string name;
    std::map<std::string,category*> categories;
  };

  typedef std::map<int,char*>::iterator char_iter;
//  typedef std::map<int,inport_node*>::iterator inport_iter;
//  typedef std::map<int,outport_node*>::iterator outport_iter;
  typedef std::map<std::string,int>::iterator string_iter;
  typedef std::map<std::string,category*>::iterator category_iter;
  typedef std::map<std::string,ModuleInfo*>::iterator module_iter;
  typedef std::map<int,package*>::iterator package_iter;
}

using namespace SCIRun;
  
PackageDB::PackageDB()
  : db_(new Packages), 
    packageList_(0)
{
}

PackageDB::~PackageDB()
{ 
  delete db_; 
}

typedef void (*pkgInitter)(const std::string& tclPath);

static void log(const std::string &m)
{
#ifdef DEBUG
  std::cerr << m;
  TCLInterface::post_message(m);
  TCLInterface::execute("update idletasks");
#endif
}

bool
PackageDB::findMaker(ModuleInfo* moduleInfo)
{
  std::string cat_bname, pak_bname;
  if(moduleInfo->package_name_ == "SCIRun") {
    cat_bname = "Dataflow_Modules_";
    pak_bname = "Dataflow";
  } else {
    cat_bname = "Packages_" + moduleInfo->package_name_ + "_Dataflow_Modules_";
    pak_bname = "Packages_" + moduleInfo->package_name_ + "_Dataflow";
  }
  std::string errstr;
  std::string errorstr;

  std::string libname = pak_bname+lib_ext;
#ifndef _WIN32
  // don't add the "lib" under windows
  libname = std::string("lib") + libname;
#endif
  // try the large version of the shared library
  LIBRARY_HANDLE package_so = findLib(libname,errorstr);
#ifdef DEBUG
//  if(NULL == package_so) {
//    log("could not load lib '" + libname + "': " + errorstr + "\n");
//    errstr = std::string(" - ")+errorstr+std::string("\n");
//  }
#endif

  // If package is FieldsChoose, FieldsCreate, FieldsHandle Fields Packages,
  // or UnuA-M, UnuN-Z
  std::string cat_name = moduleInfo->category_name_;
  if((cat_name.find("Field") != std::string::npos || 
      cat_name.find("Mesh") != std::string::npos) &&
      cat_name.find("FieldArray") == std::string::npos &&
      moduleInfo->package_name_ == "SCIRun")
  { cat_name = "Fields"; }
  else if (cat_name.substr(0, 7) == "UnuAtoM") { cat_name = "Unu"; }
  else if (cat_name.substr(0, 7) == "UnuNtoZ") { cat_name = "Unu"; }

  // try the small version of the shared library
  libname = cat_bname+cat_name+lib_ext;
#ifndef _WIN32
  // don't add the "lib" under windows
  libname = std::string("lib") + libname;
#endif
  LIBRARY_HANDLE category_so = findLib(libname,errorstr);
  
  if(NULL == category_so) {
    log("could not load category lib '" + libname + "': " + errorstr + "\n");
    errstr = std::string(" - ")+errorstr+std::string("\n");
  }

  LIBRARY_HANDLE executable = GetLibraryHandle("");

  if (!category_so && !package_so && !executable) {
    printMessage("Unable to load all of package '" + moduleInfo->package_name_ +
		 "' (category '" + moduleInfo->category_name_ + "' failed) :\n" 
		 + errstr);
    std::cerr << ("Unable to load all of package '" + moduleInfo->package_name_ +
             "' (category '" + moduleInfo->category_name_ + "' failed) :\n"
             + errstr);
    return false;
  }

  std::string makename = "make_" + moduleInfo->module_name_;
  std::vector<LIBRARY_HANDLE> libraries;
  libraries.push_back(category_so);
  libraries.push_back(package_so);
  libraries.push_back(executable);

  std::string error;
  std::ostringstream msg;
  for(std::vector<LIBRARY_HANDLE>::const_iterator lib = libraries.begin();
      lib != libraries.end(); ++lib) {
    if(*lib) {
      moduleInfo->maker_ = reinterpret_cast<ModuleMaker>
                             (GetHandleSymbolAddress(*lib, makename, error));
      if(moduleInfo->maker_) {
        return true;
      }
      msg.clear(); msg.str("");
      msg << "Unable to load symbol '" << makename << "' from module '"
          << moduleInfo->module_name_ << "': " << error << std::endl;
      log(msg.str());
    } else {
      // Unfortunately we don't know *which* library it is, at this point.
      log("Could not load library.\n");
    }
  }
  return false;
}


void
PackageDB::loadPackage(bool resolve)
{
  std::string loadPackage;
  std::string result;
  std::map<int,package*> packages;
  package* new_package = 0;
  category* new_category = 0;
  ModuleInfo* new_module = 0;
  module_iter mi;
  category_iter ci;
  package_iter pi;
  std::string packageElt;
  int mod_count = 0;
  std::string notset("notset");
  std::string packagePath;

  printMessage("Loading packages, please wait...");

  LIBXML_TEST_VERSION;

  // the format of PACKAGE_PATH is a colon seperated list of paths to the
  // root(s) of package source trees.
  const char *srcdir = sci_getenv("SCIRUN_SRCDIR");
  ASSERT(srcdir);

  std::string dtd = std::string(srcdir) +
    std::string("/Dataflow/XML/component.dtd");

  xmlInitializeCatalog();
  xmlCatalogAdd(XMLUtil::char_to_xmlChar("public"),
                XMLUtil::char_to_xmlChar("-//SCIRun/Component DTD"),
                XMLUtil::char_to_xmlChar(dtd.c_str()));

  packagePath = srcdir + std::string("/Packages");

  // if the user specififes it, build the complete package path
  const char *packpath = sci_getenv("SCIRUN_PACKAGE_SRC_PATH");
  if (packpath) packagePath = packagePath + ":" + std::string(packpath);

  // the format of LOAD_PACKAGE is a comma seperated list of package names.
  // build the complete list of packages to load
  ASSERT(sci_getenv("SCIRUN_LOAD_PACKAGE"));
  loadPackage = std::string(sci_getenv("SCIRUN_LOAD_PACKAGE"));

  while(loadPackage!="") {
    // Strip off the first element, leave the rest for the next
    // iteration.
    const unsigned int firstComma = loadPackage.find(',');
    if(firstComma < loadPackage.size()) {
      packageElt=loadPackage.substr(0,firstComma);
      loadPackage=loadPackage.substr(firstComma+1);
    } else {
      packageElt=loadPackage;
      loadPackage="";
    }

    std::string tmpPath = packagePath;
    std::string pathElt;
    std::string finalPath;

    for (;tmpPath!="";) 
    {
      if (packageElt=="SCIRun") 
      {
        tmpPath = "found";
        break;
      }
#ifdef _WIN32
      // don't find the drive letter name's ':'...
      const unsigned int firstColon = tmpPath.find(':',2);
#else
      const unsigned int firstColon = tmpPath.find(':');
#endif
      if(firstColon < tmpPath.size()) 
      {
        pathElt=tmpPath.substr(0,firstColon);
        tmpPath=tmpPath.substr(firstColon+1);
      } 
      else 
      {
        pathElt=tmpPath;
        tmpPath="";
      }
      
      finalPath = pathElt+"/"+packageElt+"/Dataflow";
      if (validDir(finalPath)) {
        tmpPath = "found";
        break;
      }
      // external Package - do Package/src/Dataflow
      finalPath = pathElt+"/"+packageElt+"/src/Dataflow";
      if (validDir(finalPath)) {
        tmpPath = "found";
        break;
      }
      // We may be building an external package with SCIRun compatible modules.
      // Assume that SCIRun_PACKAGE_SRC_PATH points to the right place and
      // simply look for /Dataflow
      finalPath = pathElt+"/Dataflow";
      if (validDir(finalPath)) {
        tmpPath = "found";
        break;
      }
    }
    
    if (tmpPath=="") {
      printMessage("Unable to load package " + packageElt +
		   ":\n - Can't find " + packageElt +
		   " directory in package path");
      continue;
    }

    std::string xmldir;
    
    if(packageElt == "SCIRun") 
    {
      xmldir = std::string(srcdir) + "/Dataflow/XML";
      TCLInterface::execute("lappend auto_path \"" + std::string(srcdir) + "/Dataflow/GUI\"");
    } 
    else 
    {
      xmldir = finalPath+"/XML";
      TCLInterface::execute(std::string("lappend auto_path \"") + finalPath + "/GUI\"");
    }
    std::map<int,char*>* files;
    files = GetFilenamesEndingWith((char*)xmldir.c_str(),".xml");

    // WARNING... looks like the 'files' variable is memory leaked...
    // both in the failure case directly below and in the success case.

    if (!files || files->size() == 0) {
      printMessage("Unable to load package " + packageElt +
		   ":\n - Couldn't find *.xml in " + xmldir );
      continue;
    }

    new_package = new package;
    new_package->name = packageElt;
    packages.insert(std::pair<int, package*>(packages.size(), new_package));

    mod_count += files->size();

    for (char_iter i=files->begin(); i != files->end();  i++) 
    {
      new_module = new ModuleInfo;
      new_module->package_name_ = packageElt;
      new_module->maker_ = 0;
      if (! read_component_file(*new_module, 
				(xmldir+"/"+(*i).second).c_str())) 
      {
        printMessage("Unable to read or validate " + 
		    xmldir+"/"+(*i).second + "\n  Module not loaded.\n");
        continue;
      }

      std::string cat_name = new_module->category_name_;

      ci = new_package->categories.find(cat_name);
      if (ci==new_package->categories.end()) 
      {
        new_category = new category;
        new_category->name = cat_name;
        new_package->categories.insert(std::pair<std::string,
				       category*>(cat_name, 
						  new_category));
        ci = new_package->categories.find(std::string(new_category->name));
      }
      
      mi = (*ci).second->modules.find(new_module->module_name_);
      if (mi==(*ci).second->modules.end()) 
      {
        (*ci).second->modules.insert(std::pair<std::string,
				     ModuleInfo*>(new_module->module_name_,
						  new_module));
      }
    }
  } // end while(loadPackage!="")

  // Update the progress bar to know how many (more) operations will
  // be taking place so it can accurately display the loading
  // progress.
  TCLInterface::execute("addProgressSteps " + to_string(mod_count));

  int index = 0;
  int numreg;
  
  for (pi = packages.begin();
       pi!=packages.end();
       pi++) 
  {

    numreg = 0;
    
    std::string pname = (*pi).second->name;

    printMessage("Loading package '" + pname + "'");
    TCLInterface::execute("setProgressText {Loading package: " + pname + " }");

    for (ci = (*pi).second->categories.begin();
          ci!=(*pi).second->categories.end();
          ci++) 
    {
      for (mi = (*ci).second->modules.begin();
           mi!=(*ci).second->modules.end();
           mi++) 
      {
        if(resolve)
        {
          if(findMaker((*mi).second))
          {
            registerModule((*mi).second);
            numreg++;
          } 
          else  
          {
            std::string mname = (*mi).second->module_name_;
            if (! ((*mi).second)->optional_) {
              printMessage("Unable to load module '" + mname +
               "' :\n - can't find symbol 'make_" + mname + "'");
            }
          }
        } 
        else 
        {
          numreg++;
          registerModule((*mi).second);
        }
        TCLInterface::execute("incrProgress");
      }
    }
    
    if (numreg) 
    {
      TCLInterface::execute("createPackageMenu " + to_string(index++));
    } 
    else 
    {
      printMessage("Unable to load package " + pname + ":\n"
		   " - could not find any valid modules.");
    }
  }

  printMessage("\nFinished loading packages.");
}
  
void
PackageDB::registerModule(ModuleInfo* info) 
{
  Package* package;
  if(!db_->lookup(info->package_name_,package))
  {
    db_->insert(info->package_name_,package=new Package);
    packageList_.push_back( info->package_name_ );
  }
  
  std::string cat_name = info->category_name_;
  if((cat_name.substr(0, 13) == "Conglomerate_"))
  { cat_name.erase(0,13); }

  Category* category;
  if(!package->lookup(cat_name,category))
    package->insert(cat_name,category=new Category);
  
  ModuleInfo* moduleInfo;
  if(!category->lookup(info->module_name_,moduleInfo)) {
    moduleInfo=new ModuleInfo;
    category->insert(info->module_name_,info);
  } else std::cerr << "WARNING: Overriding multiply registered module "
	      << info->package_name_ << "." << info->category_name_ << "."
	      << info->module_name_ << "\n";  
}
 
Module*
PackageDB::instantiateModule(const std::string& packageName,
                             const std::string& categoryName,
                             const std::string& moduleName,
                             const std::string& instanceName)
{
  Package* package;
  if(!db_->lookup(packageName,package)) 
  {
    std::cerr << "ERROR: Instantiating from nonexistant package " << packageName 
         << "\n";
    return 0;
  }
  
  Category* category;
  if(!package->lookup(categoryName,category)) 
  {
    std::cerr << "ERROR: Instantiating from nonexistant category " << packageName
         << "." << categoryName << "\n";
    return 0;
  }
  
  ModuleInfo* moduleInfo;
  if(!category->lookup(moduleName,moduleInfo)) 
  {
    std::cerr << "ERROR: Instantiating nonexistant module " << packageName 
         << "." << categoryName << "." << moduleName << "\n";
    return 0;
  }

  if(!moduleInfo->maker_){
    if(!findMaker(moduleInfo)){
      std::cerr << "ERROR: Cannot find maker for module: " << packageName 
	   << "." << categoryName << "." << moduleName << "\n";
      return 0;
    }
  }

  GuiContext* module_context = TCLInterface::createContext(instanceName);
  Module *module = (moduleInfo->maker_)(module_context);
  if(!module)
    return 0;
  
  module->set_version(moduleInfo->module_version_);
  
  // Some modules may already know their package and category.
  // If this module doesn't, then set its package and category here.
  std::string unknown("unknown");
  if (unknown == module->package_name_)
    module->package_name_ = packageName;
  if (unknown == module->category_name_)
    module->category_name_ = categoryName;
  
  // copy other fields 
  module->lastportdynamic_ = moduleInfo->last_port_dynamic_;
  
  return module;
}
 
bool
PackageDB::haveModule(const std::string& packageName,
                      const std::string& categoryName,
                      const std::string& moduleName) const
{
  Package* package;
  if(!db_->lookup(packageName,package))
    return false;
  
  Category* category;
  if(!package->lookup(categoryName,category))
    return false;
  
  ModuleInfo* moduleInfo;
  if(!category->lookup(moduleName,moduleInfo))
    return false;

  return true;
}


class ModuleEntry {
  public:

    std::string full_module_name_;
    double      low_version_;
    double      high_version_;

    std::string replacement_module_;
    double      replacement_version_;
};


bool
PackageDB::replaceDeprecatedModule(const std::string &pname,
                                   const std::string &cname,
                                   const std::string &mname,
                                   const std::string &vname,
                                   std::string &npname,
                                   std::string &ncname,
                                   std::string &nmname,
                                   std::string &nvname) const
{
  std::string full_module_name = pname + "::" + cname + "::" + mname;
  double version;
  from_string(vname,version);

  static bool loaded = false;
  static std::multimap<std::string, ModuleEntry > replacement_map;

  if (!loaded)
  {
    std::string rmfname = std::string(sci_getenv("SCIRUN_SRCDIR"))
      + "/scripts/module-remapping.txt";
    std::ifstream ifile(rmfname.c_str());
    char buffer[4096];
    while (ifile)
    {
      ifile.getline(buffer, 4096);
      std::string line(buffer);
      if (line.size() < 10) break;
      std::string before = line.substr(0, line.find_first_of(" \t\n"));
      std::string after = line.substr(line.find_last_of(" \t\n", line.size()-1)+1);
      
      ModuleEntry me;
      double low = 0.0; 
      double high = 1.0;
      double repver = 1.0;
      
      std::string::size_type verloc = before.find("#");
      if (verloc != std::string::npos)
      {
        std::string version = before.substr(verloc+1);
        before = before.substr(0,verloc);
        
        std::string::size_type dashloc = version.find("-");
        if (dashloc != std::string::npos)
        {
          std::string lowver = version.substr(0,dashloc);
          std::string highver = version.substr(dashloc+1);
          from_string(lowver,low);
          from_string(highver,high);
        }
        else
        {
          from_string(version,low);
          high = low;
        }
      }      
      

      verloc = after.find("#");
      if (verloc != std::string::npos)
      {
        std::string version = after.substr(verloc+1);
        after = after.substr(0,verloc);
        from_string(version,repver);
      } 
      
      me.full_module_name_ = before;
      me.low_version_ =  low;
      me.high_version_ = high;
      me.replacement_module_ = after;
      me.replacement_version_ = repver;
            
      replacement_map.insert(std::make_pair(before,me));
    }
    ifile.close();
    loaded = true;
  }
  
  //TODO: typedef
  std::pair<std::multimap<std::string, ModuleEntry>::iterator,std::multimap<std::string, ModuleEntry>::iterator> range;

  range = replacement_map.equal_range(full_module_name);
  double oldver;
  from_string(vname,oldver);

  while (range.first != range.second)
  {
    if ((((*(range.first)).second).low_version_ <= oldver) && (((*(range.first)).second).high_version_ >= oldver))
    {
      const std::string &replacement = ((*(range.first)).second).replacement_module_;
      size_t l0 = replacement.find_first_of(':', 0);
      size_t l1 = replacement.find_first_of(':', l0+2);
      npname = replacement.substr(0, l0);
      ncname = replacement.substr(l0+2, l1-l0-2);
      nmname = replacement.substr(l1+2);
      nvname = to_string(((*(range.first)).second).replacement_version_);
      return true;
    }
    range.first++;
  }
  
  npname = pname;
  ncname = cname;
  nmname = mname;
  nvname = vname;
  
  return (false);
}


std::vector<std::string>
PackageDB::packageNames(void) const
{
  // packageList_ is used to keep a list of the packages 
  // that are in this PSE IN THE ORDER THAT THEY ARE SPECIFIED
  // by the user in the Makefile (for main.cc) or in their
  // environment.
  
  return packageList_;
}

std::vector<std::string>
PackageDB::categoryVNames(const std::string& packageName) const
{
  Package* package;
  if(!db_->lookup(packageName, package))
  {
    std::cerr << "WARNING: Unknown package " << packageName << std::endl;
    std::vector<std::string> result(0);
    return result;
  }
  std::vector<std::string> result;
  PackageIter iter(package);
  
  for(iter.first();iter.ok();++iter)
  {
    // Do not list category if it does not contain any modules
    std::vector<std::string> test = moduleVNames(packageName,iter.get_key());
    if (test.size()) result.push_back(iter.get_key());
  }
  return result;
}

std::vector<std::string>
PackageDB::categoryNames(const std::string& packageName) const
{
  Package* package;
  if(!db_->lookup(packageName, package))
  {
    std::cerr << "WARNING: Unknown package " << packageName << std::endl;
    std::vector<std::string> result(0);
    return result;
  }
  std::vector<std::string> result;
  PackageIter iter(package);
  
  for(iter.first();iter.ok();++iter)
  {
    // Do not list category if it does not contain any modules
    result.push_back(iter.get_key());
  }
  return result;
}


std::string
PackageDB::getCategoryName(const std::string &packName,
			   const std::string &catName,
			   const std::string &modName)
{
  Package *package;
  if (!db_->lookup(packName, package)){
    std::cerr << "WARNING: Unknown package " << packName << "\n";
    return catName;
  }

  Category *category;
  ModuleInfo* modinfo;
  if (package->lookup(catName, category) &&
      category->lookup(modName, modinfo))
  {
    // Original category was fine, just return that.
    return catName;
  }

  // Look up the package name somewhere else.  Find a remapping.
  PackageIter iter(package);
  for (iter.first(); iter.ok();++iter)
  {
    if (iter.get_data()->lookup(modName, modinfo))
    {
      return iter.get_key();
    }
  }
  return catName;
}
 

std::vector<std::string>
PackageDB::moduleNames(const std::string& packageName,
		       const std::string& categoryName) const
{
  Package* package;
  if(!db_->lookup(packageName, package)){
    std::cerr << "WARNING: Unknown package " << packageName << "\n";
    std::vector<std::string> result(0);
    return result;
  }

  Category* category;
  if(!package->lookup(categoryName, category)){
    std::cerr << "WARNING: Unknown category " << packageName << "."
	 << categoryName << "\n";
    std::vector<std::string> result(0);
    return result;
  }
  
  std::vector<std::string> result;
  CategoryIter iter(category);
  
  for(iter.first();iter.ok();++iter) 
  {
    result.push_back(iter.get_key());
  }
  return result;
}

std::vector<std::string>
PackageDB::moduleVNames(const std::string& packageName,
		       const std::string& categoryName) const
{
  Package* package;
  if(!db_->lookup(packageName, package)){
    std::cerr << "WARNING: Unknown package " << packageName << "\n";
    std::vector<std::string> result(0);
    return result;
  }

  Category* category;
  if(!package->lookup(categoryName, category)){
    std::cerr << "WARNING: Unknown category " << packageName << "."
	 << categoryName << "\n";
    std::vector<std::string> result(0);
    return result;
  }
  
  std::vector<std::string> result;
  CategoryIter iter(category);
  
  for(iter.first();iter.ok();++iter) 
  {
    if (iter.get_data()->hide_ == false)
    {
      result.push_back(iter.get_key());
    }
  }
  return result;
}

ModuleInfo*
PackageDB::GetModuleInfo(const std::string& name,
			 const std::string& catname,
			 const std::string& packname)
{
  Package* package;
  if (!db_->lookup(packname,package))
    return 0;

  Category* category;
  if (!package->lookup(catname,category))
    return 0;

  ModuleInfo* info;
  if (category->lookup(name,info))
    return info;
  return 0;
}


void
PackageDB::printMessage(const std::string &msg) 
{
  TCLInterface::post_message(msg);
  TCLInterface::execute("update idletasks");
}

  
