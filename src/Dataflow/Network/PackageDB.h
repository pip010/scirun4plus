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


// PackageDB.h - Interface to module-finding and loading mechanisms

#ifndef PSE_Dataflow_PackageDB_h
#define PSE_Dataflow_PackageDB_h 1

#include <Core/Containers/AVLTree.h>
#include <Core/Util/soloader.h>
#include <Dataflow/GuiInterface/TCLInterface.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/share.h>

#include <vector>
#include <string>

namespace SCIRun {

typedef struct {
  std::string name;
  std::string datatype;
  iport_maker maker;
} IPortInfo;

typedef struct {
  std::string name;
  std::string datatype;
  oport_maker maker;
} OPortInfo;

struct ModuleInfo 
{
  std::string                       package_name_;
  std::string                       category_name_;
  std::string                       module_name_;
  std::string                       module_version_;
  bool                              optional_;
  bool                              hide_;
  bool                              dynamic_;
  std::vector<std::string>          authors_;
  std::string                       summary_;
  ModuleMaker                       maker_;
  std::vector<IPortInfo*>      iports_;
  std::vector<OPortInfo*>      oports_;
  bool                              last_port_dynamic_;
  bool                              has_gui_node_;
};

typedef AVLTree<std::string,ModuleInfo*> Category;
typedef AVLTree<std::string,Category*> Package;
typedef AVLTree<std::string,Package*> Packages;
    
typedef AVLTreeIter<std::string,ModuleInfo*> CategoryIter;
typedef AVLTreeIter<std::string,Category*> PackageIter;
typedef AVLTreeIter<std::string,Package*> PackagesIter;

class SCISHARE PackageDB {
public:
  PackageDB();
  ~PackageDB();

  void		loadPackage(bool resolve=true);
  Module*		instantiateModule(const std::string& packageName,
					  const std::string& categoryName,
					  const std::string& moduleName,
					  const std::string& instanceName);
  bool		haveModule(const std::string& packageName,
			   const std::string& categoryName,
			   const std::string& moduleName) const;

  bool          replaceDeprecatedModule(const std::string &packageName,
                                        const std::string &categoryName,
                                        const std::string &moduleName,
                                        const std::string &moduleVersion,
                                        std::string &newPackage,
                                        std::string &newCategory,
                                        std::string &newModule,
                                        std::string &newVersion) const;

  std::vector<std::string>	packageNames () const;
  std::vector<std::string>	categoryNames(const std::string& packageName) const;
  std::vector<std::string>	categoryVNames(const std::string& packageName) const;
  std::vector<std::string>	moduleNames  (const std::string& packageName, const std::string& categoryName) const;
  std::vector<std::string>	moduleVNames  (const std::string& packageName, const std::string& categoryName) const;

  ModuleInfo*	GetModuleInfo(const std::string& name,
			      const std::string& catname,
			      const std::string& packname);
  // Used if the module has changed categories.
  std::string		getCategoryName(const std::string &packName,
					const std::string &catName,
					const std::string &modName);
private:

  bool		findMaker(ModuleInfo* info);
  void		registerModule(ModuleInfo* info);
  void		printMessage(const std::string&);

  Packages *        db_;
  std::vector<std::string>    packageList_;
};

// PackageDB is intended to be a singleton class, but nothing will break
// if you instantiate it many times.  This is the singleton instance,
// on which you should invoke operations:
#if defined(_WIN32) && !defined(BUILD_Dataflow_Network)
__declspec(dllimport) PackageDB* packageDB;
#elif defined(_WIN32)
    __declspec(dllexport) extern PackageDB* packageDB;
#else
extern PackageDB* packageDB;
#endif

} // End namespace SCIRun

#endif // PSE_Dataflow_PackageDB_h
