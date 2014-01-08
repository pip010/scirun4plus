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


// Resources.h - 

#ifndef SCIRun_Resources_h
#define SCIRun_Resources_h 

#include <string>
#include <vector>
#include <map>

namespace SCIRun {

/*
 * Port
 */

class PortInfo {
public:
  std::string type_;
  std::string package_;
  std::string datatype_;
  std::string imaker_;
  std::string omaker_;
  std::vector<std::string> libs_;
};

/*
 * Package
 */

class PackageInfo {
public:
  std::string name_;
  std::string path_;
  std::string lib_path_;
  std::string ui_path_;
  int level_;
};


/*
 * Module
 */

class ModulePortInfo {
public:
  std::string name_;
  std::string type_;
};

class ModuleInfo {
public:
  std::string package_;
  std::string name_;
  std::string id_;
  std::vector<std::string> categories_;
  std::string maker_;
  std::string ui_;
  vector<ModulePortInfo*> iports_;
  vector<ModulePortInfo*> oports_;
  bool has_dynamic_port_;

  std::vector<std::string> libs_;
};


/* 
 * Resources
 */ 

class Resources {
public:
  Resources(void);
  ~Resources(void);
  
  // General Info
  // for compatibility with current NetworkEditor
  std::vector<std::string> get_packages_names();
  std::vector<std::string> get_categories_names( const std::string & );
  std::vector<std::string> get_modules_names( const std::string &, const std::string & );

  // Packages
  PackageInfo *get_package_info( const std::string & );
  std::string get_package_ui_path( const std::string & );

  // Modules
  ModuleInfo *get_module_info( const std::string & );
  const std::vector<std::string> &get_module_libs( const std::string & );
  std::string get_module_maker( const std::string & );
  std::string get_module_ui( const std::string &);

  // Ports
  PortInfo *get_port_info( const std::string & );
  const std::vector<std::string> &get_port_libs( const std::string & );
  std::string get_iport_maker( const std::string & );
  std::string get_oport_maker( const std::string & );

 
  void read( const std::string & );


private:
  typedef std::map<std::string, PackageInfo *> PackagesList;
  typedef std::map<std::string,ModuleInfo *> ModulesList;
  typedef std::map<std::string,PortInfo *> PortsList;

  PackagesList packages_;
  ModulesList modules_;
  PortsList ports_;

  std::string data_path_;

  friend class ResourcesParser;
  friend class PackageParser;
  friend class ModuleParser;
  friend class PortParser;
};

// Resources is intended to be a singleton class, but nothing will break
// if you instantiate it many times.  This is the singleton instance,
// on which you should invoke operations:

extern Resources resources;

} // End namespace SCIRun

#endif // SCIRun_Resources_h
