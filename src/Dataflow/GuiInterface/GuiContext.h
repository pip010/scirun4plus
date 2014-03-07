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
 *  GuiContext.h:
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 2002
 *
 */

#ifndef DATAFLOW_GUIINTERFACE_GUICONTEXT_H
#define DATAFLOW_GUIINTERFACE_GUICONTEXT_H

#include <string>
#include <vector>
#include <iosfwd>

#include <Dataflow/GuiInterface/share.h>

namespace SCIRun {

class SCISHARE GuiContext {
  public:
    GuiContext(const std::string& name, 
               bool save=true,
               GuiContext *parent = 0);
    ~GuiContext();
    
    GuiContext* subVar(const std::string& name, bool save=true);
    GuiContext* get_parent() { return parent_; }

    GuiContext* find_child(const std::string &name);
    
    void			async_get(std::string* value);
    bool			get(std::string& value);
    void			async_set(const std::string& value);
    void			set(const std::string& value);
    
    void			async_get(double* value);
    bool			get(double& value);
    void			async_set(double value);
    void			set(double value);
    
    void			async_get(int* value);
    bool			get(int& value);
    void			async_set(int value);
    void			set(int value);
    
    void      synchronize();
    
    std::string getfullname();

    void			dontSubstituteDatadir();
    void			doSubstituteDatadir();

    void			setIsFilename();
    void			unsetIsFilename();

    void			dontSave();
    void			doSave();
    
    void			dontCache(); // always query GUI for value
    void			doCache();  // only query GUI if not cached already
    void			reset(); // resets the cache

    void      doReplaceEnv();
    void      dontReplaceEnv();
    bool      replace_env();

    bool      is_active() { return (active_); }
    void      inactivate();
    
  private:  
    void			tcl_setVarStates();  
    
    enum  {
      SAVE_E			= 1 << 0,
      CACHE_E			= 1 << 1,
      CACHED_E	    = 1 << 2,
      SUBSTITUTE_DATADIR_E	= 1 << 3,
      IS_FILENAME_E     = 1 << 4,
      SUBSTITUTE_ENVIRONMENT_E = 1<<5
    };
    
    GuiContext *              parent_;
    std::string               name_;
    std::vector<GuiContext*>	children_;
    unsigned int              context_state_;
    
    bool active_;
    
    GuiContext(const GuiContext&);
    GuiContext& operator=(const GuiContext&);
    };

}

#endif

