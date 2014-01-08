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

#ifndef DATAFLOW_GUIINTERFACE_TCLINTERFACE_H
#define DATAFLOW_GUIINTERFACE_TCLINTERFACE_H

#include <Core/Thread/Thread.h>
#include <Core/Thread/Mutex.h>
#include <Core/Thread/ConditionVariable.h>
#include <Core/Thread/RecursiveMutex.h>

#include <Dataflow/GuiInterface/GuiCallback.h>

#include <Dataflow/GuiInterface/share.h>

namespace SCIRun {

class GuiContext;
class GuiCallback;

class SCISHARE TCLInterface {

  public:
    static void start_processing_events();

    static void execute(const std::string& str, GuiContext* ctx = 0);
    static void async_execute(const std::string& str, GuiContext* ctx = 0);

    static int eval(const std::string& str, std::string& result, GuiContext* ctx = 0);
    static std::string eval(const std::string &, GuiContext* ctx = 0);
    static void source_once(const std::string&);
    static void add_command(const std::string&, GuiCallback*, void*);
    static void delete_command( const std::string& command );
    
    static void lock();
    static void unlock();
    static bool try_lock();

    static GuiContext* createContext(const std::string& name);
    static void post_message(const std::string& errmsg, bool err = false);
    // Get TCL array var, which resembles a STL a map<string, string>

    static bool get(const std::string& name, std::string& value, GuiContext* ctx = 0);
    static bool get(const std::string& name, double& value, GuiContext* ctx = 0);
    static bool get(const std::string& name, int& value, GuiContext* ctx = 0);

    static void async_get(const std::string& name, std::string* value, GuiContext* ctx = 0);
    static void async_get(const std::string& name, double* value, GuiContext* ctx = 0);
    static void async_get(const std::string& name, int* value, GuiContext* ctx = 0);

    static bool set(const std::string& name, const std::string& value, GuiContext* ctx = 0);
    static bool set(const std::string& name, const double& value, GuiContext* ctx = 0);
    static bool set(const std::string& name, const int& value, GuiContext* ctx = 0);

    static void async_set(const std::string& name, const std::string& value, GuiContext* ctx = 0);
    static void async_set(const std::string& name, const double& value, GuiContext* ctx = 0);
    static void async_set(const std::string& name, const int& value, GuiContext* ctx = 0);

    // Force ansynchronous calls to be completed
    static void synchronize();

    //! Obtain the tcl pause mutex, which will force TCL to lock until
    //! it is released. This functionality is needed for OpenGL synchronization
    static void obtain_tcl_pause();
    static void release_tcl_pause();    

    //! Mark the current thread as the thread having access to TCL
    static void mark_as_tcl_thread();
    
    //! Check whether we are inside the TCL Thread
    static bool is_tcl_thread();
    
    //! Check whether TCL is in use
    static bool tcl_in_use();

    static void print_stats();

    static int OK;

  private:
    // Get an element of regular tcl list
    static bool extract_element_from_list(const std::string& list_contents, 
                                          const std::vector<int> &indexes, 
                                          std::string& value);

    static bool set_element_in_list(std::string& list_contents, 
                                    const std::vector<int> &indexes, 
                                    const std::string& value);

    static bool	string_is_map(const std::string &str);
    static bool	string_is_list_element(const std::string &str);

    static std::string get_map_key_from_string(const std::string &varname); 
    static std::string get_map_name_from_string(const std::string &varname);
    static int pop_last_list_index_from_string(std::string &str);
    static std::string convert_double_for_display(double value, GuiContext* ctx);
};

}

#endif

