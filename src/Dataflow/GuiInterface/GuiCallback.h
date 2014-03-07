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
 *  GuiCallback.h: Interface to user interface
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 2002
 *
 */

#ifndef SCIRun_Core_GuiInterface_GuiCallback_h
#define SCIRun_Core_GuiInterface_GuiCallback_h

#include <vector>
#include <string>
#include <Dataflow/GuiInterface/share.h>

namespace SCIRun {

  class GuiVar;

  class SCISHARE GuiArgs {
    std::vector<std::string> args_;
  public:
    bool have_error_;
    bool have_result_;
    std::string string_;
    
    GuiArgs(int argc, const char* argv[]);
    GuiArgs(int argc, char* argv[]);
    ~GuiArgs();
    int count();
    std::string operator[](int i);
    std::string get_string(int i);
    int get_int(int i);
    double get_double(int i);
    
    void error(const std::string&);
    void result(const std::string&);
    void append_result(const std::string&);
    void append_element(const std::string&);

    static std::string make_list(const std::string&, const std::string&);
    static std::string make_list(const std::string&, const std::string&, const std::string&);
    static std::string make_list(const std::vector<std::string>&);
  };

  class SCISHARE GuiCallback {
  public:
    GuiCallback();
    virtual ~GuiCallback();
    virtual void tcl_command(GuiArgs&, void*)=0;
  };
} // End namespace SCIRun


#endif
