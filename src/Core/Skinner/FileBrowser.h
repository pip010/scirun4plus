//
//  For more information, please see: http://software.sci.utah.edu
//
//  The MIT License
//
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//
//  License for the specific language governing rights and limitations under
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
//    File   : FileBrowser.cc
//    Author : McKay Davis
//    Date   : Sun Feb 11 01:47:22 2007

#ifndef SKINNER_FILEBROWSER_H
#define SKINNER_FILEBROWSER_H

#include <Core/Skinner/Parent.h>

namespace SCIRun {
namespace Skinner {

class FileBrowser : public Parent {
public:
  FileBrowser(Variables *variables);
  CatcherFunction_t                 process_event;
  virtual int                       get_signal_id(const string &) const;

protected:
  CatcherFunction_t                 rescan;
  CatcherFunction_t                 do_KeyEvent;
  CatcherFunction_t                 file_selected;

  Var<string>                       filename_;
  Var<string>                       directory_;
  string                            cached_filename_;
};


}
}

#endif
