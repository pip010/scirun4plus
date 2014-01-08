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
 *  GuiView.h   Structure that provides for easy access of view information.
 *              The view information is interactively provided by the user.
 *
 *  Written by:
 *   Steven Parker
 *   Department of Computer Science
 *   University of Utah
 *
 *   separated from the Viewer code by me (Aleksandra)
 *   in May 1996
 *
 */

#ifndef SCI_project_GuiView_h
#define SCI_project_GuiView_h 1

#include <Core/Geom/View.h>
#include <Dataflow/GuiInterface/GuiGeom.h>
#include <Dataflow/GuiInterface/GuiVar.h>

#include <Dataflow/GuiInterface/share.h>

namespace SCIRun {


class SCISHARE GuiView : public GuiVar {
  public:
    GuiView(GuiContext* ctx);
    ~GuiView();
    GuiView(const GuiView&);

    View get();
    View get_cached();
    void request();
    
    void set(const View&);
    void send(const View&);

  private:
    GuiPoint  eyep_;
    GuiPoint  lookat_;
    GuiVector up_;
    GuiDouble fov_;
    GuiVector eyep_offset_;
};

class GuiExtendedView : public GuiVar {
  private:
    GuiPoint  eyep_;
    GuiPoint  lookat_;
    GuiVector up_;
    GuiDouble fov_;
    GuiVector eyep_offset_;
    GuiInt    xres_;
    GuiInt    yres_;
    GuiColor  bg_;

  public:
    GuiExtendedView(GuiContext* ctx);
    ~GuiExtendedView();
    GuiExtendedView(const GuiExtendedView&);

    ExtendedView get();
    ExtendedView get_cached();
    void request();
    
    void set(const ExtendedView&);
    void send(const ExtendedView&);
  };

} // End namespace SCIRun


#endif
