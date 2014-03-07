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
 *  GeomPick.h: Picking information for Geometry objects
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#ifndef SCI_Geom_Pick_h
#define SCI_Geom_Pick_h 1

#include <Core/Datatypes/Material.h>

#include <Core/Geom/GeomContainer.h>
#include <Core/Geometry/Vector.h>

#include <vector>

#include <Core/Geom/share.h>

namespace SCIRun {

class BaseWidget;
class MessageBase;
class ViewWindow;
class WidgetPickable;
class BState;

class SCISHARE GeomPick : public GeomContainer
{
  public:
    explicit GeomPick(GeomHandle);
    GeomPick(GeomHandle, WidgetPickable*, int widget_data);
    GeomPick(GeomHandle, const Vector&);
    GeomPick(GeomHandle, const Vector&, const Vector&);
    GeomPick(GeomHandle, const Vector&, const Vector&, const Vector&);
    GeomPick(GeomHandle, const Array1<Vector>&);
    virtual ~GeomPick();
    virtual GeomObj* clone();

    int nprincipal();
    const Vector &principal(int i);
    void set_principal();
    void set_principal(const Vector&);
    void set_principal(const Vector&, const Vector&);
    void set_principal(const Vector&, const Vector&, const Vector&);
    void set_highlight(const MaterialHandle& matl);
    void set_module_data(void*);
    void set_widget_data(int);
    
    void set_picked_obj(GeomHandle);
    void pick(ViewWindow* viewwindow, const BState& bs);
    void moved(int axis, double distance, const Vector& delta, const BState& bs,
         const Vector &pick_offset);
    void release(const BState& bs);
    
    void ignore_until_release();
    
    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  private:
    void*             cbdata_;
    GeomHandle        picked_obj_;
    std::vector<Vector>    directions_;
    WidgetPickable*   widget_;
    int               widget_data_;
    bool              selected_;
    bool              ignore_;
    MaterialHandle    highlight_;
    bool              draw_only_on_pick_;
    
    GeomPick(const GeomPick&);
};

typedef LockingHandle<GeomPick> GeomPickHandle;
  
} // End namespace SCIRun

#endif /* SCI_Geom_Pick_h */
