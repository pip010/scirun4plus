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
 *  Pick.cc: Picking information for Geometry objects
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#include <Core/Util/Debug.h>
#include <Core/Geom/GeomPick.h>
#include <Core/Geom/DrawInfoOpenGL.h>

#include <Core/GeomInterface/Pickable.h>

#include <sci_defs/bits_defs.h>

namespace SCIRun {

GeomPick::GeomPick(GeomHandle obj)
  : GeomContainer(obj),
    cbdata_(0),
    picked_obj_(0),
    directions_(6),
    widget_(0),
    selected_(false),
    ignore_(false),
    draw_only_on_pick_(false)
{
  directions_[0]=Vector(1,0,0);
  directions_[1]=Vector(-1,0,0);
  directions_[2]=Vector(0,1,0);
  directions_[3]=Vector(0,-1,0);
  directions_[4]=Vector(0,0,1);
  directions_[5]=Vector(0,0,-1);
  DEBUG_CONSTRUCTOR("GeomPick")
}

GeomPick::GeomPick(GeomHandle obj ,
		   WidgetPickable* widget, int widget_data)
  : GeomContainer(obj),
    cbdata_(0),
    picked_obj_(0),
    directions_(6),
    widget_(widget),
    widget_data_(widget_data),
    selected_(false),
    ignore_(false),
    draw_only_on_pick_(false)
{
  directions_[0]=Vector(1,0,0);
  directions_[1]=Vector(-1,0,0);
  directions_[2]=Vector(0,1,0);
  directions_[3]=Vector(0,-1,0);
  directions_[4]=Vector(0,0,1);
  directions_[5]=Vector(0,0,-1);
  DEBUG_CONSTRUCTOR("GeomPick")
}

GeomPick::GeomPick(GeomHandle obj, const Vector& v1)
  : GeomContainer(obj),
    cbdata_(0),
    picked_obj_(0),
    directions_(2),
    widget_(0),
    selected_(false),
    ignore_(false),
    draw_only_on_pick_(false)
{
  directions_[0]=v1;
  directions_[1]=-v1;
  DEBUG_CONSTRUCTOR("GeomPick")
}

GeomPick::GeomPick(GeomHandle obj,
		   const Vector& v1, const Vector& v2)
  : GeomContainer(obj),
    cbdata_(0),
    picked_obj_(0),
    directions_(4),
    widget_(0),
    selected_(false),
    ignore_(false),
    draw_only_on_pick_(false)
{
  directions_[0]=v1;
  directions_[1]=-v1;
  directions_[2]=v2;
  directions_[3]=-v2;
  DEBUG_CONSTRUCTOR("GeomPick")
}

GeomPick::GeomPick(GeomHandle obj,
		   const Vector& v1, const Vector& v2, const Vector& v3)
  : GeomContainer(obj),
    cbdata_(0),
    picked_obj_(0),
    directions_(6),
    widget_(0),
    selected_(false),
    ignore_(false),
    draw_only_on_pick_(false)
{
  directions_[0]=v1;
  directions_[1]=-v1;
  directions_[2]=v2;
  directions_[3]=-v2;
  directions_[4]=v3;
  directions_[5]=-v3;
  DEBUG_CONSTRUCTOR("GeomPick")
}

GeomPick::GeomPick(const GeomPick& copy)
  : GeomContainer(copy),
    cbdata_(copy.cbdata_), 
    picked_obj_(copy.picked_obj_),
    directions_(copy.directions_),
    widget_(copy.widget_),
    selected_(copy.selected_),
    ignore_(copy.ignore_),
    highlight_(copy.highlight_),
    draw_only_on_pick_(copy.draw_only_on_pick_)
{
  DEBUG_CONSTRUCTOR("GeomPick")
}

GeomPick::~GeomPick()
{
  DEBUG_DESTRUCTOR("GeomPick")
}

GeomObj* GeomPick::clone()
{
  return new GeomPick(*this);
}

void
GeomPick::set_highlight(const MaterialHandle& matl)
{
  highlight_ = matl;
}

void
GeomPick::set_module_data(void* cbdata)
{
  cbdata_ = cbdata;
}

void
GeomPick::set_widget_data(int wd)
{
  widget_data_ = wd;
}

void
GeomPick::set_picked_obj(GeomHandle object)
{
  picked_obj_ = object;
}

void
GeomPick::ignore_until_release()
{
  ignore_ = true;
}

void
GeomPick::pick(ViewWindow* viewwindow, const BState& bs )
{
  selected_=true;
  ignore_=false;
  if(widget_)
  {
    widget_->geom_pick(this, viewwindow, widget_data_, bs);
  }
}

void
GeomPick::release(const BState& bs)
{
  selected_=false;
  if(widget_)
  {
    widget_->geom_release(this, widget_data_, bs);
  }
}

void
GeomPick::moved(int axis, double distance, const Vector& delta,
		const BState& bs, const Vector &pick_offset)
{
  if(ignore_) { return; }
  if(widget_)
  {
    widget_->geom_moved(this, axis, distance, delta,
			widget_data_, bs, pick_offset);
  }
}

int
GeomPick::nprincipal()
{
  return directions_.size();
}

const Vector&
GeomPick::principal(int i)
{
  return directions_[i];
}

void
GeomPick::set_principal()
{
  directions_.clear();
}

void
GeomPick::set_principal(const Vector& v1)
{
  directions_.clear();
  directions_.push_back(v1);
  directions_.push_back(-v1);
}

void
GeomPick::set_principal(const Vector& v1, const Vector& v2)
{
  directions_.clear();
  directions_.push_back(v1);
  directions_.push_back(-v1);
  directions_.push_back(v2);
  directions_.push_back(-v2);
}

void
GeomPick::set_principal(const Vector& v1, const Vector& v2, const Vector& v3)
{
  directions_.clear();
  directions_.push_back(v1);
  directions_.push_back(-v1);
  directions_.push_back(v2);
  directions_.push_back(-v2);
  directions_.push_back(v3);
  directions_.push_back(-v3);
}

void
GeomPick::draw(DrawInfoOpenGL* di, Material* matl, double time)
{
  if (draw_only_on_pick_ && !di->pickmode_) return;
  if (di->pickmode_)
  {
    ++di->npicks_;
#ifdef SCI_64BITS
    /// @todo SKETCHY CODE!  Fixme.
    unsigned long long o=reinterpret_cast<unsigned long long>(this);
    unsigned int o1=(o>>32)&0xffffffff;
    unsigned int o2=o&0xffffffff;
    glPushName(o1);
    glPushName(o2);
#else
    glPushName(reinterpret_cast<GLuint>(this));
#endif
    di->pickchild_ =1;
  }
  if (selected_ && highlight_.get_rep())
  {
    di->set_material(highlight_.get_rep());
    int old_ignore=di->ignore_matl_;
    di->ignore_matl_=1;
    child_->draw(di, highlight_.get_rep(), time);
    di->ignore_matl_=old_ignore;
  }
  else
  {
    child_->draw(di, matl, time);
  }
  if (di->pickmode_)
  {
#ifdef SCI_64BITS
    glPopName();
    glPopName();
#else
    glPopName();
#endif
        
    if ((--di->npicks_)<1)
    { // Could have multiple picks in stack.
      di->pickchild_=0;
    }
  }
}

} // End namespace SCIRun

