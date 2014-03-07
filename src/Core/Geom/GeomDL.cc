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
 *  GeomDL.cc: Create a display list for its child
 *
 *  Written by:
 *   Author Yarden Livnat
 *   Department of Computer Science
 *   University of Utah
 *   Date July 2000
 *
 */

#include <Core/Util/Debug.h>

#include <Core/Geom/GeomDL.h>
#include <Core/Geom/DrawInfoOpenGL.h>

#include <list>
#include <algorithm>

namespace SCIRun {

GeomDL::GeomDL(GeomHandle obj)
  : GeomContainer(obj)
{
  DEBUG_CONSTRUCTOR("GeomDL")
}

GeomDL::GeomDL(const GeomDL &copy)
  : GeomContainer(copy)
{
  DEBUG_CONSTRUCTOR("GeomDL")
}

GeomObj*
GeomDL::clone()
{
  return new GeomDL(*this);
}

void
GeomDL::dl_register(DrawInfoOpenGL *info)
{
  drawinfo_.push_back(info);
}

void
GeomDL::dl_unregister(DrawInfoOpenGL *info)
{
  drawinfo_.erase(std::remove(drawinfo_.begin(), drawinfo_.end(), info),
		  drawinfo_.end());
}

GeomDL::~GeomDL()
{
  std::list<DrawInfoOpenGL *>::iterator itr = drawinfo_.begin();
  while (itr != drawinfo_.end())
  {
    (*itr)->dl_remove(this);
    ++itr;
  }
  DEBUG_DESTRUCTOR("GeomDL")  
}

void
GeomDL::reset_bbox()
{
  GeomContainer::reset_bbox();

  // Update all display lists to bad state, forces redraw.
  std::list<DrawInfoOpenGL *>::iterator itr = drawinfo_.begin();
  while (itr != drawinfo_.end())
  {
    (*itr)->dl_update(this, 0xffffffff);
    ++itr;
  }
}

void
GeomDL::fbpick_draw(DrawInfoOpenGL* di, Material *m, double time)
{
  if ( !child_.get_rep() ) return;
  child_->fbpick_draw(di,m,time);  // do not use display list
}

void
GeomDL::draw(DrawInfoOpenGL* di, Material *m, double time)
{
  if ( !child_.get_rep() ) return;

  if ( !pre_draw(di, m, 0) ) return;

  if ( !di->display_list_p_ )
  {
    child_->draw(di,m,time);  // do not use display list
  }
  else
  {
    // Compute current state.
    const unsigned int current_state =
      (di->get_drawtype() << 1) | di->lighting_;

    unsigned int state, display_list;
    if (di->dl_lookup(this, state, display_list))
    {
      if (state != current_state)
      {
        di->dl_update(this, current_state);

        const int pre_polycount_ = di->polycount_; // remember poly count
        
        // Fill in the display list.
        // Don't use COMPILE_AND_EXECUTE as it is slower (NVidia linux).
        glNewList(display_list, GL_COMPILE);
        child_->draw(di,m,time);
        glEndList();
        glCallList(display_list);

        // Update poly count;
        polygons_ = di->polycount_ - pre_polycount_;
      }
      else
      {
        // Display the child using the display_list.
        glCallList(display_list);
        di->polycount_ += polygons_;
      }
    }
    else if (di->dl_addnew(this, current_state, display_list))
    {
      const int pre_polycount_ = di->polycount_; // remember poly count
        
      // Fill in the display list.
      // Don't use COMPILE_AND_EXECUTE as it is slower (NVidia linux).
      glNewList(display_list, GL_COMPILE);
      child_->draw(di,m,time);
      glEndList();
      glCallList(display_list);

      // Update poly count;
      polygons_ = di->polycount_ - pre_polycount_;
    }
    else
    {
      child_->draw(di,m,time);  // Do not use display list.
    }
  }

  post_draw(di);
}

} // End namespace SCIRun

