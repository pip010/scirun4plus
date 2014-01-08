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
 *  GeomObj.h: Displayable Geometry
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#ifndef SCI_Geom_GeomObj_h
#define SCI_Geom_GeomObj_h 1

#include <Core/Containers/Array1.h>
#include <Core/Containers/LockingHandle.h>
#include <Core/Thread/UsedWithLockingHandle.h>
#include <Core/Geometry/BBox.h>

#include <iosfwd>

#include <Core/Geom/share.h>

namespace SCIRun {

struct DrawInfoOpenGL;
class  Material;
class  Hit;

class SCISHARE GeomObj : public UsedWithLockingHandle<Mutex&>
{
  public:
    GeomObj();
    GeomObj(const GeomObj&);
    virtual ~GeomObj();
    virtual GeomObj* clone() = 0;

    virtual void reset_bbox();
    virtual void get_bounds(BBox&) = 0;

    // For OpenGL
    int pre_draw(DrawInfoOpenGL*, Material*, int lit);
    virtual void draw(DrawInfoOpenGL*, Material*, double time)=0;
    virtual void fbpick_draw(DrawInfoOpenGL*, Material*, double /*time*/) {}
    int post_draw(DrawInfoOpenGL*);

    virtual void get_view(DrawInfoOpenGL*);

    virtual void get_triangles( Array1<float> &v );

    // we want to return false if value is the default value
    virtual void setId(int id) { id_int_ = id; }
    virtual void setId(long long id) { id_longlong_ = id; }

    virtual bool getId( int &id );
    virtual bool getId( long long &id );

  private:

    int       id_int_;
    long long  id_longlong_;
};

typedef LockingHandle<GeomObj> GeomHandle;

} // End namespace SCIRun

#endif // ifndef SCI_Geom_GeomObj_h




