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
 *  Geom.cc: Displayable Geometry
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#include <Core/Util/Debug.h>
#include <Core/Geom/GeomObj.h>
#include <Core/Geom/DrawInfoOpenGL.h>

#include <Core/Geometry/Vector.h>
#include <Core/Thread/MutexPool.h>

#include <sci_defs/bits_defs.h>

namespace SCIRun {

MutexPool lock_pool("GeomObj pool", 250);

GeomObj::GeomObj() :
  UsedWithLockingHandle<Mutex&>(lock_pool.getMutex()),
  id_int_(0x1234567),
  id_longlong_((long long)(0x1234567))
{
  DEBUG_CONSTRUCTOR("GeomObj")
}


GeomObj::GeomObj(const GeomObj& obj) :
  UsedWithLockingHandle<Mutex&>(lock_pool.getMutex()),
  id_int_(obj.id_int_),
  id_longlong_(obj.id_longlong_)
{
  DEBUG_CONSTRUCTOR("GeomObj")
}

GeomObj::~GeomObj()
{
  DEBUG_DESTRUCTOR("GeomObj")
}


bool
GeomObj::getId( int& id )
{
  if ( id_int_ == 0x1234567)
    return false;
  else {
    id = id_int_;
    return true;
  }
}

bool
GeomObj::getId( long long& id )
{
  if ( id_longlong_ == (long long)(0x1234567))
    return false;
  else {
    id = id_longlong_;
    return true;
  }
}

void GeomObj::get_triangles( Array1<float> &)
{
  std::cerr << "GeomObj::get_triangles - no triangles" << std::endl;
}

void 
GeomObj::reset_bbox()
{
    // Nothing to do, by default.
}

int
GeomObj::pre_draw(DrawInfoOpenGL* di, Material* matl, int lit)
{
  /* All primitives that get drawn must check the return value of this
     function to determine if they get drawn or not */
  if ((!di->pickmode_)||(di->pickmode_&&di->pickchild_))
  {
    if (lit && di->lighting_ && !di->currently_lit_)
    {
      di->currently_lit_=1;
      glEnable(GL_LIGHTING);
      switch(di->get_drawtype())
      {
      case DrawInfoOpenGL::WireFrame:
        gluQuadricNormals(di->qobj_, (GLenum)GLU_SMOOTH);
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        break;
      case DrawInfoOpenGL::Flat:
        gluQuadricNormals(di->qobj_, (GLenum)GLU_FLAT);
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        break;
      case DrawInfoOpenGL::Gouraud:
        gluQuadricNormals(di->qobj_, (GLenum)GLU_SMOOTH);
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        break;
      }

    }
    if ((!lit || !di->lighting_) && di->currently_lit_)
    {
      di->currently_lit_=0;
      glDisable(GL_LIGHTING);
      gluQuadricNormals(di->qobj_, (GLenum)GLU_NONE);
      switch(di->get_drawtype())
      {
      case DrawInfoOpenGL::WireFrame:
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        break;
      case DrawInfoOpenGL::Flat:
      case DrawInfoOpenGL::Gouraud:
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        break;
      }
    }
    di->set_material(matl);
#ifdef SCI_64BITS
    /// @todo SKETCHY CODE!  Fixme.
    unsigned long long o= reinterpret_cast<unsigned long long>(this);
    unsigned int o1=(o>>32)&0xffffffff;
    unsigned int o2=o&0xffffffff;
    glPushName(o1);
    glPushName(o2);
#else
    glPushName(reinterpret_cast<GLuint>(this));
#endif
    glPushName(0x12345678);
    return 1;
  }
  else return 0;
}


int
GeomObj::post_draw(DrawInfoOpenGL* di)
{
  if (di->pickmode_ && di->pickchild_)
  {
#ifdef SCI_64BITS
    glPopName();
    glPopName();
#else
    glPopName();//pops the face index once the obj is rendered
#endif
    glPopName();
  }
  return 1;  // needed to quiet visual c++
}

void 
GeomObj::get_view(DrawInfoOpenGL* di) 
{
  GLdouble matrix[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, matrix);

  const double lvx = fabs(matrix[2]);
  const double lvy = fabs(matrix[6]);
  const double lvz = fabs(matrix[10]);

  if (lvx >= lvy && lvx >= lvz)
  {
    di->axis_ = 0;
    if (matrix[2] > 0) { di->dir_ = 1; }
    else { di->dir_ = -1; }

  }
  else if (lvy >= lvx && lvy >= lvz)
  {
    di->axis_ = 1;
    if (matrix[6] > 0) { di->dir_ = 1; }
    else { di->dir_ = -1; }
  }
  else if (lvz >= lvx && lvz >= lvy)
  {
    di->axis_ = 2;
    if (matrix[10] > 0) { di->dir_ = 1; }
    else { di->dir_ = -1; }
  }
}


} // End namespace SCIRun
