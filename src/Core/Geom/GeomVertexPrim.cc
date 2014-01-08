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
 *  GeomVertexPrim.cc: Base class for primitives that use the Vertex class
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   February 1995
 *
 */

#include <Core/Util/Debug.h>
#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomVertexPrim.h>

namespace SCIRun {

GeomVertexPrim::GeomVertexPrim()
{
  DEBUG_CONSTRUCTOR("GeomVertexPrim")
}

GeomVertexPrim::GeomVertexPrim(const GeomVertexPrim& copy)
: GeomObj(copy), verts(copy.verts.size())
{
  DEBUG_CONSTRUCTOR("GeomVertexPrim")
  for(int i=0;i<verts.size();i++)
    verts[i]=copy.verts[i]->clone();
}

GeomVertexPrim::~GeomVertexPrim()
{
  for(int i=0;i<verts.size();i++)
    delete verts[i];
  DEBUG_DESTRUCTOR("GeomVertexPrim")
}

void 
GeomVertexPrim::get_bounds(BBox& bb)
{
  for(int i=0;i<verts.size();i++)
    bb.extend(verts[i]->p);
}

void 
GeomVertexPrim::add(const Point& p)
{
  verts.add(new GeomVertex(p));
}

void 
GeomVertexPrim::add(const Point& p, const Vector& normal)
{
  verts.add(new GeomNVertex(p, normal));
}

void 
GeomVertexPrim::add(const Point& p, const MaterialHandle& matl)
{
  verts.add(new GeomMVertex(p, matl));
}

void 
GeomVertexPrim::add(const Point& p, const Color& clr)
{
  verts.add(new GeomCVertex(p, clr));
}

void 
GeomVertexPrim::add(const Point& p, const Vector& normal,
			 const MaterialHandle& matl)
{
  verts.add(new GeomNMVertex(p, normal, matl));
}

void 
GeomVertexPrim::add(GeomVertex* vtx)
{
  verts.add(vtx);
}

void 
GeomVertexPrim::add(double t, GeomVertex* vtx)
{
  times.add(t);
  verts.add(vtx);
}

GeomVertex::GeomVertex(const Point& p)
: p(p)
{
  DEBUG_CONSTRUCTOR("GeomVertex")
}

GeomVertex::GeomVertex(const GeomVertex& copy)
: Persistent(copy),
  p(copy.p)
{
  DEBUG_CONSTRUCTOR("GeomVertex")
}

GeomVertex::~GeomVertex()
{
  DEBUG_DESTRUCTOR("GeomVertex")
}

GeomVertex* 
GeomVertex::clone()
{
  return new GeomVertex(*this);
}

void
GeomVertex::emit_all(DrawInfoOpenGL*)
{
  glVertex3d(p.x(), p.y(), p.z());
}


void
GeomVertex::emit_point(DrawInfoOpenGL*)
{
  glVertex3d(p.x(), p.y(), p.z());
}


void
GeomVertex::emit_material(DrawInfoOpenGL*)
{
  // Do nothing
}


void
GeomVertex::emit_normal(DrawInfoOpenGL*)
{
  // Do nothing
}

GeomNVertex::GeomNVertex(const Point& p, const Vector& normal)
: GeomVertex(p), normal(normal)
{
  DEBUG_CONSTRUCTOR("GeomNVertex")
}

GeomNVertex::GeomNVertex(const GeomNVertex& copy)
: GeomVertex(copy), normal(copy.normal)
{
  DEBUG_CONSTRUCTOR("GeomNVertex")
}

GeomVertex* 
GeomNVertex::clone()
{
  return new GeomNVertex(*this);
}

GeomNVertex::~GeomNVertex()
{
  DEBUG_DESTRUCTOR("GeomNVertex")
}

void
GeomNVertex::emit_all(DrawInfoOpenGL*)
{
  glNormal3d(normal.x(), normal.y(), normal.z());
  glVertex3d(p.x(), p.y(), p.z());
}

void
GeomNVertex::emit_normal(DrawInfoOpenGL*)
{
  glNormal3d(normal.x(), normal.z(), normal.z());
}

GeomNMVertex::GeomNMVertex(const Point& p, const Vector& normal,
			   const MaterialHandle& matl)
: GeomNVertex(p, normal), matl(matl)
{
  DEBUG_CONSTRUCTOR("GeomNMVertex")
}

GeomNMVertex::GeomNMVertex(const GeomNMVertex& copy)
: GeomNVertex(copy), matl(copy.matl)
{
  DEBUG_CONSTRUCTOR("GeomNMVertex")
}

GeomVertex* 
GeomNMVertex::clone()
{
  return new GeomNMVertex(*this);
}

GeomNMVertex::~GeomNMVertex()
{
  DEBUG_DESTRUCTOR("GeomNMVertex")
}

void
GeomNMVertex::emit_all(DrawInfoOpenGL* di)
{
  di->set_material(matl.get_rep());
  glNormal3d(normal.x(), normal.y(), normal.z());
  glVertex3d(p.x(), p.y(), p.z());
}

void
GeomNMVertex::emit_material(DrawInfoOpenGL* di)
{
  di->set_material(matl.get_rep());
}

GeomMVertex::GeomMVertex(const Point& p, const MaterialHandle& matl)
: GeomVertex(p), matl(matl)
{
  DEBUG_CONSTRUCTOR("GeomMVertex")
}

GeomMVertex::GeomMVertex(const GeomMVertex& copy)
: GeomVertex(copy), matl(matl)
{
  DEBUG_CONSTRUCTOR("GeomMVertex")
}

GeomVertex* 
GeomMVertex::clone()
{
  return new GeomMVertex(*this);
}

GeomMVertex::~GeomMVertex()
{
  DEBUG_DESTRUCTOR("GeomMVertex")
}

void
GeomMVertex::emit_all(DrawInfoOpenGL* di)
{
  di->set_material(matl.get_rep());
  glVertex3d(p.x(), p.y(), p.z());
}

void
GeomMVertex::emit_material(DrawInfoOpenGL* di)
{
  di->set_material(matl.get_rep());
}

GeomCVertex::GeomCVertex(const Point& p, const Color& clr)
: GeomVertex(p), color(clr)
{
  DEBUG_CONSTRUCTOR("GeomCVertex")
}

GeomCVertex::GeomCVertex(const GeomCVertex& copy)
: GeomVertex(copy), color(copy.color)
{
  DEBUG_CONSTRUCTOR("GeomCVertex")
}

GeomVertex* GeomCVertex::clone()
{
  return new GeomCVertex(*this);
}

GeomCVertex::~GeomCVertex()
{
  DEBUG_DESTRUCTOR("GeomCVertex")
}

void
GeomCVertex::emit_all(DrawInfoOpenGL* /*di*/)
{
  glColor3f(color.r(),color.g(),color.b());
  glVertex3d(p.x(), p.y(), p.z());
}

void
GeomCVertex::emit_material(DrawInfoOpenGL* /*di*/)
{
  glColor3f(color.r(),color.g(),color.b());
}

} // End namespace SCIRun



