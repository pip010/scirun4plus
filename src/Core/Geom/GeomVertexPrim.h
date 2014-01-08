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
 *  GeomVertexPrim.h: Base class for primitives that use the Vertex class
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   February 1995
 *
 */

#ifndef CORE_GEOM_GEOMVEXTEXPRIM_H
#define CORE_GEOM_GEOMVEXTEXPRIM_H 1

#include <Core/Geometry/Point.h>
#include <Core/Geom/GeomObj.h>

#include <Core/Datatypes/Material.h>
#include <Core/Datatypes/Color.h>
#include <Core/Containers/Array1.h>

#include <stdlib.h>

#include <Core/Geometry/share.h>

namespace SCIRun {

class GeomVertex : public Persistent
{
  public:
    Point p;
    GeomVertex(const Point& p);
    GeomVertex(const GeomVertex&);
    virtual ~GeomVertex();
    virtual GeomVertex* clone();

    virtual void emit_all(DrawInfoOpenGL* di);
    void emit_point(DrawInfoOpenGL* di);
    virtual void emit_material(DrawInfoOpenGL* di);
    virtual void emit_normal(DrawInfoOpenGL* di);

};

class GeomNVertex : public GeomVertex
{
  public:
    Vector normal;
    GeomNVertex(const Point& p, const Vector& normal);
    GeomNVertex(const GeomNVertex&);
    virtual GeomVertex* clone();
    virtual ~GeomNVertex();

    virtual void emit_all(DrawInfoOpenGL* di);
    virtual void emit_normal(DrawInfoOpenGL* di);
};

class GeomNMVertex : public GeomNVertex
{
  public:
    MaterialHandle matl;
    GeomNMVertex(const Point& p, const Vector& normal,
           const MaterialHandle& matl);
    GeomNMVertex(const GeomNMVertex&);
    virtual GeomVertex* clone();
    virtual ~GeomNMVertex();

    virtual void emit_all(DrawInfoOpenGL* di);
    virtual void emit_material(DrawInfoOpenGL* di);
};

class GeomMVertex : public GeomVertex 
{
  public:
    MaterialHandle matl;
    GeomMVertex(const Point& p, const MaterialHandle& matl);
    GeomMVertex(const GeomMVertex&);
    virtual GeomVertex* clone();
    virtual ~GeomMVertex();

    virtual void emit_all(DrawInfoOpenGL* di);
    virtual void emit_material(DrawInfoOpenGL* di);
};

class GeomCVertex : public GeomVertex 
{
  public:
    Color color;
    GeomCVertex(const Point& p, const Color& clr);
    GeomCVertex(const GeomCVertex&);
    virtual GeomVertex* clone();
    virtual ~GeomCVertex();

    virtual void emit_all(DrawInfoOpenGL* di);
    virtual void emit_material(DrawInfoOpenGL* di);
};

class GeomVertexPrim : public GeomObj
{
  public:
    Array1<double> times;
    Array1<GeomVertex*> verts;

    GeomVertexPrim();
    GeomVertexPrim(const GeomVertexPrim&);
    virtual ~GeomVertexPrim();

    virtual void get_bounds(BBox&);
      
    void add(const Point&);
    void add(const Point&, const Vector&);
    void add(const Point&, const MaterialHandle&);
    void add(const Point&, const Color&);
    void add(const Point&, const Vector&, const MaterialHandle&);
    void add(GeomVertex*);
    void add(double time, GeomVertex*);
};

} // End namespace SCIRun

#endif /* SCI_Geom_VertexPrim_h */
