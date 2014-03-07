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
 *  GeomGrid.cc: Grid object
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   June 1995
 *
 */

#include <Core/Util/Debug.h>

#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomGrid.h>

namespace SCIRun {

GeomGrid::GeomGrid(int nu, int nv, const Point& corner,
		   const Vector& u, const Vector& v, int image)
: verts(nu, nv), corner(corner), u(u), v(v), image(image)
{
  DEBUG_CONSTRUCTOR("GeomGrid")
  have_matls=0;
  have_normals=0;
  adjust();
}

GeomGrid::GeomGrid(const GeomGrid& copy)
: GeomObj(copy)
{
  DEBUG_CONSTRUCTOR("GeomGrid")
}

GeomGrid::~GeomGrid()
{
  DEBUG_DESTRUCTOR("GeomGrid")
}

void 
GeomGrid::adjust()
{
  w=Cross(u, v);
  w.normalize();
}

void 
GeomGrid::set(int i, int j, double v)
{
  verts(i, j)=v;
}

double 
GeomGrid::get(int i, int j)
{
  double v = verts(i, j);
  return v;
}

void 
GeomGrid::set(int i, int j, double v, const Vector& normal)
{
  if(!have_normals)
  {
    normals.resize(verts.dim1(), verts.dim2());
    have_normals=1;
  }
  verts(i, j)=v;
  normals(i, j)=normal;
}

void 
GeomGrid::set(int i, int j, double v, const MaterialHandle& matl)
{
  if(!have_matls)
  {
    matls.resize(verts.dim1(), verts.dim2());
    have_matls=1;
  }
  verts(i, j)=v;
  matls(i, j)=matl;
}

void 
GeomGrid::set(int i, int j, double v, const Vector& normal,
		   const MaterialHandle& matl)
{
  if(!have_matls)
  {
    matls.resize(verts.dim1(), verts.dim2());
    have_matls=1;
  }
 
  if(!have_normals)
  {
    normals.resize(verts.dim1(), verts.dim2());
    have_normals=1;
  }
  verts(i, j)=v;
  matls(i, j)=matl;
  normals(i, j)=normal;
}

void 
GeomGrid::get_bounds(BBox& bb)
{
  int nu=verts.dim1();
  int nv=verts.dim2();
  if (!image) 
  {
    Vector uu(u/(nu-1));
    Vector vv(v/(nv-1));
    Point rstart(corner);
    for(int i=0;i<nu;i++)
    {
	    Point p(rstart);
	    for(int j=0;j<nv;j++)
      {
        Point pp(p+w*verts(i, j));
        bb.extend(pp);
        p+=vv;
	    }
	    rstart+=uu;
    }
  } 
  else 
  {
    bb.extend(corner-Vector(0.001,0.001,0.001));
    bb.extend(corner+u+v+Vector(0.001,0.001,0.001));
  }
}

GeomObj*
GeomGrid::clone()
{
  return new GeomGrid(*this);
}

void
GeomGrid::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  int nu=verts.dim1();
  int nv=verts.dim2();

  if (image)
  {
    di->polycount_+=2*nu*nv;
    Vector uu(u/nu);
    Vector vv(v/nv);
    if (!pre_draw(di,matl,0)) return;
    Point rstart(corner);
    glBegin(GL_QUADS);
    for (int i=0; i<nu; i++)
    {
      Point p1(rstart);
      Point p2(rstart+uu);
      Point p3(rstart+uu+vv);
      Point p4(rstart+vv);
      for (int j=0; j<nv; j++)
      {
        if (have_matls)
        {
          di->set_material(matls(i,j).get_rep());
        }
        glVertex3d(p1.x(), p1.y(), p1.z());
        glVertex3d(p2.x(), p2.y(), p2.z());
        glVertex3d(p3.x(), p3.y(), p3.z());
        glVertex3d(p4.x(), p4.y(), p4.z());
        p1+=vv;
        p2+=vv;
        p3+=vv;
        p4+=vv;
      }
      rstart+=uu;
    }
    glEnd();
    return;
  }

  if (!pre_draw(di, matl, 1)) return;
  di->polycount_+=2*(nu-1)*(nv-1);
  Vector uu(u/(nu-1));
  Vector vv(v/(nv-1));
  switch(di->get_drawtype())
  {
  case DrawInfoOpenGL::WireFrame:
    {
      Point rstart(corner);
      for (int i=0;i<nu;i++)
      {
        Point p1(rstart);
        glBegin(GL_LINE_STRIP);
        for (int j=0;j<nv;j++)
        {
          Point pp1(p1+w*verts(i, j));
          if (have_matls)
            di->set_material(matls(i, j).get_rep());
          if (have_normals)
          {
            Vector normal(normals(i, j));
            glNormal3d(normal.x(), normal.y(), normal.z());
          }
          glVertex3d(pp1.x(), pp1.y(), pp1.z());

          p1+=vv;
        }
        glEnd();
        rstart+=uu;
      }
      rstart=corner;
      for (int j=0;j<nv;j++)
      {
        Point p1(rstart);
        glBegin(GL_LINE_STRIP);
        for (int i=0;i<nu;i++)
        {
          Point pp1(p1+w*verts(i, j));
          if (have_matls)
            di->set_material(matls(i, j).get_rep());
          if (have_normals)
          {
            Vector normal(normals(i, j));
            glNormal3d(normal.x(), normal.y(), normal.z());
          }
          glVertex3d(pp1.x(), pp1.y(), pp1.z());

          p1+=uu;
        }
        glEnd();
        rstart+=vv;
      }
    }
    break;
  case DrawInfoOpenGL::Flat:
  case DrawInfoOpenGL::Gouraud:
    {
      if (have_matls)
        di->set_material(matls(0,0).get_rep());
      Point rstart(corner);
      if (have_normals && have_matls)
      {
        for (int i=0;i<nu-1;i++)
        {
          Point p1(rstart);
          Point p2(rstart+uu);
          rstart=p2;
          glBegin(GL_QUAD_STRIP);
          for (int j=0;j<nv;j++)
          {
            Point pp1(p1+w*verts(i, j));
            Point pp2(p2+w*verts(i+1, j));
            float c[4];
            matls(i,j)->diffuse.get_color(c);
            glColor3fv(c);
            Vector& normal = normals(i, j);
            glNormal3d(normal.x(), normal.y(), normal.z());
            glVertex3d(pp1.x(), pp1.y(), pp1.z());

            matls(i+1, j)->diffuse.get_color(c);
            glColor3fv(c);
            Vector& normal2 = normals(i+1, j);
            glNormal3d(normal2.x(), normal2.y(), normal2.z());
            glVertex3d(pp2.x(), pp2.y(), pp2.z());
            p1+=vv;
            p2+=vv;
          }
          glEnd();
        }
      }
      else if (have_matls)
      {
        glNormal3d(w.x(), w.y(), w.z());
        for (int i=0;i<nu-1;i++)
        {
          Point p1(rstart);
          Point p2(rstart+uu);
          rstart=p2;
          glBegin(GL_QUAD_STRIP);
          for (int j=0;j<nv;j++)
          {
            Point pp1(p1+w*verts(i, j));
            Point pp2(p2+w*verts(i+1, j));
            float c[4];
            matls(i,j)->diffuse.get_color(c);
            glColor3fv(c);
            glVertex3d(pp1.x(), pp1.y(), pp1.z());

            matls(i+1, j)->diffuse.get_color(c);
            glColor3fv(c);
            glVertex3d(pp2.x(), pp2.y(), pp2.z());
            p1+=vv;
            p2+=vv;
          }
          glEnd();
        }
      }
      else if (have_normals)
      {
        for (int i=0;i<nu-1;i++)
        {
          Point p1(rstart);
          Point p2(rstart+uu);
          rstart=p2;
          glBegin(GL_QUAD_STRIP);
          for (int j=0;j<nv;j++)
          {
            Point pp1(p1+w*verts(i, j));
            Point pp2(p2+w*verts(i+1, j));
            Vector& normal = normals(i, j);
            glNormal3d(normal.x(), normal.y(), normal.z());
            glVertex3d(pp1.x(), pp1.y(), pp1.z());

            Vector& normal2 = normals(i+1, j);
            glNormal3d(normal2.x(), normal2.y(), normal2.z());
            glVertex3d(pp2.x(), pp2.y(), pp2.z());
            p1+=vv;
            p2+=vv;
          }
          glEnd();
        }
      }
      else
      {
        glNormal3d(w.x(), w.y(), w.z());
        for (int i=0;i<nu-1;i++)
        {
          Point p1(rstart);
          Point p2(rstart+uu);
          rstart=p2;
          glBegin(GL_QUAD_STRIP);
          for (int j=0;j<nv;j++)
          {
            Point pp1(p1+w*verts(i, j));
            Point pp2(p2+w*verts(i+1, j));
            glVertex3d(pp1.x(), pp1.y(), pp1.z());

            glVertex3d(pp2.x(), pp2.y(), pp2.z());
            p1+=vv;
            p2+=vv;
          }
          glEnd();
        }
      }
    }
    break;
  }
  post_draw(di);
}

} // End namespace SCIRun

