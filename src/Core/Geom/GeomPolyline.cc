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
 *  GeomPolyline.cc: Polyline object
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#include <Core/Util/Debug.h>

#include <Core/Geom/GeomPolyline.h>
#include <Core/Geom/GeomLine.h>
#include <Core/Geom/DrawInfoOpenGL.h>

namespace SCIRun {

GeomPolyline::GeomPolyline()
{
  DEBUG_CONSTRUCTOR("GeomPolyline")
}

GeomPolyline::GeomPolyline(const GeomPolyline& copy)
: GeomVertexPrim(copy)
{
  DEBUG_CONSTRUCTOR("GeomPolyline")
}

GeomPolyline::~GeomPolyline() 
{
  DEBUG_DESTRUCTOR("GeomPolyline")
}

GeomObj* GeomPolyline::clone()
{
  return new GeomPolyline(*this);
}

void
GeomPolyline::draw(DrawInfoOpenGL* di, Material* matl, double currenttime)
{
  if (!pre_draw(di, matl, 0)) return;
  di->polycount_+=verts.size()-1;
  glBegin(GL_LINE_STRIP);
  if (times.size() == verts.size())
  {
    for (int i=0;i<verts.size() && currenttime >= times[i];i++)
    {
      verts[i]->emit_all(di);
    }
  }
  else
  {
    for (int i=0;i<verts.size();i++)
    {
      verts[i]->emit_all(di);
    }
  }
  glEnd();
  post_draw(di);
}

GeomPolylineTC::GeomPolylineTC(int drawmode, double drawdist)
: drawmode(drawmode), drawdist(drawdist)
{
}

GeomPolylineTC::~GeomPolylineTC()
{
}

void GeomPolylineTC::add(double t, const Point& p, const Color& c)
{
  int s=data.size();
  data.grow(7);
  data[s]=t;
  data[s+1]=c.r();
  data[s+2]=c.g();
  data[s+3]=c.b();
  data[s+4]=p.x();
  data[s+5]=p.y();
  data[s+6]=p.z();
  bbox.extend(p);
}

GeomPolylineTC::GeomPolylineTC(const GeomPolylineTC& copy)
  : GeomObj(copy),
    data(copy.data)
{
}

GeomObj* GeomPolylineTC::clone()
{
  return new GeomPolylineTC(*this);
}

void GeomPolylineTC::get_bounds(BBox& box)
{
  box.extend(bbox);
}

void
GeomPolylineTC::draw(DrawInfoOpenGL* di, Material* matl, double currenttime)
{
  if (data.size() == 0)
    return;
  if (!pre_draw(di, matl, 0)) return;
  float* d=&data[0];
  float* dend=d+data.size();
  if (drawmode < 1 || drawmode > 3)
  {
   std:: cerr << "Bad drawmode: " << drawmode << std::endl;
  }
  if (drawmode==1)
  {
    glBegin(GL_LINE_STRIP);
    while (d<dend && *d <= currenttime)
    {
      glColor3fv(d+1);
      glVertex3fv(d+4);
      d+=7;
    }
    di->polycount_+=(d-&data[0])/7-1;
    glEnd();
  }
  else
  {
    // Find the start and end points.
    int n=(dend-d)/7;
    int l=0;
    int h=n-1;
    while (l<h-1)
    {
      int m=(l+h)/2;
      if (currenttime < d[7*m])
      {
        h=m;
      }
      else
      {
        l=m;
      }
    }
    int iend=l;
    l=0;
    // Leave h - it still bounds us on the top
    double begtime=Max(0.0, currenttime-drawdist);
    while (l<h-1)
    {
      int m=(l+h)/2;
      if (begtime < d[7*m])
      {
        h=m;
      }
      else
      {
        l=m;
      }
    }
    int istart=l;
    if (istart==iend) return;
    d=&data[7*istart];
    dend=&data[7*iend]+7;
    di->polycount_+=(dend-d)/7-1;
    glBegin(GL_LINE_STRIP);
    if (drawmode == 2)
    {
      while (d<dend)
      {
        glColor3fv(d+1);
        glVertex3fv(d+4);
        d+=7;
      }
    }
    else if (drawmode == 3)
    {
      while (d<dend)
      {
        float s=(*d-begtime)/drawdist;
        glColor3f(d[1]*s, d[2]*s, d[3]*s);
        glVertex3fv(d+4);
        d+=7;
      }
    }
    glEnd();
  }
  post_draw(di);
}

} // End namespace SCIRun

