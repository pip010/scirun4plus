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


#ifndef CORE_GEOMETRY_SLIVRCONVERT_H
#define CORE_GEOMETRY_SLIVRCONVERT_H 1

#include <slivr/Transform.h>
#include <slivr/Vector.h>
#include <slivr/Point.h>
#include <slivr/Plane.h>
#include <slivr/BBox.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Plane.h>
#include <Core/Geometry/Transform.h>
#include <Core/Geometry/BBox.h>

#include <Core/Geometry/share.h>

// These inline functions cast from SLIVR geometry to SCIRun geometry

// Need to generate more elegant versions
// of these converters

inline SLIVR::Transform TO_SLIVR_TRANSFORM(const SCIRun::Transform& t)
{
  SLIVR::Transform st;
  for (int i=0;i<4;i++) for (int j=0;j<4;j++)
    st.set_mat_val(i,j,t.get_mat_val(i,j));
  return (st);
}

inline SCIRun::Transform FROM_SLIVR_TRANSFORM(const SLIVR::Transform& t)
{
  SCIRun::Transform st;
  for (int i=0;i<4;i++) for (int j=0;j<4;j++)
    st.set_mat_val(i,j,t.get_mat_val(i,j));
  return (st);
}

inline SLIVR::Point TO_SLIVR_POINT(const SCIRun::Point& p)
{
  return (SLIVR::Point(p.x(),p.y(),p.z()));
}

inline SCIRun::Point FROM_SLIVR_POINT(const SLIVR::Point& p)
{
  return (SCIRun::Point(p.x(),p.y(),p.z()));
}

inline SLIVR::Vector TO_SLIVR_VECTOR(const SCIRun::Vector& p)
{
  return (SLIVR::Vector(p.x(),p.y(),p.z()));
}

inline SCIRun::Vector FROM_SLIVR_VECTOR(const SLIVR::Vector& p)
{
  return (SCIRun::Vector(p.x(),p.y(),p.z()));
}

inline SCIRun::BBox FROM_SLIVR_BBOX(const SLIVR::BBox& b)
{
  return (SCIRun::BBox(FROM_SLIVR_POINT(b.min()),
                      FROM_SLIVR_POINT(b.max())));
}

inline SLIVR::BBox TO_SLIVR_BBOX(const SCIRun::BBox& b)
{
  return (SLIVR::BBox(TO_SLIVR_POINT(b.min()),
                     TO_SLIVR_POINT(b.max())));
}

#endif
