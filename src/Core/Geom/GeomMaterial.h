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
 *  Material.h:  Material Properities for Geometry
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#ifndef CORE_GEOM_GEOMMATERIAL_H
#define CORE_GEOM_GEOMMATERIAL_H 1

//! Get the base class
#include <Core/Datatypes/Material.h>

//! A container for including Material into SceneGraph
#include <Core/Geom/GeomContainer.h>

//! For Windows
#include <Core/Geom/share.h>

namespace SCIRun {

class SCISHARE GeomMaterial : public GeomContainer
{
  public:
    GeomMaterial(GeomHandle, const MaterialHandle&);
    GeomMaterial(const GeomMaterial&);
    void setMaterial(const MaterialHandle&);
    MaterialHandle getMaterial();
    virtual ~GeomMaterial();
    virtual GeomObj* clone();

    virtual void draw(DrawInfoOpenGL*, Material*, double time);
  private:
    MaterialHandle material_;
};

} // End namespace SCIRun

#endif
