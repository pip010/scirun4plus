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
 *  GeomColorMap.cc: Set colormap for indexed color primitives.
 *
 *  Written by:
 *   Michael Callahan
 *   Department of Computer Science
 *   University of Utah
 *   August 2003
 *
 */

#include <Core/Util/Debug.h>

#include <Core/Geom/GeomColorMap.h>
#include <Core/Geom/DrawInfoOpenGL.h>

namespace SCIRun {

GeomColorMap::GeomColorMap(GeomHandle obj, ColorMapHandle cmap)
  : GeomContainer(obj), cmap_(cmap)
{
  DEBUG_CONSTRUCTOR("GeomColorMap")
}

GeomColorMap::GeomColorMap(const GeomColorMap &copy)
  : GeomContainer(copy), cmap_(copy.cmap_)
{
  DEBUG_CONSTRUCTOR("GeomColorMap")
}

GeomColorMap::~GeomColorMap()
{
  DEBUG_DESTRUCTOR("GeomColorMap")
}

GeomObj*
GeomColorMap::clone()
{
  return new GeomColorMap(*this);
}

void
GeomColorMap::draw(DrawInfoOpenGL* di, Material *m, double time)
{
  if ( !pre_draw(di, m, 0) ) return;

  if (!cmap_.get_rep())
  {
    child_->draw(di, m, time);
  }
  else
  {
    unsigned int old_cmtexture = di->cmtexture_;

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, (GLuint *)&(di->cmtexture_));
    
    // Send Cmap
    glBindTexture(GL_TEXTURE_1D, di->cmtexture_);
    glTexImage1D(GL_TEXTURE_1D, 0, 4, 256, 0, GL_RGBA, GL_FLOAT,
                 cmap_->get_rgba());

#ifdef GL_CLAMP_TO_EDGE
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
#else
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
#endif
    if (cmap_->resolution() == 256)
    {
      glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    }
    else
    {
      glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    }

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    // Do Cmap transform to min-max
    glMatrixMode(GL_TEXTURE);
    glPushMatrix();

    const double r = cmap_->resolution() / 256.0;
    glScaled(r / (cmap_->getMax() - cmap_->getMin()), 1.0, 1.0);
    glTranslated(-cmap_->getMin(), 0.0, 0.0);

    glMatrixMode(GL_MODELVIEW);

    // Draw child
    di->using_cmtexture_++;
    child_->draw(di,m,time);
    di->using_cmtexture_--;

    glMatrixMode(GL_TEXTURE);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    
    glDeleteTextures(1,(GLuint *)&(di->cmtexture_));
    di->cmtexture_ = old_cmtexture;
  }

  post_draw(di);
}

} // End namespace SCIRun


