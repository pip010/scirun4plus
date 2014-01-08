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
 *  TexSquare.cc: ?
 *
 *  Written by:
 *   Author: ?
 *   Department of Computer Science
 *   University of Utah
 *   Date: ?
 *
 */

#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/TexSquare.h>

#include <string.h>

namespace SCIRun {

TexSquare::TexSquare() :
  GeomObj(),
  normal_(1.0, 0.0, 0.0),
  texture(0),
  numcolors(0),
  width_(2),
  height_(2),
  texname_(0),
  alpha_cutoff_(0.0)
{
}

TexSquare::TexSquare( const TexSquare &copy ) : GeomObj(copy)
{
}

TexSquare::~TexSquare()
{
}

void
TexSquare::set_coords(float *tex, float *coords)
{
  memcpy(tex_coords_, tex, 8*sizeof(float));
  memcpy(pos_coords_, coords, 12*sizeof(float));
}

void 
TexSquare::set_texture( unsigned char *tex, int num, int w, int h)
{
  width_ = w;
  height_ = h;
  numcolors = num;
  const int count = numcolors*width_*height_;
  texture = new unsigned char[count];
  memcpy(texture, tex, count);
}

void
TexSquare::set_texname(unsigned int texname)
{
  texname_ = texname;
}

void
TexSquare::set_normal(Vector &normal)
{
  normal_ = normal;
}

void
TexSquare::set_alpha_cutoff(double alpha)
{
  alpha_cutoff_ = alpha;
}

GeomObj* 
TexSquare::clone()
{
  return new TexSquare( *this );
}

void 
TexSquare::get_bounds( BBox& bb )
{
  for (int i = 0; i < 4; ++i) 
  {
    bb.extend(Point(pos_coords_[i*3+0],pos_coords_[i*3+1],pos_coords_[i*3+2]));
  }
}

void
TexSquare::draw(DrawInfoOpenGL* di, Material* matl, double)
{
#ifdef _WIN32
  // unwind the GL error stack.  In BioImage, the check for glGetError below
  // is blaming some error on the glBindTexture (or other tex ops)
  // for some previous error, and if the textures executed correctly,
  // their operations need to be executed
  // Ideally, the errors would be caught where they are thrown.....
  GLenum err;
  while((err = glGetError()) != GL_NO_ERROR);
#endif
  double old_ambient = di->ambient_scale_;
  di->ambient_scale_ = 15;
  if (!pre_draw(di, matl, 1)) {
    di->ambient_scale_ = old_ambient;
    return;
  }
  glEnable(GL_TEXTURE_2D);
  bool bound = glIsTexture(texname_);
  if (!bound)
    glGenTextures(1, (GLuint *)&texname_);
  glBindTexture(GL_TEXTURE_2D, texname_);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  if (!bound && texture)
  {
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelTransferi(GL_MAP_COLOR,0);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width_, height_, 0,
                 GL_RGBA, GL_UNSIGNED_BYTE, texture);
  }
  if (GL_NO_ERROR == glGetError())
  {
    glAlphaFunc(GL_GEQUAL, alpha_cutoff_);
    glEnable(GL_ALPHA_TEST);
    glDisable(GL_BLEND);
    glColor4d(1., 1., 1., 1.);
    glBegin( GL_QUADS );
    for (int i = 0; i < 4; i++)
    {
      glNormal3d(normal_.x(), normal_.y(), normal_.z());
      glTexCoord2fv(tex_coords_+i*2);
      glVertex3fv(pos_coords_+i*3);
    }
    glEnd();
    glDisable(GL_ALPHA_TEST);
  }
  glBindTexture(GL_TEXTURE_2D, 0);
  glDisable(GL_TEXTURE_2D);
  di->ambient_scale_ = old_ambient;
  post_draw(di);
}

} // End namespace SCIRun

