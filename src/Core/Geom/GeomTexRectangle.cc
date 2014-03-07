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


#include <Core/Util/Debug.h>

#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomTexRectangle.h>

#include <slivr/ShaderProgramARB.h>

#include <cstring>


namespace SCIRun {

GeomTexRectangle::GeomTexRectangle()
  : GeomObj(),
    texture_(0),
    sId_(0),
    fId_(0),
    numcolors_(0),
    width_(2),
    height_(2),
    texname_(0),
    alpha_cutoff_(0.0),
    interp_(false),
    trans_(false),
    use_normal_(false),
    shader_(0),
    fog_shader_(0)
{
  DEBUG_CONSTRUCTOR("GeomTexRectangle")
}

GeomTexRectangle::GeomTexRectangle( const GeomTexRectangle &copy ) : GeomObj(copy) 
{
  DEBUG_CONSTRUCTOR("GeomTexRectangle")
}

GeomTexRectangle::~GeomTexRectangle()
{
  if(shader_) delete shader_;
  if(fog_shader_) delete fog_shader_;
  
  if (texture_) delete texture_;
  DEBUG_DESTRUCTOR("GeomTexRectangle")
}

void
GeomTexRectangle::set_coords(float *tex, float *coords)
{
  memcpy(tex_coords_, tex, 8*sizeof(float));
  memcpy(pos_coords_, coords, 12*sizeof(float));
}


void 
GeomTexRectangle::set_texture( unsigned char *tex, int num, int w, int h) 
{
  width_ = w;
  height_ = h;
  numcolors_ = num;
  const int count = numcolors_*width_*height_;
  if( texture_ ) delete [] texture_;
  texture_ = new unsigned char[count];
  memcpy(texture_, tex, count);
}

void
GeomTexRectangle::set_texname(unsigned int texname) 
{
  texname_ = texname;
}

void GeomTexRectangle::set_normal( float *norm )
{
  memcpy(normal_, norm, 3*sizeof(float));
  use_normal_ = true;
}

void
GeomTexRectangle::set_alpha_cutoff(double alpha) 
{
  alpha_cutoff_ = alpha;
}


GeomObj* 
GeomTexRectangle::clone() 
{
  return new GeomTexRectangle( *this );
}

void 
GeomTexRectangle::get_bounds( BBox& bb ) 
{
  for (int i = 0; i < 4; ++i)
  {
    bb.extend(Point(pos_coords_[i*3+0],pos_coords_[i*3+1],pos_coords_[i*3+2]));
  }
}

#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
#define FRAG \
"!!ARBfp1.0 \n" \
"ATTRIB t = fragment.texcoord[1]; \n" \
"TEMP c; \n" \
"ATTRIB cf = fragment.color; \n" \
"TEX c, t, texture[1], 2D; \n" \
"MUL c, c, cf; \n" \
"MOV result.color, c; \n" \
"END"

#define FOG \
"!!ARBfp1.0 \n" \
"PARAM fc = state.fog.color; \n" \
"PARAM fp = state.fog.params; \n" \
"TEMP fctmp; \n" \
"ATTRIB t = fragment.texcoord[1]; \n" \
"ATTRIB tf = fragment.texcoord[2]; \n" \
"TEMP c; \n" \
"TEMP v; \n" \
"ATTRIB cf = fragment.color; \n" \
"TEX c, t, texture[1], 2D; \n" \
"MUL c, c, cf; \n" \
"SUB v.x, fp.z, tf.x; \n" \
"MUL_SAT v.x, v.x, fp.w; \n" \
"MUL fctmp, c.w, fc; \n" \
"LRP c.xyz, v.x, c.xyzz, fctmp.xyzz; \n" \
"MOV result.color, c; \n" \
"END"
#endif


void
GeomTexRectangle::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;
  GLboolean use_fog = glIsEnabled(GL_FOG);
#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
  if ( !shader_ || !fog_shader_ )
  {
    shader_ = new SLIVR::FragmentProgramARB( FRAG );
    shader_->create();
    fog_shader_ = new SLIVR::FragmentProgramARB( FOG );
    fog_shader_->create();
  }

  if (use_fog)
  {
#ifdef _WIN32
    if (glActiveTexture)
#endif
    {
      // enable texture unit 2 for fog
      glActiveTexture(GL_TEXTURE2);
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
      glEnable(GL_TEXTURE_3D);
    }
  }

#ifdef _WIN32
  if (glActiveTexture)
#endif
    glActiveTexture(GL_TEXTURE1);
#endif

  bool bound = glIsTexture(texname_);

  if (!bound)
  {
    glGenTextures(1, (GLuint *)&texname_);
  }

  glBindTexture(GL_TEXTURE_2D, texname_);

  if (!bound)
  {
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width_, height_, 0,
                 GL_RGBA, GL_UNSIGNED_BYTE, texture_ );
  }

  if (GL_NO_ERROR == glGetError())
  {
    glEnable(GL_TEXTURE_2D);

    if (interp_)
    {
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    }
    else
    {
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    }

#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
    if (use_fog)
    {
      fog_shader_->bind();
    }
    else
    {
      shader_->bind();
    }
#endif

    GLfloat ones[4] = {1.0, 1.0, 1.0, 1.0};
    glColor4fv(ones);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glDepthMask(GL_TRUE);

    if (trans_)
    {
#if !defined(GL_ARB_fragment_program) && !defined(GL_ATI_fragment_shader)
      glAlphaFunc(GL_GREATER, alpha_cutoff_);
      glEnable(GL_ALPHA_TEST);
#endif
      glEnable(GL_BLEND);

      // Workaround for old bad nvidia headers.
#if defined(GL_FUNC_ADD)
      glBlendEquation(GL_FUNC_ADD);
#else
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
#endif
    }

    float mvmat[16];
    if (use_fog)
    {
      glGetFloatv(GL_MODELVIEW_MATRIX, mvmat);
    }

    glBegin( GL_QUADS );
    if ( use_normal_ )
    {
      glNormal3fv( normal_ );
    }

    for (int i = 0; i < 4; i++)
    {
#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
#  ifdef _WIN32
      if (glMultiTexCoord2fv && glMultiTexCoord3f)
      {
#  endif
        glMultiTexCoord2fv(GL_TEXTURE1,tex_coords_+i*2);
        if (use_fog)
        {
          float *pos = pos_coords_+i*3;
          float vz = mvmat[2]* pos[0]
            + mvmat[6]*pos[1]
            + mvmat[10]*pos[2] + mvmat[14];
          glMultiTexCoord3f(GL_TEXTURE2, -vz, 0.0, 0.0);
        }
#  ifdef _WIN32
      }
      else
      {
        glTexCoord2fv(tex_coords_+i*2);
      }
#  endif // _WIN32
#else
      glTexCoord2fv(tex_coords_+i*2);
#endif
      glVertex3fv(pos_coords_+i*3);
    }
    glEnd();

    glFlush();
    if ( trans_ )
    {
      glDisable(GL_ALPHA_TEST);
    }
    glDisable(GL_BLEND);
    glDisable(GL_TEXTURE_2D);

    glBindTexture(GL_TEXTURE_2D, 0);
#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
    if ( use_fog )
    {
      fog_shader_->release();
#ifdef _WIN32
      if (glActiveTexture)
#endif
        glActiveTexture(GL_TEXTURE2);
      glDisable(GL_TEXTURE_3D);
    }
    else
    {
      shader_->release();
    }
#endif

  }
  else
  {
    std::cerr<<"Some sort of texturing error\n";
  }

#if defined(GL_ARB_fragment_program) || defined(GL_ATI_fragment_shader)
#ifdef _WIN32
  if (glActiveTexture)
#endif
    glActiveTexture(GL_TEXTURE0);
#endif
  di->polycount_++;
  post_draw(di);
}

} // End namespace SCIRun

