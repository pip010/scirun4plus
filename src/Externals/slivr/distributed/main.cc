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

#include <OGLContext.h>
#include <slivr/Texture.h>
#include <slivr/ShaderProgramARB.h>
#include <slivr/VideoCardInfo.h>
#include <slivr/VolumeRenderer.h>

using namespace std;
using namespace SLIVR_DIST;
using namespace SLIVR;

bool
compute_depth(const Point& veyep, 
              const Point &vlookat, 
              double& znear, double& zfar,
              const BBox &bb)
{
  znear = DBL_MAX;
  zfar =- DBL_MAX;
  //  BBox bb = input_bounds_;
  if (bb.valid())
  {
    // We have something to draw.
    Point min(bb.min());
    Point max(bb.max());
    Point eyep(veyep);
    Vector dir(vlookat-eyep);
    const double dirlen2 = dir.length2();
    if (dirlen2 < 1.0e-6 || dirlen2 != dirlen2)
    {
      return false;
    }
    dir.safe_normalize();
    const double d = -Dot(eyep, dir);
    for (int i=0;i<8;i++)
    {
      const Point p((i&1)?max.x():min.x(),
                    (i&2)?max.y():min.y(),
                    (i&4)?max.z():min.z());
      const double dist = Dot(p, dir) + d;
      znear = Min(znear, dist);
      zfar = Max(zfar, dist);
    }
    znear *= 0.99;
    zfar  *= 1.01;

    if (znear <= 0.0)
    {
      if (zfar <= 0.0)
      {
        // Everything is behind us - it doesn't matter what we do.
        znear = 1.0;
        zfar = 2.0;
      }
      else
      {
        znear = zfar * 0.001;
      }
    }
    return true;
  }
  return false;
}

inline 
double 
degrees_to_radians(double d) {
  return d*(M_PI / 180.);
}


bool
viewall(const BBox &bbox)
{
  //BBox bbox = input_bounds_;

  Point lookat(bbox.center());

  // Move forward/backwards until entire view is in scene.
  // change this a little, make it so that the FOV must be 20 deg.
  double myfov = 20.0;

  Vector diag(bbox.diagonal());
    
  double w = diag.length();
  if( w < 0.000001 ){
    BBox bb;
    bb.reset();
    Vector epsilon(0.001, 0.001, 0.001 );
    bb.extend( bbox.min() - epsilon );
    bb.extend( bbox.max() + epsilon );
    w = bb.diagonal().length();
  }
    
  Point home_eyep(-0.5, -1.5, 0);
  Point home_lookat(0.0, 0.0, 0.0);
  Vector home_up(0.0, 0.0, 1.0);
  Vector lookdir(home_lookat - home_eyep); 
  lookdir.safe_normalize();
  double r = degrees_to_radians(myfov / 2.0);
  const double scale = 1.0 / (2. * tan(r));
  double length = w * scale;
  Point eyep(lookat - lookdir * length);

  Plane upplane(eyep, lookdir);
  Vector up(upplane.project(home_up));


  gluLookAt(eyep.x(), eyep.y(), eyep.z(),
            lookat.x(), lookat.y(), lookat.z(),
            up.x(), up.y(), up.z());

  cerr << lookat << endl;

  return true;

}

SLIVR::VolumeRenderer*
init(const BBox &bounds, Texture *tex)
{
  ShaderProgramARB::init_shaders_supported();

  float w = 256;
  float h = 256;
  //cerr << "w: " << w << " h: " << h << endl;
  glViewport(0, 0, (GLint)w, (GLint)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  Point home_eyep(22.1, 17.6, 601.5);
  Point home_lookat(50.0, 50.0, 50.0);
  Vector home_up(0.0, 1.0, 0.0);
  double near, far;
  if (compute_depth(home_eyep, home_lookat, near, far, bounds)) {
    gluPerspective(20, w / h, near, far);
  }
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  viewall(bounds);


  float rgba[256*4];
  double spacing = 360.0/255.0;
  double hue = 0;
  double dhue = spacing;
  for (int i = 0; i < 256; ++i) {
    
    Color col(HSVColor(hue, 1.0, 1.0));
    //                       (i % 3) / 3.0 + 1.0/3.0, 
    //(i % 4) / 4.0 + 1.0/4.0 ));
    hue += dhue;
    
    rgba[i*4+0] = col.r();
    rgba[i*4+1] = col.g();
    rgba[i*4+2] = col.b();
    if (i < 90) {
      rgba[i*4+3] = 0.0;
    } else {
      rgba[i*4+3] = 0.052;
    }
  }
  ColorMap *cm =  new ColorMap(rgba);
  vector<Plane*> planes;
  // Create a colormap2.
  //cmap2_.push_back(create_cmap2());
  
  int vcm = video_card_memory_size();
  vector<ColorMap2*> empty(0);
  VolumeRenderer *vr = 
    new VolumeRenderer(tex, cm, empty, planes, vcm*1024*1024);
  vr->set_sampling_rate(4.4);
  vr->set_material(0.5, 0.5, .388, 24);
  vr->set_slice_alpha(0.0);

  return vr;
}



void 
draw(VolumeRenderer *vr)
{
  /* clear color and depth buffers */
  //glDrawBuffer(GL_BACK);
  glDisable(GL_DEPTH_TEST);
  //glClearColor(0.0, 0.0, 0.0, 0.0);
  //glClearDepth(1.0); 
  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  vr->draw(false, false);

  glFlush();

}

SLIVR::Texture*
build_texture(string in, BBox &inbounds)
{
  Nrrd *nin = nrrdNew();

  /* read in the nrrd from file */
  char *err;
  if (nrrdLoad(nin, in.c_str(), NULL)) {
    err = biffGetDone(NRRD);
    fprintf(stderr, "%s: trouble reading \"%s\":\n%s", 
	    "voldist texture", in.c_str(), err);
    free(err);
    return 0;
  }

  if (nin->dim != 3) {
    cerr << "Error must input a Nrrd of dim == 3" << endl;
    return 0;
  }
  
  Texture *tex = new Texture();
  cerr << "video_card_memory_size():" << video_card_memory_size() << endl;
  tex->build(nin, 0, 0, 256, 0, 0, video_card_memory_size());

  BBox bounds;
  Point max(nin->axis[0].size, nin->axis[1].size, nin->axis[2].size);
  bounds.extend(Point(0.0, 0.0,0.0));
  bounds.extend(max);
  //cerr << bounds << endl;
  inbounds = bounds;

  return tex;
}



int main(int argc, char *argv[])
{

  if (argc != 3) {
    cerr << "Usage: voldist <path to data nrrd> <path to viewing matrix>" 
         << endl;
    return 0;
  }

  string inrrd_fn = argv[1];
  string view_fn = argv[2];

  BBox bb;
  Texture *tex = build_texture(inrrd_fn, bb);
  if (! tex) {
    cerr << "Could not build texture from: " << inrrd_fn << " , exiting." 
         << endl;
    exit(1);
  }

  OGLContext ogl;
  ogl.create();
  ogl.make_current();
  VolumeRenderer *vr = init(bb, tex);


  XEvent event;

 
  XSelectInput(ogl.display(), ogl.window(), 0xffff);
  XMapWindow(ogl.display(), ogl.window());
  XFlush(ogl.display());
  glXWaitX();

//   if (!tex_) {
//     build_texture();
//   }
//   if (!init_done_) {
//     init();
//   }


  while (1)
  {
    ogl.draw();
    draw(vr);
    XNextEvent(ogl.display(), &event);
    cerr << "Event recieved is " << event.type << endl;
  }

}
