//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : OGLCanvas.cc
//    Author : Martin Cole
//    Date   : Wed Jul 11 13:47:51 2007
//    Modified from James Bigler's wxtest code.


#include <slivr/Color.h>
#include <slivr/VolumeRenderer.h>
#include <slivr/VideoCardInfo.h>
#include <slivr/ShaderProgramARB.h>
#include <slivr/BBox.h>
#include <slivr/ColorMap.h>
#include <slivr/ColorMap2.h>
#include <slivr/CM2Widget.h>
#include <slivr/Plane.h>
#include <OGLCanvas.h>
#include <vector>
#include <iostream>

#ifdef __APPLE__
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

namespace SLIVR_EXAMPLES {

using namespace std;
//needed to make mac happy with the included forward declares...
using SLIVR::HSVColor; 
using SLIVR::Point;
using namespace SLIVR;

BEGIN_EVENT_TABLE(OGLCanvas, wxGLCanvas)
  EVT_PAINT(OGLCanvas::draw)
END_EVENT_TABLE()

const char*s = "OGLCanvas";
const wxString can_name(&s[0], &s[7]);


OGLCanvas::OGLCanvas(wxWindow* parent, int* attribs) : 
  wxGLCanvas(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, 
	     0, can_name, attribs), 
  init_done_(false),
  bounds_set_(false),
  use_2d_cmap_(false),
  last_use_2d_cmap_(false),
  tex_(0),
  vr_(0),
  spin_p_(false),
  window_(parent)
{
  //  this->Connect(wxID_ANY, wxEVT_MOTION,
  //		wxMoveEventHandler(OGLCanvas::OnMouseMove));
}

void OGLCanvas::draw(wxPaintEvent& event)
{
  wxPaintDC dc(this);
  if (! window_->IsShownOnScreen()) return;

  SetCurrent();  

  if (!tex_) {
    build_texture();
  }
  if (!init_done_) {
    init();
  }

  /* clear color and depth buffers */
  glDrawBuffer(GL_BACK);
  glDisable(GL_DEPTH_TEST);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClearDepth(1.0); 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (!vr_) {
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
    cm_ =  new ColorMap(rgba);

    // Create a colormap2.
    cmap2_.push_back(create_cmap2());

    int vcm = video_card_memory_size();
    vector<ColorMap2*> empty(0);
    vr_ = new VolumeRenderer(tex_, cm_, empty, planes_, vcm*1024*1024);
    vr_->set_sampling_rate(4.4);
    vr_->set_material(0.5, 0.5, .388, 24);
    vr_->set_slice_alpha(0.0);
  }
  
  if (last_use_2d_cmap_ != use_2d_cmap_) {
    if (use_2d_cmap_) {
      vr_->set_colormap1(0);
      vr_->set_colormap2(cmap2_);
    } else {
      vr_->set_colormap2_dirty();
      vector<ColorMap2*> empty(0);
      vr_->set_colormap2(empty);
      vr_->set_colormap1(cm_);
    }
    last_use_2d_cmap_ = use_2d_cmap_;
  }

  vr_->draw(false, false);

  glFlush();
  SwapBuffers();

  if (spin_p()) {
    glMatrixMode(GL_MODELVIEW);
    glTranslatef(51.5, 47, 80.5);
    glRotatef(13.0, 0, 0, 1); 
    glTranslatef(-51.5, -47, -80.5);
  }

}

#ifndef CALLBACK
#define CALLBACK
#endif


void CALLBACK errorCallback(GLenum errorCode)
{
  const GLubyte* estring;

  estring = gluErrorString(errorCode);
  wxLogError(wxT("Quadric Error: %s\n"), estring);
}

void OGLCanvas::init()
{
  SetCurrent();
  // Must call init_shaders_supported, it also inits glew.
  if (! init_done_) {
    ShaderProgramARB::init_shaders_supported();
  }

  float w = GetSize().x;
  float h = GetSize().y;
  //cerr << "w: " << w << " h: " << h << endl;
  glViewport(0, 0, (GLint)w, (GLint)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  Point home_eyep(22.1, 17.6, 601.5);
  Point home_lookat(50.0, 50.0, 50.0);
  Vector home_up(0.0, 1.0, 0.0);
  double near, far;
  if (compute_depth(home_eyep, home_lookat, near, far)) {
    gluPerspective(20, w / h, near, far);
  }
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  viewall();
  init_done_ = true;
}


ColorMap2*
OGLCanvas::create_cmap2()
{
  static ColorMap2 *cm2 = 0;
  if (cm2) return cm2;

  cm2 = new ColorMap2();

  RectangleCM2Widget *r = new RectangleCM2Widget(CM2_RECTANGLE_1D, 0.116191,
						 0.34375, 0.311544, 0.34375,
						 0.474659);
  r->set_color(Color(.8, .17, .17));
  r->set_alpha(1.0);
  cm2->widgets().push_back(r);

  TriangleCM2Widget *t = new TriangleCM2Widget(0.482499, 0.0573953, 0.422751,
					       -0.181896, 0.453518);
  
  t->set_color(Color(.8, .8, .33));
  t->set_alpha(0.6);
  cm2->widgets().push_back(t);


  r = new RectangleCM2Widget(CM2_RECTANGLE_1D, 0.664062, 0.167969, 0.201172,
			     0.246094, 0.515559);  
  r->set_color(Color(.4, .4, 1.0));
  r->set_alpha(0.8);
  cm2->widgets().push_back(r);

  r = new RectangleCM2Widget(CM2_RECTANGLE_1D, 0.470703, 0.535156, 0.378906,
			     0.328125, 0.497439);
  r->set_color(Color(1.0, 1.0, 1.0));
  r->set_alpha(0.8);
  cm2->widgets().push_back(r);

  return cm2;
}


bool
OGLCanvas::compute_depth(const Point& veyep, 
			 const Point &vlookat, 
			 double& znear, double& zfar)
{
  znear = DBL_MAX;
  zfar =- DBL_MAX;
  BBox bb = input_bounds_;
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
  else
  {
    return false;
  }
}

inline 
double 
degrees_to_radians(double d) {
  return d*(M_PI / 180.);
}


bool
OGLCanvas::viewall()
{
  BBox bbox = input_bounds_;
  if (bounds_set_)
  {

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
  return false;
}
#define DATAPATH(x) 
void
OGLCanvas::build_texture()
{
  Nrrd *nin = nrrdNew();

  /* read in the nrrd from file */
  string fn = EX_DATA_PATH;
  string filename = fn + "tooth.nhdr";
  //char *filename = "/scratch/mjc/SRData/1.25.0/volume/tooth.nhdr";
  char *err;
  if (nrrdLoad(nin, filename.c_str(), NULL)) {
    err = biffGetDone(NRRD);
    fprintf(stderr, "%s: trouble reading \"%s\":\n%s", 
	    "test wxOGL volume canvas", filename.c_str(), err);
    free(err);
    return;
  }

  if (nin->dim != 3) {
    cerr << "Error must input a Nrrd of dim == 3" << endl;
    exit(0);
  }
  
  Nrrd *nv = 0;
  Nrrd *gm = 0;
  build_gm_nrrds(nin, nv, gm);
  

  //Quantize the gradient magnitude nrrd.
  NrrdRange *range = nrrdRangeNewSet(gm, nrrdBlind8BitRangeState);
  cerr << "gm min:max:" << range->min << ":" << range->max << endl;
  Nrrd *gm8 = nrrdNew();
  if (nrrdQuantize(gm8, gm, range, 8))
  {
    char *err = biffGetDone(NRRD);
    cerr << "Trouble quantizing: " << err << endl;
    free(err);
    return;
  }


  BBox bounds;
  Point max(nin->axis[0].size, nin->axis[1].size, nin->axis[2].size);
  bounds.extend(Point(0.0, 0.0,0.0));
  bounds.extend(max);
  //cerr << bounds << endl;
  input_bounds_ = bounds;
  bounds_set_ = true;
  //cerr << nin->axis[0].max << endl;


  tex_ = new Texture();
  //cerr << "video_card_memory_size():" << video_card_memory_size() << endl;
  tex_->build(nv, gm8, 0, 256, 0, 0, video_card_memory_size());

}

//! Build 2 nrrds from the input.
//! 1) 4 dim nrrd with 3 values for the gradient, 1 for the original value.
//! 2) 1 dim scalar with the gradient magnitude 
//! quantize the output to 8 bit values.
void
OGLCanvas::build_gm_nrrds(Nrrd *nin, Nrrd *&nvout, Nrrd *&gmout)
{
  size_t nvsize[NRRD_DIM_MAX];
  size_t gmsize[NRRD_DIM_MAX];
  size_t dim = 0;
  // Create a local array of axis sizes, so we can allocate the output Nrrd
  for (dim=0; dim < nin->dim; dim++)
  {
    gmsize[dim] = nin->axis[dim].size;
    nvsize[dim+1]=nin->axis[dim].size;
  }
  nvsize[0] = 4;

  nvout = nrrdNew();
  nrrdAlloc_nva(nvout, nrrdTypeUChar, nin->dim+1, nvsize);
  
  // Set axis info for (new) axis 0
  nvout->axis[0].kind = nrrdKind4Vector;
  nvout->axis[0].center=nin->axis[0].center;
  nvout->axis[0].spacing=AIR_NAN;
  nvout->axis[0].min=AIR_NAN;
  nvout->axis[0].max=AIR_NAN;

  for (int i = 0; i < NRRD_SPACE_DIM_MAX; i++) {
    nvout->axis[0].spaceDirection[i] = AIR_NAN;
    nvout->spaceOrigin[i] = nin->spaceOrigin[i];
  }

  // Copy the other axes' info
  for (dim = 0; dim < nin->dim; dim++)
  {
    nvout->axis[dim+1].kind = nin->axis[dim].kind;
    nvout->axis[dim+1].center = nin->axis[dim].center;
    nvout->axis[dim+1].spacing = nin->axis[dim].spacing;
    nvout->axis[dim+1].min = nin->axis[dim].min;
    nvout->axis[dim+1].max = nin->axis[dim].max;
    for (int i=0; i<NRRD_SPACE_DIM_MAX; i++) {
      nvout->axis[dim+1].spaceDirection[i] = nin->axis[dim].spaceDirection[i];
    }
  }

  // Allocate the nrrd's data, set the size of each axis
  gmout = nrrdNew();
  nrrdAlloc_nva(gmout, nrrdTypeFloat, nin->dim, gmsize);
  
  // Copy the other axes' info
  for (dim=0; dim<nin->dim; dim++)
  {
    gmout->axis[dim].kind = nin->axis[dim].kind;
    gmout->axis[dim].center = nin->axis[dim].center;
    gmout->axis[dim].spacing = nin->axis[dim].spacing;
    gmout->axis[dim].min = nin->axis[dim].min;
    gmout->axis[dim].max = nin->axis[dim].max;
  }

  Transform transform;
  transform.load_identity();
  
  // Reconstruct the axis aligned transform.
//   const Point nmin(nin->axis[0].min, nin->axis[1].min, nin->axis[2].min);
//   const Point nmax(nin->axis[0].max, nin->axis[1].max, nin->axis[2].max);
//   transform.pre_scale(nmax - nmin);
//   transform.pre_translate(nmin.asVector());

//   // Add the matrix size into the canonical transform.
//   transform.post_scale(Vector(1.0 / (nin->axis[0].size - 1.0),
//                               1.0 / (nin->axis[1].size - 1.0),
//                               1.0 / (nin->axis[2].size - 1.0)));

  NrrdRange *range = nrrdRangeNewSet(nin, nrrdBlind8BitRangeState);
  float *gmoutdata = (float *)(gmout->data);
  unsigned char *nvoutdata = (unsigned char *)nvout->data;
  unsigned short *nindata = (unsigned short *)nin->data;

  compute_data(nindata, nvoutdata, gmoutdata,
	       nin->axis[0].size, nin->axis[1].size, nin->axis[2].size,
	       transform, range->min, range->max);
}

static inline
unsigned char
NtoUC(double a)
{
  return (unsigned char)(a * 127.5 + 127.5);
}


static inline
unsigned char
VtoUC(double v, double dmin, double dmaxplus)
{
  v = (v - dmin) * dmaxplus;
  if (v > 255.0) return 255;
  if (v < 0.0) return 0;
  return (unsigned char)(v);
}


void
OGLCanvas::compute_data(unsigned short *nindata, 
			unsigned char *nvoutdata, float *gmoutdata,
			const int ni, const int nj, const int nk, 
			Transform &transform, double dmin, double dmax)
{
  int i, j, k;
  const unsigned int nji = nj * ni;
  const double dmaxplus = 255.0 / (dmax - dmin) + 1e-9;  

  for (k = 0; k < nk; k++)
  {
    double zscale = 0.5;
    int k0 = k-1;
    int k1 = k+1;
    if (k == 0)   { k0 = k; zscale = 1.0; }
    if (k1 == nk) { k1 = k; zscale = 1.0; }
    for (j = 0; j < nj; j++)
    {
      double yscale = 0.5;
      int j0 = j-1;
      int j1 = j+1;
      if (j == 0)   { j0 = j; yscale = 1.0; }
      if (j1 == nj) { j1 = j; yscale = 1.0; }
      for (i = 0; i < ni; i++)
      {
	double xscale = 0.5;
	int i0 = i-1;
	int i1 = i+1;
	if (i == 0)   { i0 = i; xscale = 1.0; }
	if (i1 == ni) { i1 = i; xscale = 1.0; }

	const unsigned int idx = k * nji + j * ni + i;
	const double val = (double)nindata[idx];
	const double x0 = nindata[k * nji + j * ni + i0];
	const double x1 = nindata[k * nji + j * ni + i1];
	const double y0 = nindata[k * nji + j0 * ni + i];
	const double y1 = nindata[k * nji + j1 * ni + i];
	const double z0 = nindata[k0 * nji + j * ni + i];
	const double z1 = nindata[k1 * nji + j * ni + i];

	const Vector g((x1 - x0)*xscale, (y1 - y0)*yscale, (z1 - z0)*zscale);
	Vector gradient = transform.project_normal(g);
	float gm = gradient.safe_normalize();
	if (gm < 1.0e-5) gm = 0.0;
	
	nvoutdata[idx * 4 + 0] = NtoUC(gradient.x());
	nvoutdata[idx * 4 + 1] = NtoUC(gradient.y());
	nvoutdata[idx * 4 + 2] = NtoUC(gradient.z());
	nvoutdata[idx * 4 + 3] = VtoUC(val, dmin, dmaxplus);
	gmoutdata[k * nji + j * ni + i] = gm;
      }
    }
  }
}

}

