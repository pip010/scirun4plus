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
//    File   : OGLCanvas.h
//    Author : Martin Cole
//    Date   : Wed Jul 11 13:46:55 2007
//    Modified from James Bigler's wxtest code.

#if !defined(OGLCanvas_h)
#define OGLCanvas_h

#include <slivr/VolumeRenderer.h>
#include <slivr/BBox.h>
#include <slivr/Point.h>

#include <wx/wx.h>
#include <wx/glcanvas.h>

namespace SLIVR_EXAMPLES {

using namespace SLIVR;

class OGLCanvas: public wxGLCanvas
{
public:
  OGLCanvas(wxWindow* parent, int* attribs);
  virtual ~OGLCanvas() {}

  // Called when there is need to draw.
  void draw(wxPaintEvent& event);


  void set_blend_num_bits_32() {
    vr_->set_blend_num_bits(32);
  }
  void set_blend_num_bits_8() {
    vr_->set_blend_num_bits(8);
  }
  void set_shading(bool b) {
    vr_->set_shading(b);
  }
  void set_blend_mode_over(bool b) {
    if (b) {
      vr_->set_mode(VolumeRenderer::MODE_OVER);
    } else {
      vr_->set_mode(VolumeRenderer::MODE_MIP);
    }
  }

  void set_use_2d_cmap(bool b) { use_2d_cmap_ = b; }

  void set_sampling_rate(float r) {
    vr_->set_sampling_rate(r);
  }

  void set_spin_p(bool b) { spin_p_ = b; }
  bool spin_p() const { return spin_p_; }
  
private:
  void init();
  void build_texture();
  bool viewall();
  ColorMap2* create_cmap2();
  bool compute_depth(const SLIVR::Point &veyep, 
		     const SLIVR::Point &vlookat, 
		     double& znear, double& zfar);

  void build_gm_nrrds(Nrrd *nin, Nrrd *&nv, Nrrd *&gm);
  void compute_data(unsigned short *nindata, 
		    unsigned char *nvoutdata, float *gmoutdata,
		    const int ni, const int nj, const int nk, 
		    Transform &transform, double dmin, double dmax);

  bool                              init_done_;
  bool                              bounds_set_;
  bool                              use_2d_cmap_;
  bool                              last_use_2d_cmap_;
  BBox                              input_bounds_;
  ColorMap                         *cm_;
  vector<ColorMap2*>                cmap2_;
  vector<Plane*>                    planes_;

  Texture                          *tex_;
  VolumeRenderer                   *vr_;
  
  bool                              spin_p_;
  wxWindow                         *window_;
  
  DECLARE_EVENT_TABLE();
};

}
#endif // #ifndef OGLCanvas_h
