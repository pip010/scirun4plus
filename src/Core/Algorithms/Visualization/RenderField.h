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
 *  RenderFieldData.h: Rendering alogrithms for data
 *
 *  Written by:
 *   Allen R. Sanderson
 *   Martin Cole
 *   SCI Institute
 *   University of Utah
 *   April 2007
 *
 */

#if !defined(Visualization_RenderField_h)
#define Visualization_RenderField_h

#include <Core/Algorithms/Visualization/RenderFieldState.h>
#include <Core/Basis/QuadBilinearLgn.h>
#include <Core/Util/StringUtil.h>
#include <Core/Datatypes/Field.h>
#include <Core/Geom/GeomMaterial.h>
#include <Core/Geom/GeomGroup.h>
#include <Core/Geom/GeomSwitch.h>
#include <Core/Geom/GeomGlyph.h>
#include <Core/Geom/GeomLine.h>
#include <Core/Geom/GeomPoint.h>
#include <Core/Geom/GeomTriangles.h>
#include <Core/Geom/GeomQuads.h>
#include <Core/Geom/GeomText.h>
#include <Core/Geom/GeomDL.h>
#include <Core/Geom/GeomColorMap.h>
#include <Core/Geom/GeomTexRectangle.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Tensor.h>
#include <Core/Math/MiscMath.h>
#include <Core/Math/TrigTable.h>
#include <Core/Util/TypeDescription.h>
#include <sci_hash_map.h>
#include <sstream>
#include <iomanip>
#include <iostream>

#if defined(_WIN32) && !defined(uint)
  // for some reason, this isn't defined...
#define uint unsigned
#endif


#include <Core/Algorithms/Visualization/share.h>

namespace SCIRun {

//! RenderFieldBase supports the dynamically loadable algorithm concept.
//! when dynamically loaded the user will dynamically cast to a 
//! RenderFieldBase from the DynamicAlgoBase they will have a pointer to.
class SCISHARE RenderFieldBase :
  public RenderStateBase, public UsedWithLockingHandle<Mutex>
{
  public:

#if !defined(_MSC_VER) && !defined(__ECC)

  struct equal_str
  {
    bool operator()(const std::string &s1, const std::string &s2) const
    {
      return s1 == s2;
    }
  };
  
  struct str_hasher
  {
    size_t operator()(const std::string& s) const
    {
      hash<const char*> H;
      return H(s.c_str());
    }
  };
#endif

  virtual GeomHandle render_mesh_nodes(FieldHandle field_handle, 
				       ColorMapHandle colormap_handle,
				       const std::string &display_mode,
				       double scale,
				       int resolution,
				       bool use_def_color,
				       bool use_transparency)
  {
    unsigned int render_state = 0;
    switch_flag( render_state, IS_ON, 1 );
    switch_flag( render_state, USE_DEFAULT_COLOR, use_def_color );
    switch_flag( render_state, USE_TRANSPARENCY,  use_transparency );
    switch_flag( render_state, USE_COLORMAP,      colormap_handle.get_rep() );

    if (field_handle->basis_order() == -1 || colormap_handle.get_rep() == 0)
      switch_flag( render_state, USE_DEFAULT_COLOR, true );
    
    return render_mesh_nodes( field_handle,
			      display_mode, scale, resolution, render_state);
  }
  
  virtual GeomHandle render_mesh_edges(FieldHandle field_handle,
				       ColorMapHandle colormap_handle,
				       const std::string &display_mode,
				       double scale,
				       int resolution,
				       bool use_def_color,
				       bool use_transparency,
				       unsigned int approx_div)
  {
    unsigned int render_state = 0;
    switch_flag( render_state, IS_ON, 1 );
    switch_flag( render_state, USE_DEFAULT_COLOR, use_def_color );
    switch_flag( render_state, USE_TRANSPARENCY,  use_transparency );
    switch_flag( render_state, USE_COLORMAP,      colormap_handle.get_rep() );

    if (field_handle->basis_order() == -1 || colormap_handle.get_rep() == 0)
      switch_flag( render_state, USE_DEFAULT_COLOR, true );

    return render_mesh_edges( field_handle, display_mode,
			      scale, resolution, render_state, approx_div);

  }

  virtual GeomHandle render_mesh_faces(FieldHandle field_handle, 
				       ColorMapHandle colormap_handle,
				       bool use_normals,
				       bool use_def_color,
				       bool use_transparency,
				       bool use_texture,
				       unsigned int approx_div)
  {
    unsigned int render_state = 0;
    switch_flag( render_state, IS_ON, 1 );
    switch_flag( render_state, USE_DEFAULT_COLOR, use_def_color );
    switch_flag( render_state, USE_TRANSPARENCY,  use_transparency );
    switch_flag( render_state, USE_COLORMAP,      colormap_handle.get_rep() );
    switch_flag( render_state, USE_NORMALS, use_normals );
    switch_flag( render_state, USE_TEXTURE, use_texture );

    if (field_handle->basis_order() == -1 || colormap_handle.get_rep() == 0)
      switch_flag( render_state, USE_DEFAULT_COLOR, true );

    return render_mesh_faces(field_handle, colormap_handle,
			     render_state, approx_div);
  }


  virtual GeomHandle render_mesh_nodes(FieldHandle field_handle, 
				       const std::string &display_mode,
				       double scale, int resolution,
				       unsigned int render_state) = 0;

  virtual GeomHandle render_mesh_edges(FieldHandle field_handle,
				       const std::string &display_mode,
				       double scale, int resolution,
				       unsigned int render_state,
				       unsigned int approx_div) = 0;

  virtual GeomHandle render_mesh_faces(FieldHandle field_handle, 
				       ColorMapHandle color_handle,
				       unsigned int render_state,
				       unsigned int approx_div) = 0;

  virtual GeomHandle render_text(FieldHandle field_handle,
				 ColorMapHandle color_handle,
				 bool use_colormap,
				 bool use_default_color,
				 bool backface_cull_p,
				 int  fontsize,
				 int  precision,
				 bool render_locations,
				 bool render_data,
				 bool render_nodes,
				 bool render_edges,
				 bool render_faces,
				 bool render_cells,
				 bool always_visible) = 0;

  RenderFieldBase();
  virtual ~RenderFieldBase();
  
  GeomHandle node_switch_;
  GeomHandle edge_switch_;
  GeomHandle face_switch_;
};


class SCISHARE RenderFieldV : public RenderFieldBase
{
public:
  //! virtual interface. 
  virtual GeomHandle render_mesh_nodes(FieldHandle field_handle, 
				       const std::string &display_mode,
				       double scale,
				       int resolution,
				       unsigned int render_state );

  virtual GeomHandle render_mesh_edges(FieldHandle field_handle,
				       const std::string &display_mode,
				       double scale,
				       int resolution,
				       unsigned int render_state,
				       unsigned int approx_div);

  virtual GeomHandle render_mesh_faces(FieldHandle field_handle, 
				       ColorMapHandle color_handle,
				       unsigned int render_state,
				       unsigned int approx_div);

  virtual GeomHandle render_text(FieldHandle field_handle,
				 ColorMapHandle color_handle,
				 bool use_color_map,
				 bool use_default_color,
				 bool backface_cull_p,
				 int fontsize,
				 int precision,
				 bool render_locations,
				 bool render_data,
				 bool render_nodes,
				 bool render_edges,
				 bool render_faces,
				 bool render_cells,
				 bool always_visible);

protected:
  GeomHandle render_nodes(FieldHandle field_handle,
			  const std::string &display_mode,
			  double scale,
			  int resolution,
			  unsigned int render_state);


  GeomHandle render_edges(FieldHandle field_handle,
			  const std::string &display_mode,
			  double scale,
			  int resolution,
			  unsigned int render_state,
			  unsigned int approx_div);

  GeomHandle render_edges_linear(FieldHandle field_handle,
				 const std::string &display_mode,
				 double scale,
				 int resolution,
				 unsigned int render_state);

  GeomHandle render_faces(FieldHandle field_handle,
			  unsigned int render_state,
			  unsigned int approx_div);

  GeomHandle render_faces_linear(FieldHandle field_handle,
				 unsigned int render_state);         

  GeomHandle render_text_data(FieldHandle field_handle,
			      ColorMapHandle color_handle,
			      bool use_color_map,
			      bool use_default_color,
			      bool backface_cull_p,
			      int fontsize,
			      int precision,
			      bool always_visible);

  GeomHandle render_text_data_nodes(FieldHandle field_handle,
				    ColorMapHandle color_handle,
				    bool use_color_map,
				    bool use_default_color,
				    bool backface_cull_p,
				    int fontsize,
				    int precision,
				    bool always_visible);

  GeomHandle render_text_nodes(FieldHandle field_handle,
			       ColorMapHandle color_handle,
			       bool use_color_map,
			       bool use_default_color,
			       bool backface_cull_p,
			       int fontsize,
			       int precision,
			       bool render_locations,
			       bool always_visible);

  GeomHandle render_text_edges(FieldHandle field_handle,
			       ColorMapHandle color_handle,
			       bool use_color_map,
			       bool use_default_color,
			       int fontsize,
			       int precision,
			       bool render_locations,
			       bool always_visible);

  GeomHandle render_text_faces(FieldHandle field_handle,
			       ColorMapHandle color_handle,
			       bool use_color_map,
			       bool use_default_color,
			       int fontsize,
			       int precision,
			       bool render_locations,
			       bool always_visible);

  GeomHandle render_text_cells(FieldHandle field_handle,
			       ColorMapHandle color_handle,
			       bool use_color_map,
			       bool use_default_color,
			       int fontsize,
			       int precision,
			       bool render_locations,
			       bool always_visible);

  GeomHandle render_faces_texture(FieldHandle field_handle, 
					  ColorMapHandle colormap_handle,
					  unsigned int render_state);

  void add_edge_geom(GeomLines *lines,
		     GeomGlyphBase *glyphs,
		     const Point &p0, const Point &p1,
		     const double radius,
		     const int resolution,
		     bool cylinders_p,
		     unsigned int color_scheme,
		     double scol0,
		     double scol1,
		     MaterialHandle& vcol0,
		     MaterialHandle& vcol1);

  void add_face_geom(GeomFastTriangles *faces,
		     GeomFastQuads *qfaces,
		     const std::vector<Point>  &points,
		     const std::vector<Vector> &normals,
		     bool with_normals,
		     unsigned int color_scheme,
		     std::vector<double> &scols,
		     std::vector<MaterialHandle> &vcols );
};


//#define DEBUG_PRINT
#if defined(DEBUG_PRINT)
inline
void
print_coords(const char *pre, std::vector<double> &c0)
{
  cout << pre;
  std::vector<double>::iterator dbg_iter = c0.begin();
  while (dbg_iter != c0.end()) {
    cout << *dbg_iter++ << " ";
  }
  cout << endl;
}
#endif


} // end namespace SCIRun

#endif // Visualization_RenderField_h
