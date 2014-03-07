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
 *  RenderFieldGlyphs.cc: Rendering alogrithms for data
 *
 *  Written by:
 *   Allen R. Sanderson
 *   SCI Institute
 *   University of Utah
 *   April 2007
 *
 */

#include <Core/Algorithms/Visualization/RenderFieldGlyphs.h>

namespace SCIRun {

void
RenderScalarField::add_axis(const Point &p0, double scale,
				GeomLines *lines)
{
  static const Vector x(1., 0., 0.);
  static const Vector y(0., 1., 0.);
  static const Vector z(0., 0., 1.);

  Point p1 = p0 + x * scale;
  Point p2 = p0 - x * scale;
  lines->add(p1, p2);

  p1 = p0 + y * scale;
  p2 = p0 - y * scale;
  lines->add(p1, p2);

  p1 = p0 + z * scale;
  p2 = p0 - z * scale;
  lines->add(p1, p2);
}


void 
RenderScalarField::add_axis(const Point &p0, double scale,
				GeomLines *lines,
				const MaterialHandle &vcol)
{
  static const Vector x(1., 0., 0.);
  static const Vector y(0., 1., 0.);
  static const Vector z(0., 0., 1.);

  Point p1 = p0 + x * scale;
  Point p2 = p0 - x * scale;
  lines->add(p1, vcol, p2, vcol);

  p1 = p0 + y * scale;
  p2 = p0 - y * scale;
  lines->add(p1, vcol, p2, vcol);

  p1 = p0 + z * scale;
  p2 = p0 - z * scale;
  lines->add(p1, vcol, p2, vcol);
}


void 
RenderScalarField::add_axis(const Point &p0, double scale,
				GeomLines *lines,
				double cindex)
{
  static const Vector x(1., 0., 0.);
  static const Vector y(0., 1., 0.);
  static const Vector z(0., 0., 1.);

  Point p1 = p0 + x * scale;
  Point p2 = p0 - x * scale;
  lines->add(p1, cindex, p2, cindex);

  p1 = p0 + y * scale;
  p2 = p0 - y * scale;
  lines->add(p1, cindex, p2, cindex);

  p1 = p0 + z * scale;
  p2 = p0 - z * scale;
  lines->add(p1, cindex, p2, cindex);
}


void 
RenderScalarField::add_axis(const Point &p0, double scale,
				GeomLines *lines,
				const MaterialHandle &vcol,
				double cindex)
{
  static const Vector x(1., 0., 0.);
  static const Vector y(0., 1., 0.);
  static const Vector z(0., 0., 1.);

  Point p1 = p0 + x * scale;
  Point p2 = p0 - x * scale;
  lines->add(p1, vcol, cindex, p2, vcol, cindex);

  p1 = p0 + y * scale;
  p2 = p0 - y * scale;
  lines->add(p1, vcol, cindex, p2, vcol, cindex);

  p1 = p0 + z * scale;
  p2 = p0 - z * scale;
  lines->add(p1, vcol, cindex, p2, vcol, cindex);
}

double
RenderTensorField::map_emphasis(double old)
{
  if (old < 0.0) old = 0.0;
  else if (old > 1.0) old = 1.0;
  return tan(old * (M_PI / 2.0 * 0.999));
  // Map old 3.5 value onto new 0.825 value.
  //return tan(old * (atan(3.5) / (0.825 * 4.0)));
}


//

class SVTvalue {
  public:
    Tensor tensor;
    Vector vector;
    double scalar;

    int type;
};

inline bool get_value(VField* field,SVTvalue& val, VField::index_type idx)
{
  if (field->is_scalar())
  {
    field->get_value(val.scalar,idx);
    val.type = 1;
  }
  if (field->is_vector())
  {
    field->get_value(val.vector,idx);
    val.type = 2;
  }
  if (field->is_tensor())
  {
    field->get_value(val.tensor,idx);
    val.type = 3;
  }
  return (true);
}

inline bool to_double(SVTvalue& val,double &d)
{
  if (val.type == 1) return(to_double(val.scalar,d));
  else if (val.type == 2) return(to_double(val.vector,d));
  else return(to_double(val.tensor,d));
}

inline bool to_color(SVTvalue& val,Color &d)
{
  if (val.type == 1) return(to_color(val.scalar,d));
  else if (val.type == 2) return(to_color(val.vector,d));
  else return(to_color(val.tensor,d));
}

inline bool to_vector(SVTvalue& val,Vector &d)
{
  if (val.type == 1) return(to_vector(val.scalar,d));
  else if (val.type == 2) return(to_vector(val.vector,d));
  else return(to_vector(val.tensor,d));
}

inline bool to_tensor(SVTvalue& val,Tensor &d)
{
  if (val.type == 1) return(to_tensor(val.scalar,d));
  else if (val.type == 2) return(to_tensor(val.vector,d));
  else return(to_tensor(val.tensor,d));
}

void value_to_color(unsigned int color_scheme, SVTvalue& val, double& scol, MaterialHandle &vcol)
{
  if (val.type == 1) value_to_color(color_scheme,val.scalar,scol,vcol);
  else if (val.type == 1) value_to_color(color_scheme,val.vector,scol,vcol);
  else value_to_color(color_scheme,val.tensor,scol,vcol);
}

// emphasis was 3.5 and looked reasonable.  It should currently be
// computed via 'tan(VAL * M_PI / 2 * 0.999)', where VAL is [0,1].

GeomHandle 
RenderScalarField::render_data(FieldHandle pfld_handle,
                               FieldHandle sfld_handle,
                               FieldHandle tfld_handle,
                               const std::string &display_mode,
                               double scale,
                               double sscale,
                               double tscale,
                               double threshold,
                               int resolution,
                               unsigned int render_state,
                               unsigned int secondary_render_state,
                               unsigned int tertiary_render_state)
{
  VField *pfld = 0; if (pfld_handle.get_rep()) pfld = pfld_handle->vfield(); 
  VField *sfld = 0; if (sfld_handle.get_rep()) sfld = sfld_handle->vfield(); 
  VField *tfld = 0; if (tfld_handle.get_rep()) tfld = tfld_handle->vfield(); 
  
  bool points_p  = (display_mode == "Points");
  bool spheres_p = (display_mode == "Spheres");
  bool boxes_p   = (display_mode == "Boxes");
  bool axes_p    = (display_mode == "Axes");

  GeomPoints  *points = 0;
  GeomLines *lines = 0;
  GeomGlyphBase *glyphs = 0;
  GeomGroup *grp = new GeomGroup();
  GeomHandle data_switch = new GeomDL(grp);

  // Points when too small or when requested
  if (get_flag(render_state, USE_TRANSPARENCY) )
    points = new GeomTranspPoints();
  else
    points = new GeomPoints();
  
  grp->add(points);

  points->setPointSize(resolution/5.0);

  if (axes_p) // Axis
  {
    scale /= 2.0;

    if (get_flag(render_state, USE_TRANSPARENCY))
      lines = new GeomTranspLines();
    else
      lines = new GeomLines();

    grp->add(lines);

    lines->setLineWidth(resolution/5.0);
  }

  if (spheres_p || boxes_p) // Spheres or Boxes
  {
    if (spheres_p)
      scale /= 2.0;

    if (get_flag(render_state, USE_TRANSPARENCY))
      glyphs = new GeomTranspGlyph();
    else
      glyphs = new GeomGlyph();

    grp->add(glyphs->getObj());
  }

  Vector box_vec(0,0,1); // Needed for the box orienation

  bool normalize = get_flag(render_state, NORMALIZE_DATA);
  bool small_is_dot = get_flag(render_state,SMALL_IS_DOT);
  bool secondary_on    = get_flag(secondary_render_state, IS_ON);
  bool secondary_color = (get_flag(secondary_render_state, IS_ON|USE_COLORMAP) ||
			  get_flag(secondary_render_state, IS_ON|USE_COLOR_CONVERT));
  bool secondary_alpha = get_flag(secondary_render_state, IS_ON|USE_ALPHA) &&
    get_flag(render_state, USE_TRANSPARENCY);
  bool secondary_value = get_flag(secondary_render_state, IS_ON|USE_VALUE);
  //bool secondary_small_is_dot = get_flag(secondary_render_state,SMALL_IS_DOT);


  bool tertiary_on    = get_flag(tertiary_render_state, IS_ON);
  bool tertiary_color = (get_flag(tertiary_render_state, IS_ON|USE_COLORMAP) ||
			 get_flag(tertiary_render_state, IS_ON|USE_COLOR_CONVERT));
  bool tertiary_alpha = get_flag(tertiary_render_state, IS_ON|USE_ALPHA) &&
    get_flag(render_state, USE_TRANSPARENCY);
  bool tertiary_value = get_flag(tertiary_render_state, IS_ON|USE_VALUE);
  //bool tertiary_small_is_dot = get_flag(tertiary_render_state,SMALL_IS_DOT);


  unsigned int color_scheme = 0;
  double scol;
  MaterialHandle vcol(0);

  if( !tertiary_color && !secondary_color &&
      get_flag(render_state, USE_DEFAULT_COLOR) )
  {
    color_scheme = 0; // Default color
  }
  else if ( (!tertiary_color && !secondary_color &&
	     get_flag(render_state, USE_COLORMAP)) ||
	    ( secondary_color && get_flag(secondary_render_state, USE_COLORMAP)) ||
	    ( tertiary_color && get_flag(tertiary_render_state, USE_COLORMAP)) )
  {
    // Color map lookup using either a scalar value
    // or the vector magnitude.
    color_scheme = 1;
  }
  else if( pfld->basis_order() >= 0 ||
           (secondary_color && sfld->basis_order() >= 0) ||
           (tertiary_color  && tfld->basis_order() >= 0) )
  {
    color_scheme = 2; // Values become RGB
  }


  if( color_scheme == 2 || secondary_alpha || tertiary_alpha )
  {
    vcol = new Material(Color(1.0, 1.0, 1.0));
    if (get_flag(render_state, USE_TRANSPARENCY) )
      vcol->transparency = 0.75;
    else
      vcol->transparency = 1.0;
  }


  VMesh* mesh = pfld->vmesh();
  
  VMesh::size_type num_values = pfld->num_values();
   
  for(VMesh::index_type idx=0; idx<num_values; idx++)
  {
    Point p;
    if (pfld->basis_order() == 0)
      mesh->get_center(p,VMesh::Elem::index_type(idx));
    else
      mesh->get_center(p,VMesh::Node::index_type(idx));    
    
    double val;
    SVTvalue sval;
    SVTvalue tval;

    pfld->get_value(val, idx);
    {
      if( secondary_on )
        get_value(sfld, sval, idx);
      
      if( tertiary_on )
        get_value(tfld, tval, idx);
      
      if (color_scheme)
      {
        if( tertiary_color )
          value_to_color( color_scheme, tval, scol, vcol );
        else if( secondary_color )
          value_to_color( color_scheme, sval, scol, vcol );
        else
          value_to_color( color_scheme,  val, scol, vcol );
      }

      if( tertiary_alpha )
      {
        double atmp;
        to_double(tval, atmp);
        vcol->transparency = fabs(atmp);
      }
      else if( secondary_alpha )
      {
        double atmp;
        to_double(sval, atmp);
        vcol->transparency = fabs(atmp);
      }

      // The color is done so the scalar value can be changed
      if (normalize)
        val = 1.0;

      if (fabs(val) > threshold && !points_p)
      {
        val *= scale;

        double value = fabs(val);

        if (secondary_value)
        {
          to_vector(sval, box_vec);
        }

        double svalue;
        double tvalue;

        if (secondary_value)
        {
          double tmp;
          to_double(sval, tmp);
          svalue = fabs(tmp * sscale);
        }
        else
        {
          svalue = value;
        }

        if (tertiary_value)
        {
          double tmp;
          to_double(tval, tmp);
          tvalue = fabs(tmp * tscale);
        }
        else
        {
          tvalue = value;
        }

        if( spheres_p )
        {
          if (color_scheme == 0)
            glyphs->add_sphere(p, value, resolution, resolution);
          else if (color_scheme == 2)
            glyphs->add_sphere(p, value, vcol, resolution, resolution);
          else // if (color_scheme == 1)
          {
            if( secondary_alpha )
              glyphs->add_sphere(p, value, vcol, scol, resolution, resolution);
            else
              glyphs->add_sphere(p, value, scol, resolution, resolution);
          }
        }
        else if( boxes_p )
        {
          if (color_scheme == 0)
            glyphs->add_box(p, box_vec, value, svalue, tvalue);
          else if (color_scheme == 2)
            glyphs->add_box(p, box_vec, value, svalue, tvalue, vcol);
          else //if (color_scheme == 1)
          {
            if( secondary_alpha )
              glyphs->add_box(p, box_vec, value, svalue, tvalue, vcol, scol);
            else
              glyphs->add_box(p, box_vec, value, svalue, tvalue, scol);
          }
        }
        else if (axes_p)
        {
          if (color_scheme == 0)
            add_axis(p, value, lines);
          else if (color_scheme == 2)
            add_axis(p, value, lines, vcol);
          else //if (color_scheme == 1)
          {
            if( secondary_alpha )
              add_axis(p, value, lines, vcol, scol);
            else
              add_axis(p, value, lines, scol);
          }
        }
      }
      else // if (points_p)
      {
        if (small_is_dot)
        {
          if (color_scheme == 0)
            points->add(p);
          else if (color_scheme == 2)
            points->add(p, vcol);
          else //if (color_scheme == 1)
          {
            if( secondary_alpha )
              points->add(p, vcol, scol);
            else
              points->add(p, scol);
          }
        }
      }
    }
  }
  return data_switch;
}


GeomHandle 
RenderVectorField::render_data(FieldHandle pfld_handle,
                               FieldHandle sfld_handle,
                               FieldHandle tfld_handle,
                               MaterialHandle default_material,
                               const std::string &display_mode,
                               double scale,
                               double sscale,
                               double tscale,
                               double threshold,
                               int resolution,
                               unsigned int render_state,
                               unsigned int secondary_render_state,
                               unsigned int tertiary_render_state)
{
  VField* pfld = 0;
  VField* sfld = 0;
  VField* tfld = 0;
  if (pfld_handle.get_rep()) pfld = pfld_handle->vfield();
  if (sfld_handle.get_rep()) sfld = sfld_handle->vfield();
  if (tfld_handle.get_rep()) tfld = tfld_handle->vfield();

  const bool arrows_p  = (display_mode == "Arrows");
  const bool cones_p   = (display_mode == "Cones");
  const bool comets_p  = (display_mode == "Comets");
  const bool springs_p = (display_mode == "Springs");
  const bool disks_p   = (display_mode == "Disks");
  const bool rings_p   = (display_mode == "Rings");
  const bool needles_p = (display_mode == "Needles");
  const bool lines_p   = (display_mode == "Lines") ||
     arrows_p || cones_p;

  // Note: arrows, cones, and springs default to lines if the
  // radius is too small.

  GeomPoints  *points = 0;
  GeomLines *lines = 0;
  GeomGlyphBase *glyphs = 0;
  GeomGroup *grp = new GeomGroup();
  GeomHandle data_switch = new GeomDL(grp);

  // Points when too small
  if (get_flag(render_state, USE_TRANSPARENCY) )
    points = new GeomTranspPoints();
  else
    points = new GeomPoints();

  grp->add(points);

  points->setPointSize(resolution/5.0);

  // Note: comets and rings default to linear sections if the radius
  // is too small.
  if (lines_p || needles_p || comets_p || rings_p || springs_p)
  {
    if( needles_p || comets_p || get_flag(render_state, USE_TRANSPARENCY))
      lines = new GeomTranspLines();
    else
      lines = new GeomLines();

    grp->add(lines);

    lines->setLineWidth(resolution/5.0);
  }

  if (arrows_p || cones_p || disks_p || comets_p || rings_p || springs_p)
  {
    if( get_flag(render_state, USE_TRANSPARENCY) )
      glyphs = new GeomTranspGlyph();
    else 
      glyphs = new GeomGlyph();

    grp->add(glyphs->getObj());
  }

  // Reduce the scale by 1/2 because these are radius based glyphs
  // or in the case of the arrows 1/2 is a cone and 1/2 is a line.
  if (arrows_p || disks_p || rings_p)
    scale *= (1.0 / 2.0);

  bool normalize       = get_flag(render_state, NORMALIZE_DATA);
  bool bidirectional   = get_flag(render_state, BIDIRECTIONAL);
  bool small_is_dot    = get_flag(render_state,SMALL_IS_DOT);

  // Secondary
  bool secondary_on    = get_flag(secondary_render_state, IS_ON);
  bool secondary_color = (get_flag(secondary_render_state, IS_ON|USE_COLORMAP) ||
			  get_flag(secondary_render_state, IS_ON|USE_COLOR_CONVERT));
  bool secondary_alpha = get_flag(secondary_render_state, IS_ON|USE_ALPHA) &&
    get_flag(render_state, USE_TRANSPARENCY);
  bool secondary_value = get_flag(secondary_render_state, IS_ON|USE_VALUE);
  bool secondary_major = get_flag(secondary_render_state, IS_ON|USE_VALUE|USE_MAJOR_RADIUS);
  bool secondary_minor = get_flag(secondary_render_state, IS_ON|USE_VALUE|USE_MINOR_RADIUS);
  bool secondary_pitch = get_flag(secondary_render_state, IS_ON|USE_VALUE|USE_PITCH);
  //  bool secondary_small_is_dot = get_flag(secondary_render_state,SMALL_IS_DOT);

  // Tertiary
  bool tertiary_on    = get_flag(tertiary_render_state, IS_ON);
  bool tertiary_color = (get_flag(tertiary_render_state, IS_ON|USE_COLORMAP) ||
			 get_flag(tertiary_render_state, IS_ON|USE_COLOR_CONVERT));
  bool tertiary_alpha = get_flag(tertiary_render_state, IS_ON|USE_ALPHA) &&
    get_flag(render_state, USE_TRANSPARENCY);
  bool tertiary_value = get_flag(tertiary_render_state, IS_ON|USE_VALUE);
  bool tertiary_major = get_flag(tertiary_render_state, IS_ON|USE_VALUE|USE_MAJOR_RADIUS);
  bool tertiary_minor = get_flag(tertiary_render_state, IS_ON|USE_VALUE|USE_MINOR_RADIUS);
  bool tertiary_pitch = get_flag(tertiary_render_state, IS_ON|USE_VALUE|USE_PITCH);
  //  bool tertiary_small_is_dot = get_flag(tertiary_render_state,SMALL_IS_DOT);

  double ratio = 6.0;

  unsigned int color_scheme = 0;
  double scol;
  MaterialHandle vcol(0);

  if( !tertiary_color && !secondary_color &&
      get_flag(render_state, USE_DEFAULT_COLOR) )
  {
    color_scheme = 0; // Default color
  }
  else if ( (!tertiary_color && !secondary_color &&
	     get_flag(render_state, USE_COLORMAP)) ||
	    ( secondary_color && get_flag(secondary_render_state, USE_COLORMAP)) ||
	    ( tertiary_color && get_flag(tertiary_render_state, USE_COLORMAP)) )
  {
    // Color map lookup using either a scalar value
    // or the vector magnitude.
    color_scheme = 1;
  }
  else if( pfld->basis_order() >= 0 ||
           (secondary_color && sfld->basis_order() >= 0) ||
           (tertiary_color  && tfld->basis_order() >= 0) )
  {
    color_scheme = 2; // Values become RGB
  }

  if( color_scheme == 2 || secondary_alpha || tertiary_alpha )
  {
    vcol = new Material(Color(1.0, 1.0, 1.0));
    if (get_flag(render_state, USE_TRANSPARENCY) )
      vcol->transparency = 0.75;
    else
      vcol->transparency = 1.0;
  }

  MaterialHandle opaque, transparent;
  if (needles_p || comets_p)
  {
    if (color_scheme == 0)
    {
      opaque = new Material(default_material->diffuse);
      opaque->transparency = default_material->transparency;
      transparent = new Material(default_material->diffuse);
    }
    else
    {
      opaque = new Material(Color(1.0, 1.0, 1.0));
      opaque->transparency = 1.0;
      transparent = new Material(Color(1.0, 1.0, 1.0));
    }

    if (needles_p)
      transparent->transparency = 0.20;
    else // if(comets_p)
      transparent->transparency = 0.33;
  }

  VMesh* mesh = pfld->vmesh();

  VMesh::size_type sz = pfld->num_values();
  for (VMesh::index_type idx=0; idx<sz; idx++)
  {
    Point p;
    if (pfld->basis_order() == 0)
      mesh->get_center(p, VMesh::Elem::index_type(idx));
    else
      mesh->get_center(p, VMesh::Node::index_type(idx));
    
    Vector val;
    SVTvalue sval;
    SVTvalue tval;

    pfld->get_value(val, idx);
    {
      if( secondary_on )
        get_value(sfld, sval, idx);
      
      if( tertiary_on )
        get_value(tfld, tval, idx);
      
      if (color_scheme)
      {
        if( tertiary_color )
          value_to_color( color_scheme, tval, scol, vcol );
        else if( secondary_color )
          value_to_color( color_scheme, sval, scol, vcol );
        else
          value_to_color( color_scheme,  val, scol, vcol );
      }

      if( tertiary_alpha )
      {
        double atmp;
        to_double(tval, atmp);
        vcol->transparency = fabs(atmp);
      }
      else if( secondary_alpha )
      {
        double tmp;
        to_double(sval, tmp);
        vcol->transparency = fabs(tmp);
      }

      // The color is done so the vector value can be changed
      if (normalize)
        val.safe_normalize();

      if (val.length() > threshold)
      {
        val *= scale;

        double value = val.length();
        double svalue;
        double tvalue;

        if (secondary_value)
        {
          double tmp;
          to_double(sval, tmp);
          svalue = fabs(tmp * sscale);
        }
        else
        {
          svalue = fabs(value / ratio);
        }

        if (tertiary_value)
        {
          double tmp;
          to_double(tval, tmp);
          tvalue = fabs(tmp * tscale);
        }
        else
        {
          tvalue = fabs(value / ratio / ratio);
        }

        if (arrows_p && svalue > 1e-5)
        {
          if (color_scheme == 0)
          {
            glyphs->add_cylinder(p+val, val, svalue, 0.0, value,
               resolution);
            
            if (bidirectional)
            {
              glyphs->add_cylinder(p-val, -val, svalue, 0.0, value,
                 resolution);
              lines->add(p-val, p+val);
            }
            else 
            {
              lines->add(p, p+val);
            }
          }
          else if (color_scheme == 2)
          {
            glyphs->add_cylinder(p+val, val, svalue, 0.0, value, vcol,
               resolution);
              
            if (bidirectional)
            {
              glyphs->add_cylinder(p-val, -val, svalue, 0.0, value, vcol,
                 resolution);
              lines->add(p-val, vcol, p+val, vcol);
            }
            else 
            {
              lines->add(p, vcol, p+val, vcol);
            }
          }
          else //if (color_scheme == 1)
          {
            if( secondary_alpha )
            {
              glyphs->add_cylinder(p+val, val, svalue, 0.0, value, vcol, scol, 
                 resolution);

              if (bidirectional)
              {
                glyphs->add_cylinder(p-val, -val, svalue, 0.0, value, vcol, scol, 
                   resolution);
                lines->add(p-val, vcol, scol, p+val, vcol, scol);
              }
              else 
              {
                lines->add(p, vcol, scol, p+val, vcol, scol);
              }
            }
            else
            {
              glyphs->add_cylinder(p+val, val, svalue, 0.0, value, scol, 
                 resolution);

              if (bidirectional)
              {
                glyphs->add_cylinder(p-val, -val, svalue, 0.0, value, scol, 
                   resolution);
                lines->add(p-val, scol, p+val, scol);
              }
              else 
              {
                lines->add(p, scol, p+val, scol);
              }
            }

          }
        }
        else if (cones_p && svalue> 1e-5)
        {
          if (color_scheme == 0)
            glyphs->add_cylinder(p, val, svalue, 0.0, value,
               resolution);
          else if (color_scheme == 2)
            glyphs->add_cylinder(p, val, svalue, 0.0, value, vcol,
               resolution);
          else //if (color_scheme == 1)
          {
            if( secondary_alpha )
              glyphs->add_cylinder(p, val, svalue, 0.0, value, vcol, scol, 
                 resolution);
            else
              glyphs->add_cylinder(p, val, svalue, 0.0, value, scol, 
                 resolution);
          }

          if (bidirectional)
          {
            if (color_scheme == 0)
              glyphs->add_cylinder(p-val, val, 0.0, svalue, value,
                 resolution);
            else if (color_scheme == 2)
              glyphs->add_cylinder(p-val, val, 0.0, svalue, value, vcol,
                 resolution);
            else //if (color_scheme == 1)
            {
              if( secondary_alpha )
                glyphs->add_cylinder(p-val, val, 0.0, svalue, value, vcol, scol, 
                   resolution);
              else
                glyphs->add_cylinder(p-val, val, 0.0, svalue, value, scol, 
                   resolution);
            }
          }
        }
        else if (comets_p && svalue > 1e-5)
        {
          if (color_scheme == 0)
          {
            glyphs->add_ellipsoid(p, val, svalue, opaque,
                resolution, resolution, 1);
            glyphs->add_cylinder(p-val, val, 0.0, svalue, value,
               transparent, opaque, resolution);
          }

          else if (color_scheme == 2)
          {
            transparent->diffuse = vcol->diffuse;

            glyphs->add_ellipsoid(p, val, svalue, vcol,
                resolution, resolution, 1);
            glyphs->add_cylinder(p-val, val, 0.0, svalue, value,
               transparent, vcol,
               resolution);
          }

          else //if (color_scheme == 1)
          {
            if( secondary_alpha )
            {
              glyphs->add_ellipsoid(p, val, svalue, vcol, scol,
                  resolution, resolution, 1);
              glyphs->add_cylinder(p-val, val, 0.0, svalue, value,
                 transparent, scol, vcol, scol,
                 resolution);
            }
            else
            {
              glyphs->add_ellipsoid(p, val, svalue, opaque, scol,
                  resolution, resolution, 1);
              glyphs->add_cylinder(p-val, val, 0.0, svalue, value,
                 transparent, scol, opaque, scol,
                 resolution);
            }
          }
        }
        else if (comets_p)
        {
          transparent->transparency = 0.20;

          if (color_scheme == 0)
          {
            lines->add(p, opaque, p - val, transparent);
          }
          else if (color_scheme == 2)
          {
            transparent->diffuse = vcol->diffuse;
            
            lines->add(p, vcol, p - val, transparent);
          }
          else //if (color_scheme == 1)
          {
            if( secondary_alpha )
              lines->add(p, vcol, scol, p - val, transparent, scol);
            else
              lines->add(p, opaque, scol, p - val, transparent, scol);
          }
          
          transparent->transparency = 0.33;

        }
        else if (springs_p && (value > 1e-5 || svalue > 1e-5|| tvalue > 1e-5))
        {
          double length   = value;
          double maj_rad1 = value/ratio;
          double maj_rad2 = 0;;
          double min_rad = value/(ratio*ratio);
          unsigned int pitch = (unsigned int) ratio / 3;

          // Secondary value controls the major radius
          if( secondary_major )
            maj_rad1 = svalue;

          // Secondary value controls the minor radius
          if( secondary_minor )
            min_rad = svalue;

          // Secondary value controls the pitch
          if( secondary_pitch )
            pitch = 1 + (unsigned int) svalue;


          // Tertiary value controls the major radius
          if( tertiary_major )
            maj_rad1 = tvalue;

          // Tertiary value controls the minor radius
          if( tertiary_minor )
            min_rad = tvalue;

          // Tertiary value controls the pitch
          if( tertiary_pitch )
            pitch = 1 + (unsigned int) tvalue;


          if( length > 1e-5 && maj_rad1 > 1e-5 && min_rad > 1e-5 )
          {
            if (color_scheme == 0)
              glyphs->add_helix(p, val,
              maj_rad1, maj_rad2, min_rad, length, pitch,
              resolution);
            else if (color_scheme == 2)
              glyphs->add_helix(p, val,
              maj_rad1, maj_rad2, min_rad, length, pitch,
              vcol, resolution);
            else //if (color_scheme == 1)
            {
              if( secondary_alpha )
                glyphs->add_helix(p, val,
                  maj_rad1, maj_rad2, min_rad, length, pitch,
                  vcol, scol, resolution);
              else
                glyphs->add_helix(p, val,
                  maj_rad1, maj_rad2, min_rad, length, pitch,
                  scol, resolution);
            }

            if (bidirectional)
            {
              if (color_scheme == 0)
                glyphs->add_helix(p-val, val,
                  maj_rad2, maj_rad1, min_rad, length, pitch,
                  resolution);
              else if (color_scheme == 2)
                glyphs->add_helix(p-val, val,
                  maj_rad2, maj_rad1, min_rad, length, pitch,
                  vcol, resolution);
              else //if (color_scheme == 1)
              {
                if( secondary_alpha )
                  glyphs->add_helix(p-val, val,
                        maj_rad2, maj_rad1, min_rad, length, pitch,
                        vcol, scol, resolution);
                else
                  glyphs->add_helix(p-val, val,
                        maj_rad2, maj_rad1, min_rad, length, pitch,
                        scol, resolution);
              }
            }
          }

          // Wire frame - render as lines
          else if( length > 1e-5 && maj_rad1 > 1e-5 )
          {
            //Bring nu to expected value for shape.
            unsigned int nu = resolution + 1;
            
            SinCosTable tab1(nu, 0, 2*M_PI);
          
            Transform trans;
            Transform rotate;

            GeomGlyphBase::gen_transforms( p, val, trans, rotate );

            //Single direction only.
            double radius = maj_rad1;

            // The total number of sections is based on the pitch times the
            // number of sections in each rotation which is nu-1.
            double nsections = (double) (pitch * (nu-1));
            
            double dz =  length / nsections;
            double dr = (maj_rad2-maj_rad1) / nsections;
            
            for (unsigned int p=0; p<pitch; p++)
            {
              for (unsigned int u=0; u<nu-1; u++)
              {
                double section = p * (nu-1) + u;

                double z1 = dz *  section;
                double z2 = dz * (section+1);

                double radius1 = radius + dr *  section;
                double radius2 = radius + dr * (section+1);

                double x1 = tab1.sin(u) * radius1;
                double y1 = tab1.cos(u) * radius1;

                double x2 = tab1.sin(u+1) * radius2;
                double y2 = tab1.cos(u+1) * radius2;
              
                Point p1 = trans * Point(x1, y1, z1);
                Point p2 = trans * Point(x2, y2, z2);

                if (color_scheme == 0)
                  lines->add(p1, p2);
                else if (color_scheme == 2)
                  lines->add(p1, vcol, p2, vcol);
                else //if (color_scheme == 1)
                {
                  if( secondary_alpha )
                    lines->add(p1, vcol, scol, p2, vcol, scol);
                  else
                    lines->add(p1, scol, p2, scol);
                }
              }
            }
          }

          // No major radius - render as disks.
          else if( length > 1e-5 && min_rad > 1e-5 )
          {
            if (color_scheme == 0)
              glyphs->add_capped_cylinder(p, val, min_rad, min_rad, length,
                  resolution);
            else if (color_scheme == 2)
              glyphs->add_capped_cylinder(p, val,
                  min_rad, min_rad, length, vcol,
                  resolution);
            else //if (color_scheme == 1)
            {
              if( secondary_alpha )
                glyphs->add_capped_cylinder(p, val,
                    min_rad, min_rad, length, vcol, scol,
                    resolution);
              else
                glyphs->add_capped_cylinder(p, val,
                    min_rad, min_rad, length, scol,
                    resolution);
            }

            if (bidirectional)
            {
              if (color_scheme == 0)
                glyphs->add_capped_cylinder(p-val, val, min_rad, min_rad, length,
                    resolution);
              else if (color_scheme == 2)
                glyphs->add_capped_cylinder(p-val, val,
                    min_rad, min_rad, length, vcol,
                    resolution);
              else //if (color_scheme == 1)
              {
                if( secondary_alpha )
                  glyphs->add_capped_cylinder(p-val, val,
                      min_rad, min_rad, length, vcol, scol,
                      resolution);
                else
                  glyphs->add_capped_cylinder(p-val, val,
                      min_rad, min_rad, length, scol,
                      resolution);
              }
            }

          }

          // No Length - render as rings
          else if( maj_rad1 > 1e-5 && min_rad > 1e-5 )
          {
            if (color_scheme == 0)
              glyphs->add_torus(p, val, maj_rad1, min_rad, resolution);
            else if (color_scheme == 2)
              glyphs->add_torus(p, val, maj_rad1, min_rad, vcol, resolution);
            else //if (color_scheme == 1)
            {
              if( secondary_alpha )
                glyphs->add_torus(p, val, maj_rad1, min_rad, vcol, scol, resolution);
              else
                glyphs->add_torus(p, val, maj_rad1, min_rad, scol, resolution);
            }
          }

          // No major or minor radius - render as lines.
          else if( length > 1e-5 )
          {
            if (bidirectional)
            {
              if (color_scheme == 0)
                lines->add(p - val, p + val);
              else if (color_scheme == 2)
                lines->add(p - val, vcol, p + val, vcol);
              else //if (color_scheme == 1)
              {
          if( secondary_alpha )
            lines->add(p - val, vcol, scol, p + val, vcol, scol);
          else
            lines->add(p - val, scol, p + val, scol);
              }
            }
            else
            {
              if (color_scheme == 0)
                lines->add(p, p + val);
              else if (color_scheme == 2)
                lines->add(p, vcol, p + val, vcol);
              else //if (color_scheme == 1)
              {
                if( secondary_alpha )
                  lines->add(p, vcol, scol, p + val, vcol, scol);
                else
                  lines->add(p, scol, p + val, scol);
              }
            }
          }

          // No Length or minor radius - render as wireframe rings
          else if( maj_rad1 > 1e-5 )
          {
            //Bring nu to expected value for shape.
            unsigned int nu = resolution + 1;

            SinCosTable tab1(nu, 0, 2*M_PI);

            Transform trans;
            Transform rotate;

            GeomGlyphBase::gen_transforms( p, val, trans, rotate );

            double radius = maj_rad1;
            
            for (unsigned int u=0; u<nu-1; u++)
            {
              double x1 = tab1.sin(u) * radius;
              double y1 = tab1.cos(u) * radius;
              
              double x2 = tab1.sin(u+1) * radius;
              double y2 = tab1.cos(u+1) * radius;
              
              Point p1 = trans * Point(x1, y1, 0);
              Point p2 = trans * Point(x2, y2, 0);
              
              if (color_scheme == 0)
                lines->add(p1, p2);
              else if (color_scheme == 2)
                lines->add(p1, vcol, p2, vcol);
              else //if (color_scheme == 1)
              {
                if( secondary_alpha )
                  lines->add(p1, vcol, scol, p2, vcol, scol);
                else
                  lines->add(p1, scol, p2, scol);
              }
            }
          }

          // No Length or major radius - render as sphere
          else if( min_rad > 1e-5 )
          {
            if (color_scheme == 0)
              glyphs->add_sphere(p, min_rad, resolution, resolution);
            else if (color_scheme == 2)
              glyphs->add_sphere(p, min_rad, vcol, resolution, resolution);
            else // if (color_scheme == 1)
            {
              if( secondary_alpha )
                glyphs->add_sphere(p, min_rad, vcol, scol, resolution, resolution);
              else
                glyphs->add_sphere(p, min_rad, scol, resolution, resolution);
            }
          }
        }
        else if (disks_p)
        {
          if (color_scheme == 0)
            glyphs->add_capped_cylinder(p, val, value, value, svalue,
                resolution);
          else if (color_scheme == 2)
            glyphs->add_capped_cylinder(p, val, value, value, svalue, vcol,
                resolution);
          else //if (color_scheme == 1)
          {
            if( secondary_alpha )
              glyphs->add_capped_cylinder(p, val, value, value, svalue, vcol, scol,
                  resolution);
            else
              glyphs->add_capped_cylinder(p, val, value, value, svalue, scol,
                  resolution);
          }
        }
        else if (rings_p && (value > 1e-5 || svalue > 1e-5))
        {
          // Normal ring
          if( value > 1e-5 && svalue > 1e-5 )
          {
            if (color_scheme == 0)
              glyphs->add_torus(p, val, value, svalue, resolution);
            else if (color_scheme == 2)
              glyphs->add_torus(p, val, value, svalue, vcol, resolution);
            else //if (color_scheme == 1)
            {
              if( secondary_alpha )
                glyphs->add_torus(p, val, value, svalue, vcol, scol, resolution);
              else
                glyphs->add_torus(p, val, value, svalue, scol, resolution);
            }
          }

          // Wire frame ring - render a series of lines.
          else if( value > 1e-5 )
          {
            //Bring nu to expected value for shape.
            unsigned int nu = resolution + 1;

            SinCosTable tab1(nu, 0, 2*M_PI);

            Transform trans;
            Transform rotate;

            GeomGlyphBase::gen_transforms( p, val, trans, rotate );

            double radius = value;
            
            for (unsigned int u=0; u<nu-1; u++)
            {
              double x1 = tab1.sin(u) * radius;
              double y1 = tab1.cos(u) * radius;
              
              double x2 = tab1.sin(u+1) * radius;
              double y2 = tab1.cos(u+1) * radius;
              
              Point p1 = trans * Point(x1, y1, 0);
              Point p2 = trans * Point(x2, y2, 0);
              
              if (color_scheme == 0)
                lines->add(p1, p2);
              else if (color_scheme == 2)
                lines->add(p1, vcol, p2, vcol);
              else //if (color_scheme == 1)
              {
                if( secondary_alpha )
                  lines->add(p1, vcol, scol, p2, vcol, scol);
                else
                  lines->add(p1, scol, p2, scol);
              }
            }
          }

          // No radius but thickness - render as a sphere.
          else if( svalue > 1e-5)
          {
            if (color_scheme == 0)
              glyphs->add_sphere(p, svalue, resolution, resolution);
            else if (color_scheme == 2)
              glyphs->add_sphere(p, svalue, vcol, resolution, resolution);
            else // if (color_scheme == 1)
            {
              if( secondary_alpha )
                glyphs->add_sphere(p, svalue, vcol, scol, resolution, resolution);
              else
                glyphs->add_sphere(p, svalue, scol, resolution, resolution);
            }
          }
        }
        else if (lines_p)
        {
          // Readjust for the scaling of arrows 1/2 line 1/2 arrow.
          if(arrows_p)
            val *= 2.0;

          if (bidirectional)
          {
            if (color_scheme == 0)
              lines->add(p - val, p + val);
            else if (color_scheme == 2)
              lines->add(p - val, vcol, p + val, vcol);
            else //if (color_scheme == 1)
            {
              if( secondary_alpha )
                lines->add(p - val, vcol, scol, p + val, vcol, scol);
              else
                lines->add(p - val, scol, p + val, scol);
            }
          }
          else
          {
            if (color_scheme == 0)
              lines->add(p, p + val);
            else if (color_scheme == 2)
              lines->add(p, vcol, p + val, vcol);
            else //if (color_scheme == 1)
            {
              if( secondary_alpha )
                lines->add(p, vcol, scol, p + val, vcol, scol);
              else
                lines->add(p, scol, p + val, scol);
            }
          }
        }
        else  if (needles_p)
        {
          if (color_scheme == 0)
          {
            lines->add(p, opaque, p + val, transparent);
            if (bidirectional)
              lines->add(p, opaque, p - val, transparent);
          }
          else if (color_scheme == 2)
          {
            transparent->diffuse = vcol->diffuse;
            
            lines->add(p, vcol, p + val, transparent);
            if (bidirectional)
              lines->add(p, vcol, p - val, transparent);
          }
          else //if (color_scheme == 1)
          {
            if( secondary_alpha )
            {
              lines->add(p, vcol, scol, p + val, transparent, scol);
              if (bidirectional)
                lines->add(p, vcol, scol, p - val, transparent, scol);
            }
            else
            {
              lines->add(p, opaque, scol, p + val, transparent, scol);
              if (bidirectional)
                lines->add(p, opaque, scol, p - val, transparent, scol);
            }
          }
        }
        else // Too small render as a point
        {
          if (small_is_dot)
          {
            if (color_scheme == 0)
              points->add(p);
            else if (color_scheme == 2)
              points->add(p, vcol);
            else // if (color_scheme == 1)
            {
              if( secondary_alpha )
                points->add(p, vcol, scol);
              else
                points->add(p, scol);
            }
          }
        }
      }
      else // Too small render as a point
      {
        if (small_is_dot)
        {
          if (color_scheme == 0)
            points->add(p);
          else if (color_scheme == 2)
            points->add(p, vcol);
          else // if (color_scheme == 1)
          {
            if( secondary_alpha )
              points->add(p, vcol, scol);
            else
              points->add(p, scol);
          }
        }
      }
    }
  }
  return data_switch;
}


GeomHandle 
RenderTensorField::render_data(FieldHandle pfld_handle,
                               FieldHandle sfld_handle,
                               FieldHandle tfld_handle,
                               const std::string &display_mode,
                               double pscale,
                               double /*sscale*/,
                               double /*tscale*/,
                               double threshold,
                               int resolution,
                               double emphasis,
                               unsigned int render_state,
                               unsigned int secondary_render_state,
                               unsigned int tertiary_render_state)
{
  VField *pfld = 0; if (pfld_handle.get_rep()) pfld = pfld_handle->vfield(); 
  VField *sfld = 0; if (sfld_handle.get_rep()) sfld = sfld_handle->vfield(); 
  VField *tfld = 0; if (tfld_handle.get_rep()) tfld = tfld_handle->vfield(); 
  
  // For superquadrics if the emphasis is near 0 or near 1 use
  // ellipsoids or boxes respectively as they are quicker.
  const bool boxes_p      = (display_mode == "Boxes") ||
    ((display_mode == "Superquadrics") && emphasis > 0.98);
  const bool cboxes_p     = (display_mode == "Colored Boxes");
  const bool ellipsoids_p = (display_mode == "Ellipsoids") ||
    ((display_mode == "Superquadrics") && emphasis < 0.02);
  const bool squadrics_p  = (display_mode == "Superquadrics") &&
    0.02 <= emphasis && emphasis <= 0.98;

  const double mapped_emphasis = map_emphasis(emphasis);

  GeomPoints *points = 0;
  GeomLines *lines = 0;
  GeomGlyphBase *glyphs = 0;
  GeomGroup *grp = new GeomGroup();
  GeomHandle data_switch = new GeomDL(grp);

  // Points when too small
  if (get_flag(render_state, USE_TRANSPARENCY) )
    points = new GeomTranspPoints();
  else
    points = new GeomPoints();

  grp->add(points);

  points->setPointSize(resolution/5.0);

  // Lines when too small
  if( get_flag(render_state, USE_TRANSPARENCY))
    lines = new GeomTranspLines();
  else
    lines = new GeomLines();

  grp->add(lines);
  
  lines->setLineWidth(resolution/5.0);

  // Glyphs otherwise
  if( get_flag(render_state, USE_TRANSPARENCY) )
    glyphs = new GeomTranspGlyph();
  else 
    glyphs = new GeomGlyph();

  grp->add(glyphs->getObj());

  if( cboxes_p || boxes_p )
    pscale *= 2.0;

  bool normalize       = get_flag(render_state, NORMALIZE_DATA);
  bool small_is_dot    = get_flag(render_state,SMALL_IS_DOT);
  
  bool secondary_on    = get_flag(secondary_render_state, IS_ON);
  bool secondary_color = (get_flag(secondary_render_state, IS_ON|USE_COLORMAP) ||
			  get_flag(secondary_render_state, IS_ON|USE_COLOR_CONVERT));
  bool secondary_alpha = get_flag(secondary_render_state, IS_ON|USE_ALPHA) &&
    get_flag(render_state, USE_TRANSPARENCY);
  //  bool secondary_small_is_dot = get_flag(secondary_render_state,SMALL_IS_DOT);

  bool tertiary_on    = get_flag(tertiary_render_state, IS_ON);
  bool tertiary_color = (get_flag(tertiary_render_state, IS_ON|USE_COLORMAP) ||
			 get_flag(tertiary_render_state, IS_ON|USE_COLOR_CONVERT));
  bool tertiary_alpha = (get_flag(tertiary_render_state, IS_ON|USE_ALPHA) &&
			 get_flag(render_state, USE_TRANSPARENCY));
  //  bool tertiary_value = get_flag(tertiary_render_state, IS_ON|USE_VALUE);
  //  bool tertiary_small_is_dot = get_flag(tertiary_render_state,SMALL_IS_DOT);

  unsigned int color_scheme = 0;
  double scol;
  MaterialHandle vcol(0);

  if( cboxes_p ||
      (!tertiary_color && !secondary_color &&
       get_flag(render_state, USE_DEFAULT_COLOR)) )
  {
    color_scheme = 0; // Default color
  }
  else if ( (!tertiary_color && !secondary_color && get_flag(render_state, USE_COLORMAP)) ||
	    ( secondary_color && get_flag(secondary_render_state, USE_COLORMAP)) ||
	    ( tertiary_color && get_flag(tertiary_render_state, USE_COLORMAP)) )
  {
    // Color map lookup using either a scalar value
    // or the vector magnitude.
    color_scheme = 1;
  }
  else if( pfld->basis_order() >= 0 ||
           (secondary_color && sfld->basis_order() >= 0) ||
           (tertiary_color  && tfld->basis_order() >= 0) )
  {
    color_scheme = 2; // Values become RGB
  }

  if( color_scheme == 2 || secondary_alpha || tertiary_alpha )
  {
    vcol = new Material(Color(1.0, 1.0, 1.0));
    if (get_flag(render_state, USE_TRANSPARENCY) )
      vcol->transparency = 0.75;
    else
      vcol->transparency = 1.0;
  }

  VMesh* mesh = pfld->vmesh();

  VMesh::size_type num_values = pfld->num_values();

  for (VMesh::index_type idx=0; idx<num_values;idx++)
  {
    Point p;
    if (pfld->basis_order() == 0)
      mesh->get_center(p, VMesh::Elem::index_type(idx));
    else
      mesh->get_center(p, VMesh::Node::index_type(idx));
      
      
    Tensor val;
    SVTvalue sval;
    SVTvalue tval;

    pfld->get_value(val, idx);
    {
      if( secondary_on )
        get_value(sfld, sval, idx);
      
      if( tertiary_on )
        get_value(tfld, tval, idx);
      
      if (color_scheme)
      {
        if( tertiary_color )
          value_to_color( color_scheme, tval, scol, vcol );
        else if( secondary_color )
          value_to_color( color_scheme, sval, scol, vcol );
        else
          value_to_color( color_scheme,  val, scol, vcol );
      }

      if( tertiary_alpha )
      {
        double atmp;
        to_double(tval, atmp);
        vcol->transparency = fabs(atmp);
      }
      else if( secondary_alpha )
      {
        double atmp;
        to_double(sval, atmp);
        vcol->transparency = fabs(atmp);
      }
      // The color is done so the eigen values can be changed
      Vector e1, e2, e3;
      val.get_eigenvectors(e1, e2, e3);

      static const Point origin(0.0, 0.0, 0.0);
      Transform trans(origin, e1, e2, e3);

      double v1, v2, v3;
      val.get_eigenvalues(v1, v2, v3);
      Vector scales(fabs(v1), fabs(v2), fabs(v3));
      if (normalize)
        scales.safe_normalize();

      scales *= pscale;

      // Do not render tensors that are too small - because surfaces
      // are not rendered at least two of the scales must be non zero.
      if ((scales.x() > threshold && scales.y() > threshold) ||
        (scales.y() > threshold && scales.z() > threshold) ||
        (scales.z() > threshold && scales.x() > threshold))
      {
        // This gives some thickness so there is no z fighting
        if ( scales.x() < threshold ) scales.x(threshold);
        if ( scales.y() < threshold ) scales.y(threshold);
        if ( scales.z() < threshold ) scales.z(threshold);

        if( boxes_p || cboxes_p )
        {
          if (color_scheme == 0)
            glyphs->add_box(p, trans, scales.x(), scales.y(), scales.z(),
                cboxes_p);
          else if (color_scheme == 2)
            glyphs->add_box(p, trans, scales.x(), scales.y(), scales.z(),
                vcol);
          else //if (color_scheme == 1)
          {
            if( secondary_alpha )
              glyphs->add_box(p, trans, scales.x(), scales.y(), scales.z(),
                  vcol, scol);
            else
              glyphs->add_box(p, trans, scales.x(), scales.y(), scales.z(),
                  scol);
          }
        }
        else if (ellipsoids_p)
        {
          if (color_scheme == 0)
            glyphs->add_ellipsoid(p, trans, scales,
                resolution, resolution);
          else if (color_scheme == 2)
            glyphs->add_ellipsoid(p, trans, scales,
                vcol, resolution, resolution);
          else // if (color_scheme == 1)
          {
            if( secondary_alpha )
              glyphs->add_ellipsoid(p, trans, scales,
                  vcol, scol, resolution, resolution);
            else
              glyphs->add_ellipsoid(p, trans, scales,
                  scol, resolution, resolution);
          }
        }
        else if (squadrics_p)
        {
          const double cl = (v1 - v2) / (v1 + v2 + v3);
          const double cp = 2.0 * (v2 - v3) / (v1 + v2 + v3);

          double qA, qB;
          int axis;

          // This is for swapping the axis and qA and qB. But if the
          // axis is swapped and qA and qB it is the same as not doing
          // it at all.
      // 	  if (cl > cp)
      // 	  {
      // 	    axis = 0;
      // 	    qA = pow((1.0 - cp), mapped_emphasis);
      // 	    qB = pow((1.0 - cl), mapped_emphasis);
      // 	  }
      // 	  else
          {
            axis = 2;
            qA = pow((1.0 - cl), mapped_emphasis);
            qB = pow((1.0 - cp), mapped_emphasis);
          }

          if (color_scheme == 0)
            glyphs->add_superquadric(p, trans, scales, qA, qB, axis,
                   resolution);

          else if (color_scheme == 2)
            glyphs->add_superquadric(p, trans, scales, qA, qB, axis,
                   vcol, resolution);

          else //if (color_scheme == 1)
          {
            if( secondary_alpha )
              glyphs->add_superquadric(p, trans, scales, qA, qB, axis,
                     vcol, scol, resolution);
            else
              glyphs->add_superquadric(p, trans, scales, qA, qB, axis,
                     scol, resolution);
          }
        }
      }

      // Tensor as a line
      else if ((scales.x() > threshold && scales.y() < threshold && scales.z() < threshold) ||
	       (scales.x() < threshold && scales.y() > threshold && scales.z() < threshold) ||
	       (scales.x() < threshold && scales.y() < threshold && scales.z() > threshold) )
      {
        Point p1 = p + trans *  scales/ ((cboxes_p || boxes_p) ? 4.0 : 2.0 );
        Point p2 = p + trans * -scales/ ((cboxes_p || boxes_p) ? 4.0 : 2.0 );

        if (color_scheme == 0)
          lines->add(p1, p2);
        else if (color_scheme == 2)
          lines->add(p1, vcol, p2, vcol);
        else //if (color_scheme == 1)
        {
          if( secondary_alpha || tertiary_alpha )
            lines->add(p1, vcol, scol, p2, vcol, scol);
          else
            lines->add(p1, scol, p2, scol);
        }
      }

      else // Too small render as a point
      {
        if (small_is_dot)
        {
          if (color_scheme == 0)
            points->add(p);
          else if (color_scheme == 2)
            points->add(p, vcol);
          else if (color_scheme == 1)
            points->add(p, scol);
       }
      }
    }
  }
  return data_switch;
}



} // end namespace SCIRun
