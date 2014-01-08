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
 *  RenderFieldGlyphs.h: Rendering alogrithms for data
 *
 *  Written by:
 *   Allen R. Sanderson
 *   SCI Institute
 *   University of Utah
 *   April 2007
 *
 */

#if !defined(Visualization_RenderFieldGlyphs_h)
#define Visualization_RenderFieldGlyphs_h

#include <Core/Algorithms/Visualization/RenderFieldState.h>

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Datatype.h>

#include <Core/Geom/GeomDL.h>
#include <Core/Geom/GeomGroup.h>
#include <Core/Geom/GeomSwitch.h>
#include <Core/Geom/GeomTransform.h>
#include <Core/Geom/GeomGlyph.h>
#include <Core/Geom/GeomLine.h>
#include <Core/Geom/GeomPoint.h>
#include <Core/Geom/GeomText.h>
#include <Core/Geom/GeomSphere.h>
#include <Core/Geom/GeomMaterial.h>
#include <Core/Geom/GeomColorMap.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Tensor.h>

#include <Core/Algorithms/Visualization/share.h>

namespace SCIRun {


//! RenderScalarFieldBase supports the dynamically loadable algorithm concept.
//! when dynamically loaded the user will dynamically cast to a 
//! RenderScalarFieldBase from the DynamicAlgoBase they will have a pointer to.
class SCISHARE RenderScalarField : public RenderStateBase, public Datatype
{
  public:

    GeomHandle render_data(FieldHandle pfld_handle,
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
				 unsigned int tertiary_render_state);
    virtual std::string dynamic_type_name() const { return "RenderScalarField"; }

  protected:
    void add_axis(const Point &p, double scale, GeomLines *lines);
    void add_axis(const Point &p, double scale, GeomLines *lines,
      const MaterialHandle &vcol);
    void add_axis(const Point &p, double scale, GeomLines *lines,
      double cindex);
    void add_axis(const Point &p, double scale, GeomLines *lines,
      const MaterialHandle &vcol, double cindex);
};

//! RenderVectorFieldBase supports the dynamically loadable algorithm concept.
//! when dynamically loaded the user will dynamically cast to a 
//! RenderVectorFieldBase from the DynamicAlgoBase they will have a pointer to.
class SCISHARE RenderVectorField : public RenderStateBase, public Datatype
{
  public:

    GeomHandle render_data(FieldHandle pfld_handle,
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
				 unsigned int tertiary_render_state);
    virtual std::string dynamic_type_name() const { return "RenderVectorField"; }
};

//! RenderTensorFieldBase supports the dynamically loadable algorithm concept.
//! when dynamically loaded the user will dynamically cast to a 
//! RenderTensorFieldBase from the DynamicAlgoBase they will have a pointer to.
class SCISHARE RenderTensorField : public RenderStateBase, public Datatype
{
  public:

    GeomHandle render_data(FieldHandle pfld_handle,
                           FieldHandle sfld_handle,
                           FieldHandle tfld_handle,
                           const std::string &display_mode,
                           double pscale,
                           double sscale,
                           double tscale,
                           double threshold,
                           int resolution,
                           double emphasis,
                           unsigned int render_state,
                           unsigned int secondary_render_state,
                           unsigned int tertiary_render_state);
    virtual std::string dynamic_type_name() const { return "RenderTensorField"; }
  protected:
    double map_emphasis(double zero_to_one);
};

} // end namespace SCIRun

#endif // Visualization_RenderField_h
