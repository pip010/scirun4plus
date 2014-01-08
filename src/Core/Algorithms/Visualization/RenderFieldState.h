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
 *  RenderFieldState.h: Rendering alogrithms for data
 *
 *  Written by:
 *   Allen R. Sanderson
 *   SCI Institute
 *   University of Utah
 *   April 2007
 *
 */

#if !defined(Visualization_RenderFieldState_h)
#define Visualization_RenderFieldState_h

#include <Core/Algorithms/Visualization/DataConversions.h>

#include <Core/Geom/GeomMaterial.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Tensor.h>

#if defined(_WIN32) && !defined(uint)
  // for some reason, this isn't defined...
#define uint unsigned
#endif


namespace SCIRun {

template <class Val_t >
void value_to_color( unsigned int color_scheme,
         Val_t val,
         double &scol,
         MaterialHandle &vcol )
{
  // Set the color based on the color table look up
  if (color_scheme == 1)
    to_double(val, scol);
  // Set the color based on the value.
  else if (color_scheme == 2) {
    to_color(val, vcol->ambient);
    to_color(val, vcol->diffuse);
  }
}


class RenderStateBase
{
public:
  enum action_flag_e {
    CLEAR_FLAGS      = 0x0000,

    IS_ON             = 0x0001,
    HAS_DATA          = 0x0002,
    USE_TRANSPARENCY  = 0x0004,
    NORMALIZE_DATA    = 0x0008,

    USE_DEFAULT_COLOR = 0x0010,
    USE_COLORMAP      = 0x0020,
    USE_COLOR_CONVERT = 0x0040,
  
    SMALL_IS_DOT      = 0x0080,

    DIRTY             = 0x0800,

    // Node Flags

    // Edge Flags

    // Face Flags
    USE_NORMALS = 0x1000,
    USE_TEXTURE = 0x2000,

    // Text Flags
    BACKFACES_CULL = 0x1000,
    ALWAYS_VISIBLE = 0x2000,

    // Scalar Data Flags

    // Vector Data Flags
    BIDIRECTIONAL = 0x1000,

    // Tensor Data Flags

    // Secondary/Tertiary Data Flags
    USE_ALPHA        = 0x01000,
    USE_VALUE        = 0x02000,
    USE_MAJOR_RADIUS = 0x04000,
    USE_MINOR_RADIUS = 0x08000,
    USE_PITCH        = 0x10000
  };

  enum glyph_type_e {
    POINT_GLYPH     = 0x0001,
    SPHERE_GLYPH    = 0x0002,
    BOX_GLYPH       = 0x0004,
    AXIS_GLYPH      = 0x0008,
    LINE_GLYPH      = 0x0010,
    NEEDLE_GLYPH    = 0x0020,
    ARROW_GLYPH     = 0x0040,
    CONE_GLYPH      = 0x0080,
    DISK_GLYPH      = 0x0100,
    RING_GLYPH      = 0x0200
  };

protected:



  std::string clean_fieldname(const std::string& fname)
  {
    std::string result;
    int counter = 0;

    for (size_t i = 0; i < fname.size(); i++)
    {
      if (fname[i] == ':')
      {
        if (counter) {
          // do nothing
        } else {
          result += fname[i];
          counter = 1;
        }
      }
      else
      {
        result += fname[i];
        counter = 0;
      }
    }

    return result;
  }
};

#define set_flag( value, flag )	\
    ( value |= (flag) )

#define clear_flag( value, flag ) \
    ( value &= ~(flag) )

#define switch_flag( value, flag, boolean ) \
    if( boolean ) value |= (flag); else value &= ~(flag);

#define get_flag( value, flag ) \
    ( (value & (flag)) == (flag) )

} // end namespace SCIRun

#endif // Visualization_RenderState_h
