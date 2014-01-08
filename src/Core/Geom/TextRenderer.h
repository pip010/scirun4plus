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

#ifndef CORE_GEOM_TEXTRENDERER_H
#define CORE_GEOM_TEXTRENDERER_H 1

#include <Core/Geom/FreeType.h>
#include <Core/Geom/View.h>
#include <Core/Geometry/BBox.h>
#include <Core/Containers/Handle.h>
#include <Core/Geom/share.h>

#include <string>
#include <map>
#include <set>

namespace SCIRun {

class TextureObj;

class SCISHARE TextRenderer {

  public:
    // Reference counting, so the object can be destroyed when it is no longer
    // in use. NOTE: currently TextRenderers cannot be shared by multiple threads
    int ref_cnt;

    // For SCIRun, thread safe texture renderers
    TextRenderer();
    TextRenderer(const std::string& name);

    // For Seg3D
    TextRenderer(FreeTypeFace *face);
    virtual ~TextRenderer();

    enum flags_e { N  = 1,
                   E  = 2, 
                   S  = 3,
                   W  = 4,
                   NE = 5,
                   SE = 6,
                   SW = 7,
                   NW = 8,
                   C  = 9,
                   ANCHOR_MASK = 15, 
                   VERTICAL = 16,
                   SHADOW = 32,
                   REVERSE = 64,
                   EXTRUDED = 128,
                   CURSOR = 256};

    int     width(const std::string &text, int flags = 0);
    int     height(const std::string &text, int flags = 0);

    void    set_color(float, float, float, float);
    void    set_color(float color[4]);

    void    set_shadow_color(float, float, float, float);
    void    set_shadow_color(float color[4]);

    void    set_shadow_offset(int x, int y);

    void    set_default_flags(int);

    void		render(const std::string &, float x, float y, int flags);

    // rendering in a 3d scene.
    void		render(const std::string &, float x, float y, float z, 
                                 const View &view, double sf, int flags);

    void		render(const std::string &, Point& anchor, Vector& normal,
                                double sf, int flags);

    double  find_scale_factor(const std::vector<Point> &p, const View &v) const;
    
    void    set_cursor_position(unsigned int pos);
    
    void    set_size(double size);
    void    set_offset(double offset);
    
  private:
    // Info for glyphs rendered to texture
    struct GlyphInfo 
    {
      FT_Glyph_Metrics    ft_metrics_;
      FT_UInt             index_;
      TextureObj *        texture_;
      float               tex_coords_[8];
    };

    typedef std::map<wchar_t, GlyphInfo *> GlyphMap;
    
    // Info for rendering glyphs from texture to screen
    struct GlyphLayoutInfo {
      GlyphInfo *         glyph_info_;
      Point               vertices_[4];
      float               color_[4];
    };

    typedef std::vector<GlyphLayoutInfo*> GlyphLayoutVector;

    struct TextLayoutInfo {
      std::string         text_;
      GlyphLayoutVector   layout_vector_;
      BBox                bbox_;
      Vector              offset_;
      Point               anchor_;
      Vector              left_;
      Vector              up_;
      int                 flags_;
      unsigned int        total_glyphs_;    
    };

    
    FreeTypeFaceHandle    face_;
    GlyphMap              glyphs_;

    std::map<std::string, std::map<int, TextLayoutInfo> > layouts_;

    float                 color_[4];
    float                 shadow_color_[4];
    std::pair<int, int>   shadow_offset_;  
    unsigned int          cursor_position_;

    // Textures of all glyphs
    typedef std::set<TextureObj *> TextureSet;
    TextureSet            textures_;

    // Texture to render next glyph to
    TextureObj *          texture_;

    // Info on where to render next glyph
    int                   x_;
    int                   y_;
    int                   height_;
    double                offset_;

    GlyphInfo *           render_glyph_to_texture(const wchar_t &);
    void                  render_string_glyphs_to_texture(const std::string &);
    void                  layout_text(TextLayoutInfo &, bool clamp_vals = false);
    void                  move_text_to_anchor(TextLayoutInfo &, int);
    
    // Get default SCIRun font
    std::string get_default_font_name();
};

// Define a handle to this object
typedef Handle<TextRenderer> TextRendererHandle;

}


#endif
