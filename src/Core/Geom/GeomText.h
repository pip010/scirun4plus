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

#ifndef SCI_Geom_Text_h
#define SCI_Geom_Text_h 1

#include <Core/Geom/GeomObj.h>
#include <Core/Geom/TextRenderer.h>
#include <Core/Geometry/Transform.h>

#include <Core/Datatypes/Color.h>

class FTGLTextureFont;

#include <Core/Geom/share.h>

namespace SCIRun {

class SCISHARE GeomTexts : public GeomObj
{
  protected:
    int fontindex_;
    std::vector<std::string> text_;
    std::vector<Point>  location_;
    std::vector<Color>  color_;
    std::vector<float>  index_;
    std::vector<Vector> normal_;

    TextRendererHandle renderer_;
    bool               disable_depth_test_;
    bool               is_2d_;
    double             offset_;
    double             scale_;

  public:
    GeomTexts();
    GeomTexts(TextRendererHandle& renderer);

    GeomTexts(const GeomTexts &);
    virtual ~GeomTexts();
    virtual GeomObj* clone();

    void set_font_index(int);
    void set_offset(double);

    void add (const std::string &text, const Point &loc);
    void add (const std::string &text, const Point &loc, const Color &c);
    void add (const std::string &text, const Point &loc, float index);

    void add (const std::string &text, const Vector &norm, const Point &loc);
    void add (const std::string &text, const Vector &norm, const Point &loc, const Color &c);
    void add (const std::string &text, const Vector &norm, const Point &loc, float index);

    void clear();
    
    void set_always_visible() { disable_depth_test_ = true; }
    void set_is_2d(bool b) { is_2d_ = b; }

    bool is_2d_p() const { return is_2d_; }

    virtual void reset_bbox();
    virtual void get_bounds(BBox&);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);
};

class SCISHARE GeomTextsCulled : public GeomTexts
{
  protected:
    std::vector<Vector> normal_;

  public:
    GeomTextsCulled();
    GeomTextsCulled(const GeomTextsCulled &);
    virtual ~GeomTextsCulled();
    virtual GeomObj* clone();

    void add (const std::string &text, const Point &loc, const Vector &vec);
    void add (const std::string &text, const Point &loc, const Vector &vec, const Color &c);
    void add (const std::string &text, const Point &loc, const Vector &vec, float index);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);
};

class SCISHARE GeomTextTexture : public GeomObj
{
  public:
    enum anchor_e { n, e, s, w, ne, se, sw, nw, c};
  private:

    void		build_transform(Transform &);
    void		render();
    Transform	transform_;
    std::string	text_;
    Point		origin_;
    Vector	up_;
    Vector	left_;
    Color		color_;
    anchor_e	anchor_;
    bool		own_font_;
    
  public:
    GeomTextTexture(const std::string &fontfile);
    GeomTextTexture(const GeomTextTexture&);
    GeomTextTexture(const std::string &fontfile,
        const std::string &text, 
        const Point &at,
        const Vector &up,
        const Vector &left,
        const Color &c = Color(1.,1.,1.));

    virtual ~GeomTextTexture();
    virtual GeomObj* clone();
    
    virtual void get_bounds(BBox&);

    void	set_anchor(anchor_e anchor) { anchor_ = anchor; }
    
    virtual void draw(DrawInfoOpenGL*, Material*, double time);

    bool		up_hack_;
};

} // End namespace SCIRun

#endif /* SCI_Geom_Text_h */
