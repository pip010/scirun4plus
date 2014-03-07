//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  
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
//    File   : CM2Widget.cc
//    Author : Milan Ikits
//    Date   : Mon Jul  5 18:33:29 2004


#include <Core/Volume/CM2Widget.h>


#include <Core/Datatypes/Color.h>
#include <Core/Datatypes/ColorMap.h>

#include <sci_gl.h>

#include <iostream>
#include <sstream>


#include <math.h>
#include <stdlib.h>


namespace SCIRun {


using SLIVR::CM2RectangleType;

PersistentTypeID CM2Widget::type_id("CM2Widget", "Datatype", 0);

CM2Widget::CM2Widget() :
  SLIVR::CM2Widget()
{
}

CM2Widget::CM2Widget(const CM2Widget &copy) : 
  SLIVR::CM2Widget(copy), Datatype(copy)
{
}

CM2Widget::~CM2Widget()
{}

#define CM2WIDGET_VERSION 1

void
CM2Widget::io(Piostream &stream)
{  
  stream.begin_class("CM2Widget", CM2WIDGET_VERSION);

  Pio(stream, name_);
  Pio(stream, color_);
  Pio(stream, alpha_);
  Pio(stream, selected_);
  Pio(stream, shadeType_);
  Pio(stream, onState_);
  Pio(stream, faux_);
  Pio(stream, value_range_);
  stream.end_class();
}

CM2Widget* 
CM2Widget::slivr2sr(SLIVR::CM2Widget *w)
{
  
  SLIVR::ClippingCM2Widget *cw = dynamic_cast<SLIVR::ClippingCM2Widget*>(w);
  SLIVR::RectangleCM2Widget *rw = dynamic_cast<SLIVR::RectangleCM2Widget*>(w);
  SLIVR::TriangleCM2Widget *tw = dynamic_cast<SLIVR::TriangleCM2Widget*>(w);
  SLIVR::EllipsoidCM2Widget *ew = dynamic_cast<SLIVR::EllipsoidCM2Widget*>(w);
  SLIVR::ParaboloidCM2Widget *pw = dynamic_cast<SLIVR::ParaboloidCM2Widget*>(w);
  SLIVR::ColorMapCM2Widget *cmw = dynamic_cast<SLIVR::ColorMapCM2Widget*>(w);
  SLIVR::ImageCM2Widget *iw = dynamic_cast<SLIVR::ImageCM2Widget*>(w);
  SLIVR::PaintCM2Widget *ptw = dynamic_cast<SLIVR::PaintCM2Widget*>(w);

  // check in most commonly used order for conversion.
  if (rw) {
    return new RectangleCM2Widget(rw);
  } else if (tw) {
    return new TriangleCM2Widget(tw);
  } else if (cw) {
    return new ClippingCM2Widget(cw);
  } else if (ew) {
    return new EllipsoidCM2Widget(ew);
  } else if (pw) {
    return new ParaboloidCM2Widget(pw);
  } else if (cmw) {
    return new ColorMapCM2Widget(cmw);
  } else if (iw) {
    return new ImageCM2Widget(iw);
  } else if (ptw) {
    return new PaintCM2Widget(ptw);
  }
  ASSERTFAIL("CM2Widget conversion Failed.")
  return 0;
}

static Persistent* ParaboloidCM2Widget_maker()
{
  return new ParaboloidCM2Widget;
}

PersistentTypeID ParaboloidCM2Widget::type_id("ParaboloidCM2Widget",
					      "CM2Widget",
					      ParaboloidCM2Widget_maker);

#define PARABOLOIDCM2WIDGET_VERSION 1

void
ParaboloidCM2Widget::io(Piostream &stream)
{
  stream.begin_class("ParaboloidCM2Widget", PARABOLOIDCM2WIDGET_VERSION);
  CM2Widget::io(stream);

  Pio(stream, shadeType_);
  Pio(stream, onState_);

  // save/load position attributes
  Pio(stream, top_x_);
  Pio(stream, top_y_);
  Pio(stream, bottom_x_);
  Pio(stream, bottom_y_);
  Pio(stream, left_x_);
  Pio(stream, left_y_);
  Pio(stream, right_x_);
  Pio(stream, right_y_);

  stream.end_class();
}

ParaboloidCM2Widget::ParaboloidCM2Widget() : 
  SLIVR::ParaboloidCM2Widget()
{
}

ParaboloidCM2Widget::ParaboloidCM2Widget(const ParaboloidCM2Widget &copy) :
  SLIVR::CM2Widget(copy),
  SCIRun::CM2Widget(copy),
  SLIVR::ParaboloidCM2Widget(copy)
{
}

ParaboloidCM2Widget::ParaboloidCM2Widget(const SLIVR::ParaboloidCM2Widget *w) : 
  SLIVR::ParaboloidCM2Widget(*w)
{
}

ParaboloidCM2Widget::ParaboloidCM2Widget(float top_x, float top_y,
					 float bottom_x, float bottom_y,
					 float left_x, float left_y,
					 float right_x, float right_y) :
  SLIVR::ParaboloidCM2Widget(top_x, top_y,
			     bottom_x, bottom_y,
			     left_x, left_y,
			     right_x, right_y)
{
}

ParaboloidCM2Widget::~ParaboloidCM2Widget()
{}


ParaboloidCM2Widget*
ParaboloidCM2Widget::clone()
{
  return new ParaboloidCM2Widget(*this);
}




static Persistent* EllipsoidCM2Widget_maker()
{
  return new EllipsoidCM2Widget;
}

PersistentTypeID EllipsoidCM2Widget::type_id("EllipsoidCM2Widget",
					     "CM2Widget",
					     EllipsoidCM2Widget_maker);
#define ELLIPSOIDCM2WIDGET_VERSION 1

void
EllipsoidCM2Widget::io(Piostream &stream)
{
  stream.begin_class("EllipsoidCM2Widget", ELLIPSOIDCM2WIDGET_VERSION);
  CM2Widget::io(stream);

  // save/load position attributes
  Pio(stream, center_x_);
  Pio(stream, center_y_);
  Pio(stream, a_);
  Pio(stream, b_);
  Pio(stream, rot_);

  stream.end_class();
}

EllipsoidCM2Widget::EllipsoidCM2Widget() : 
  SLIVR::EllipsoidCM2Widget()
{
}

EllipsoidCM2Widget::EllipsoidCM2Widget(const EllipsoidCM2Widget &copy) :
  SLIVR::CM2Widget(copy),
  SCIRun::CM2Widget(copy),
  SLIVR::EllipsoidCM2Widget(copy)
{
}

EllipsoidCM2Widget::EllipsoidCM2Widget(const SLIVR::EllipsoidCM2Widget *w) : 
  SLIVR::EllipsoidCM2Widget(*w)
{
}

EllipsoidCM2Widget::EllipsoidCM2Widget(float x, float y, float a, float b, float rot) :
  SLIVR::EllipsoidCM2Widget(x, y, a, b, rot)
{
}

EllipsoidCM2Widget::~EllipsoidCM2Widget()
{}


EllipsoidCM2Widget*
EllipsoidCM2Widget::clone()
{
  return new EllipsoidCM2Widget(*this);
}

static Persistent* RectangleCM2Widget_maker()
{
  return new RectangleCM2Widget;
}

PersistentTypeID RectangleCM2Widget::type_id("RectangleCM2Widget", 
					     "CM2Widget",
                                             RectangleCM2Widget_maker);

#define RECTANGLECM2WIDGET_VERSION 4

void
RectangleCM2Widget::io(Piostream &stream)
{
  const int version = 
    stream.begin_class("RectangleCM2Widget", RECTANGLECM2WIDGET_VERSION);
  // base class slots need to pio (didn't in older versions)
  if (version >= 4) {
    CM2Widget::io(stream);
  }
  // Originally used "Pio(stream, (int)type_);", but this did not
  // compile on the SGI, so needed to do it this way.
  int tmp = (int)type_;
  Pio(stream, tmp);
  if (stream.reading())
  {
    type_ = (CM2RectangleType)tmp;
  }

  Pio(stream, left_x_);
  Pio(stream, left_y_);
  Pio(stream, width_);
  Pio(stream, height_);
  Pio(stream, offset_);

  if (version < 4) {
    Pio(stream, shadeType_);
    Pio(stream, onState_);
    Color c(color_.r(), color_.g(), color_.b());
    Pio(stream, c);
    color_.r(c.r());
    color_.g(c.g());
    color_.b(c.b());
    Pio(stream, alpha_);
  }

  if (version == 2) {
    Pio(stream, name_);
    double temp;
    Pio(stream, temp);
    Pio(stream, temp);
    value_range_.first = 0.0;
    value_range_.second = -1.0;
  }

  if (version == 3) {
    Pio(stream, name_);
    Pio(stream, value_range_.first);
    Pio(stream, value_range_.second);
  }

  stream.end_class();
}

RectangleCM2Widget::RectangleCM2Widget() : 
  SLIVR::RectangleCM2Widget()
{
}

RectangleCM2Widget::RectangleCM2Widget(const RectangleCM2Widget &copy) :
  SLIVR::CM2Widget(copy),
  SCIRun::CM2Widget(copy),
  SLIVR::RectangleCM2Widget(copy)
{
}

RectangleCM2Widget::RectangleCM2Widget(const SLIVR::RectangleCM2Widget *w) : 
  SLIVR::RectangleCM2Widget(*w)
{
}

RectangleCM2Widget::RectangleCM2Widget(CM2RectangleType type, 
				       float left_x, float left_y, 
				       float width, float height,
                                       float offset) : 
  SLIVR::RectangleCM2Widget(type, left_x, left_y, width, height, offset)
{
}

RectangleCM2Widget::~RectangleCM2Widget()
{}


RectangleCM2Widget*
RectangleCM2Widget::clone()
{
  return new RectangleCM2Widget(*this);
}

static Persistent* TriangleCM2Widget_maker()
{
  return new TriangleCM2Widget;
}

PersistentTypeID TriangleCM2Widget::type_id("TriangleCM2Widget", "CM2Widget",
                                            TriangleCM2Widget_maker);

#define TRIANGLECM2WIDGET_VERSION 4

void
TriangleCM2Widget::io(Piostream &stream)
{
  const int version = 
    stream.begin_class("TriangleCM2Widget", TRIANGLECM2WIDGET_VERSION);


  if (version >= 4) {
    CM2Widget::io(stream);
  }

  Pio(stream, base_);
  Pio(stream, top_x_);
  Pio(stream, top_y_);
  Pio(stream, width_);
  Pio(stream, bottom_);

  if (version < 4) {
    Pio(stream, shadeType_);
    Pio(stream, onState_);
    Color c(color_.r(), color_.g(), color_.b());
    Pio(stream, c);
    color_.r(c.r());
    color_.g(c.g());
    color_.b(c.b());
    Pio(stream, alpha_);
  }

  if (version == 2) {
    Pio(stream, name_);
    double temp;
    Pio(stream, temp);
    Pio(stream, temp);
    value_range_.first = 0.0;
    value_range_.second = -1.0;
  }

  if (version == 3) {
    Pio(stream, name_);
    Pio(stream, value_range_.first);
    Pio(stream, value_range_.second);
  }
    
  stream.end_class();
}

TriangleCM2Widget::TriangleCM2Widget() : 
  SLIVR::TriangleCM2Widget()
{
}

TriangleCM2Widget::TriangleCM2Widget(const TriangleCM2Widget &copy) :
  SLIVR::CM2Widget(copy),
  SCIRun::CM2Widget(copy), 
  SLIVR::TriangleCM2Widget(copy)
{
}

TriangleCM2Widget::TriangleCM2Widget(const SLIVR::TriangleCM2Widget *w) : 
  SLIVR::TriangleCM2Widget(*w)
{
}


TriangleCM2Widget::TriangleCM2Widget(float base, float top_x, float top_y,
                                     float width, float bottom) : 
  SLIVR::TriangleCM2Widget(base, top_x, top_y, width, bottom)
{}

TriangleCM2Widget::~TriangleCM2Widget()
{}

TriangleCM2Widget*
TriangleCM2Widget::clone()
{
  return new TriangleCM2Widget(*this);
}


// Image --
#define IMAGECM2WIDGET_VERSION 2

static Persistent* ImageCM2Widget_maker()
{
  return new ImageCM2Widget;
}

PersistentTypeID ImageCM2Widget::type_id("ImageCM2Widget", "CM2Widget",
                                         ImageCM2Widget_maker);

void
ImageCM2Widget::io(Piostream &stream)
{
  int version = stream.begin_class("ImageCM2Widget", IMAGECM2WIDGET_VERSION);

  if (version >= 2) {
    CM2Widget::io(stream);
  }

  NrrdData nd(pixels_);
  Pio(stream, nd);
  pixels_ = nd.nrrd_;
  stream.end_class();
}

ImageCM2Widget::ImageCM2Widget() : 
  SLIVR::ImageCM2Widget()
{}

ImageCM2Widget::ImageCM2Widget(const ImageCM2Widget &copy) :
  SLIVR::CM2Widget(copy),
  SCIRun::CM2Widget(copy), 
  SLIVR::ImageCM2Widget(copy)
{}

ImageCM2Widget::ImageCM2Widget(const SLIVR::ImageCM2Widget *w) : 
  SLIVR::ImageCM2Widget(*w)
{}

ImageCM2Widget::ImageCM2Widget(NrrdDataHandle d) :
  SLIVR::ImageCM2Widget(d->nrrd_)
{}

ImageCM2Widget::~ImageCM2Widget()
{}

ImageCM2Widget*
ImageCM2Widget::clone()
{
  return new ImageCM2Widget(*this);
}

#define PAINTCM2WIDGET_VERSION 3

static Persistent* PaintCM2Widget_maker()
{
  return new PaintCM2Widget;
}

PersistentTypeID PaintCM2Widget::type_id("PaintCM2Widget", "CM2Widget",
                                         PaintCM2Widget_maker);

void
PaintCM2Widget::io(Piostream &stream)
{
  const int version = 
    stream.begin_class("PaintCM2Widget", PAINTCM2WIDGET_VERSION);
  
  if (version >= 3) {
    CM2Widget::io(stream);
  }

  Pio(stream, strokes_);

  if (version <= 3) {
  Pio(stream, name_);
  }
  if (version == 2) {
    Color c(color_.r(), color_.g(), color_.b());
    Pio(stream, c);
    color_.r(c.r());
    color_.g(c.g());
    color_.b(c.b());
    Pio(stream, alpha_);
    Pio(stream, shadeType_);
    Pio(stream, onState_);
    Pio(stream, faux_);
    Pio(stream, value_range_.first);
    Pio(stream, value_range_.second);
  }

  stream.end_class();
}

PaintCM2Widget::PaintCM2Widget() : 
  SLIVR::PaintCM2Widget()
{}

PaintCM2Widget::~PaintCM2Widget()
{}

PaintCM2Widget::PaintCM2Widget(const SLIVR::PaintCM2Widget *w) : 
  SLIVR::PaintCM2Widget(*w)
{}

PaintCM2Widget*
PaintCM2Widget::clone()
{
  return new PaintCM2Widget(*this);
}

static Persistent* ColorMapCM2Widget_maker()
{
  return new ColorMapCM2Widget;
}

PersistentTypeID ColorMapCM2Widget::type_id("ColorMapCM2Widget", 
                                            "RectangleCM2Widget",
                                            ColorMapCM2Widget_maker);

#define COLORMAPCM2WIDGET_VERSION 1

void
ColorMapCM2Widget::io(Piostream &stream)
{
  stream.begin_class("ColorMapCM2Widget", COLORMAPCM2WIDGET_VERSION);

  // Originally used "Pio(stream, (int)type_);", but this did not
  // compile on the SGI, so needed to do it this way.
  int tmp = (int)type_;
  Pio(stream, tmp);
  if (stream.reading())
  {
    type_ = (CM2RectangleType)tmp;
  }

  Pio(stream, left_x_);
  Pio(stream, left_y_);
  Pio(stream, width_);
  Pio(stream, height_);
  Pio(stream, offset_);
  Pio(stream, shadeType_);
  Pio(stream, onState_);
  Color c(color_.r(), color_.g(), color_.b());
  Pio(stream, c);
  color_.r(c.r());
  color_.g(c.g());
  color_.b(c.b());
  Pio(stream, alpha_);

  Pio(stream, name_);
  Pio(stream, value_range_.first);
  Pio(stream, value_range_.second);  


  SCIRun::ColorMap cm(*colormap_);
  Pio(stream, cm);
  if (stream.reading()) {
    SLIVR::ColorMap &cmr = cm;
    *colormap_ = cmr;
  }
  stream.end_class();
}

ColorMapCM2Widget::ColorMapCM2Widget() : 
  SLIVR::ColorMapCM2Widget()
{}

ColorMapCM2Widget::ColorMapCM2Widget(const ColorMapCM2Widget &copy) :
  SLIVR::CM2Widget(copy),
  SLIVR::RectangleCM2Widget(copy),
  SCIRun::CM2Widget(copy), 
  SLIVR::ColorMapCM2Widget(copy)
{}

ColorMapCM2Widget::ColorMapCM2Widget(const SLIVR::ColorMapCM2Widget *w) : 
  SLIVR::ColorMapCM2Widget(*w)
{}

ColorMapCM2Widget::ColorMapCM2Widget(CM2RectangleType type, float left_x, 
                                     float left_y, float width, float height,
                                     float offset) : 
  SLIVR::ColorMapCM2Widget(type, left_x, left_y, width, height, offset)
{}

ColorMapCM2Widget::~ColorMapCM2Widget()
{}

ColorMapCM2Widget*
ColorMapCM2Widget::clone()
{
  return new ColorMapCM2Widget(*this);
}

static Persistent* ClippingCM2Widget_maker()
{
  return new ClippingCM2Widget;
}

PersistentTypeID ClippingCM2Widget::type_id("ClippingCM2Widget", "CM2Widget",
                                            ClippingCM2Widget_maker);

#define CLIPPINGCM2WIDGET_VERSION 1

void
ClippingCM2Widget::io(Piostream &/*stream*/)
{
  printf("WARNING: ClippingCM2Widget::io() not implemented!!!\n");
}

ClippingCM2Widget::ClippingCM2Widget() : 
  SLIVR::ClippingCM2Widget()
{}

ClippingCM2Widget::ClippingCM2Widget(const ClippingCM2Widget &copy) :
  SLIVR::CM2Widget(copy), 
  SCIRun::CM2Widget(copy), 
  SLIVR::ClippingCM2Widget(copy)
{}

ClippingCM2Widget::ClippingCM2Widget(const SLIVR::ClippingCM2Widget *w) : 
  SLIVR::ClippingCM2Widget(*w)
{}

ClippingCM2Widget::~ClippingCM2Widget()
{}

ClippingCM2Widget*
ClippingCM2Widget::clone()
{
  return new ClippingCM2Widget(*this);
}

} // namespace SCIRun
