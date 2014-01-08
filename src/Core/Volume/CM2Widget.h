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
//    File   : CM2Widget.h
//    Author : Milan Ikits
//    Date   : Mon Jul  5 18:33:12 2004

#ifndef CM2Widget_h
#define CM2Widget_h

#include <slivr/CM2Widget.h>

#include <Core/Datatypes/Datatype.h>
#include <Core/Containers/LockingHandle.h>
#include <Core/Persistent/Persistent.h>
#include <Core/Persistent/PersistentSTL.h>
#include <Core/Datatypes/NrrdData.h>

#include <string>

#include <Core/Volume/share.h>

namespace SCIRun {

class SCISHARE CM2Widget : public virtual SLIVR::CM2Widget, public Datatype
{
public:
  CM2Widget();
  CM2Widget(const CM2Widget &copy);
  virtual ~CM2Widget();  
  virtual CM2Widget*    clone() = 0;

  static CM2Widget* slivr2sr(SLIVR::CM2Widget*);
  
  // State management to/from disk
  virtual void          io(Piostream &stream);
  static PersistentTypeID type_id;
  virtual std::string dynamic_type_name() const { return type_id.type; }
};

typedef LockingHandle<CM2Widget> CM2WidgetHandle;

class SCISHARE ClippingCM2Widget : 
    public CM2Widget, public SLIVR::ClippingCM2Widget
{
public:
  ClippingCM2Widget();
  ClippingCM2Widget(const ClippingCM2Widget &copy);
  ClippingCM2Widget(const SLIVR::ClippingCM2Widget *);
  ~ClippingCM2Widget();

  virtual ClippingCM2Widget*    clone();
  
  virtual void          io(Piostream &stream);
  static PersistentTypeID type_id;
  virtual std::string dynamic_type_name() const { return type_id.type; }
};

class SCISHARE TriangleCM2Widget : 
  public CM2Widget, public SLIVR::TriangleCM2Widget
{
public:
  TriangleCM2Widget();
  TriangleCM2Widget(const TriangleCM2Widget &copy);
  TriangleCM2Widget(const SLIVR::TriangleCM2Widget *);
  TriangleCM2Widget(float base, float top_x, float top_y,
                    float width, float bottom);
  ~TriangleCM2Widget();
  virtual TriangleCM2Widget*    clone();
  
  virtual void          io(Piostream &stream);
  static PersistentTypeID type_id;
  virtual std::string dynamic_type_name() const { return type_id.type; }
};

class SCISHARE RectangleCM2Widget : 
  public CM2Widget, public SLIVR::RectangleCM2Widget
{
public:
  RectangleCM2Widget();
  RectangleCM2Widget(const RectangleCM2Widget &copy);
  RectangleCM2Widget(const SLIVR::RectangleCM2Widget *);
  RectangleCM2Widget(SLIVR::CM2RectangleType type, float left_x, float left_y,
                     float width, float height, float offset);
  ~RectangleCM2Widget();
  virtual RectangleCM2Widget*    clone();

  virtual void          io(Piostream &stream);
  static PersistentTypeID type_id;
  virtual std::string dynamic_type_name() const { return type_id.type; }
};

class SCISHARE EllipsoidCM2Widget : 
  public CM2Widget, public SLIVR::EllipsoidCM2Widget
{
public:
  EllipsoidCM2Widget();
  EllipsoidCM2Widget(const EllipsoidCM2Widget &copy);
  EllipsoidCM2Widget(const SLIVR::EllipsoidCM2Widget *);
  EllipsoidCM2Widget(float x, float y, float a, float b, float rot);
  ~EllipsoidCM2Widget();
  virtual EllipsoidCM2Widget*    clone();

  virtual void          io(Piostream &stream);
  static PersistentTypeID type_id;
  virtual std::string dynamic_type_name() const { return type_id.type; }
};

class SCISHARE ParaboloidCM2Widget : 
  public CM2Widget, public SLIVR::ParaboloidCM2Widget
{
public:
  ParaboloidCM2Widget();
  ParaboloidCM2Widget(const ParaboloidCM2Widget &copy);
  ParaboloidCM2Widget(const SLIVR::ParaboloidCM2Widget *);
  ParaboloidCM2Widget(float top_x, float top_y,
		      float bottom_x, float bottom_y,
		      float left_x, float left_y,
		      float right_x, float right_y);
  ~ParaboloidCM2Widget();
  virtual ParaboloidCM2Widget*    clone();

  virtual void          io(Piostream &stream);
  static PersistentTypeID type_id;
  virtual std::string dynamic_type_name() const { return type_id.type; }
};



class SCISHARE ColorMapCM2Widget : 
  public CM2Widget, public SLIVR::ColorMapCM2Widget
{
public:
  ColorMapCM2Widget();
  ColorMapCM2Widget(const ColorMapCM2Widget &copy);
  ColorMapCM2Widget(const SLIVR::ColorMapCM2Widget *);
  ColorMapCM2Widget(SLIVR::CM2RectangleType type, float left_x, float left_y,
                     float width, float height, float offset);
  ~ColorMapCM2Widget();
  virtual ColorMapCM2Widget*    clone();

  virtual void          io(Piostream &stream);
  static PersistentTypeID type_id;
  virtual std::string dynamic_type_name() const { return type_id.type; }
};

// The image widget cannot be manipulated, only drawn.
class SCISHARE ImageCM2Widget : 
  public CM2Widget, public SLIVR::ImageCM2Widget
{
public:
  ImageCM2Widget();
  ImageCM2Widget(const ImageCM2Widget &copy);
  ImageCM2Widget(const SLIVR::ImageCM2Widget *);
  ImageCM2Widget(NrrdDataHandle p);
  ~ImageCM2Widget();

  virtual ImageCM2Widget* clone();

  virtual void          io(Piostream &stream);
  static PersistentTypeID type_id;
  virtual std::string dynamic_type_name() const { return type_id.type; }
};

class SCISHARE PaintCM2Widget :
  public CM2Widget, public SLIVR::PaintCM2Widget
{
public:
  PaintCM2Widget();
  ~PaintCM2Widget();
  PaintCM2Widget(const SLIVR::PaintCM2Widget *);
  virtual PaintCM2Widget*    clone();
  
  virtual void          io(Piostream &stream);
  static PersistentTypeID type_id;
  virtual std::string dynamic_type_name() const { return type_id.type; }
};

} // End namespace SCIRun

#endif // CM2Widget_h
