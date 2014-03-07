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
//    File   : ColorMap2.cc
//    Author : Milan Ikits
//    Date   : Mon Jul  5 18:33:29 2004

#include <Core/Persistent/Persistent.h>
#include <Core/Persistent/PersistentSTL.h>
#include <Core/Volume/ColorMap2.h>
#include <Core/Volume/CM2Widget.h>

#include <Core/Util/Debug.h>


namespace SCIRun {

static Persistent* maker()
{
  return new ColorMap2;
}

PersistentTypeID ColorMap2::type_id("ColorMap2", "PropertyManager", maker);

#define COLORMAP2_VERSION 4

void
ColorMap2::io(Piostream &stream)
{
  const int version = stream.begin_class("ColorMap2", COLORMAP2_VERSION);

  if (version >= 2)
    PropertyManager::io(stream);

  bool faux = false;
  if (version <= 3)
    SCIRun::Pio(stream, faux);

  std::vector<CM2WidgetHandle> tmp;
  if (stream.writing()) {
    // load tmp with handles if we are writing.
    for (size_t i = 0; i  < widgets_.size(); ++i)
    {
      CM2WidgetHandle h(CM2Widget::slivr2sr(widgets_[i]));
      tmp.push_back(h);
    }
  }

  SCIRun::Pio(stream, tmp);
    
  if (stream.reading()) {
    std::vector<CM2WidgetHandle>::iterator iter = tmp.begin();
    while(iter != tmp.end()) {
      CM2WidgetHandle h = *iter;
      widgets_.push_back(h->clone());
      ++iter;
    }
  }

  if (version <= 3)
    for (unsigned int w = 0; w < widgets_.size(); ++w)
      widgets_[w]->set_faux(faux);
      
  if (version >= 3)
    SCIRun::Pio(stream, selected_);

  if (version >= 4)
    SCIRun::Pio(stream, value_range_);

  stream.end_class();
}

ColorMap2::ColorMap2() :
  SLIVR::ColorMap2()
{
  DEBUG_CONSTRUCTOR("ColorMap2")
}

ColorMap2::ColorMap2(const ColorMap2 &copy) :
  SLIVR::ColorMap2((SLIVR::ColorMap2&)copy),
  PropertyManager(copy)
{
  DEBUG_CONSTRUCTOR("ColorMap2")
}

ColorMap2::ColorMap2(const std::vector<SLIVR::CM2Widget*>& widgets,
                     bool updating,
                     bool selected,
                     std::pair<float, float> value_range) :
  SLIVR::ColorMap2(widgets, updating, selected, value_range)
{
  DEBUG_CONSTRUCTOR("ColorMap2")
}

ColorMap2::~ColorMap2()
{
  DEBUG_DESTRUCTOR("ColorMap2")
}

} // End namespace SCIRun
