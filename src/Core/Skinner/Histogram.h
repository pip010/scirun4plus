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
//
//    File   : Histogram.h
//    Author : McKay Davis
//    Date   : Tue Jun 27 13:00:11 2006

#ifndef SKINNER_HISTOGRAM_H
#define SKINNER_HISTOGRAM_H

#include <Core/Skinner/Parent.h>
#include <Core/Skinner/Color.h>
#include <Core/Skinner/Variables.h>
#include <Core/Datatypes/NrrdData.h>

#include <Core/Skinner/share.h>

namespace SCIRun {
class TextureObj;
namespace Skinner {


class SCISHARE Histogram : public Parent {
public:
  Histogram(Variables *);
  virtual ~Histogram();
  CatcherFunction_t process_event;
  static NrrdDataHandle    UnuQuantize(NrrdDataHandle,
				       double minf=0.0,
				       double maxf=0.0,
				       unsigned int nbits=8);
  static NrrdDataHandle    GenerateHistogram(NrrdDataHandle);
private:
  static NrrdDataHandle    UnuGamma(NrrdDataHandle, double);
  static NrrdDataHandle    UnuHeq(NrrdDataHandle);
  static NrrdDataHandle    Unu1op(NrrdDataHandle);
  static NrrdDataHandle    Unu2op(NrrdDataHandle);
  static NrrdDataHandle    UnuJhisto(NrrdDataHandle, NrrdDataHandle);

  static NrrdDataHandle BuildJHistoDirectly(NrrdDataHandle in);

  CatcherFunction_t do_PointerEvent;
  CatcherFunction_t redraw;

  NrrdDataHandle    histo_nrrd_;
  TextureObj *      tex_;
  Var<string>       filename_;
  Var<double>       opacity_;
  Var<double>       gamma_;
  double            old_gamma_;
  unsigned int      current_time_;
};


}
}

#endif
