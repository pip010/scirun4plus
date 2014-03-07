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
 *  GuiVar.cc: Interface to TCL variables
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   September 1994
 *
 */

#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/GuiInterface/TCLInterface.h>

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Util/Assert.h>

using namespace SCIRun;

#include <iostream>

GuiVar::GuiVar(GuiContext* ctx)
  : ref_cnt(0), ctx(ctx)
{
  ASSERT(ctx);
  ctx->doReplaceEnv();
}

GuiVar::~GuiVar()
{
  if (ctx) delete ctx;
}

void 
GuiVar::reset()
{
  ASSERT(ctx);
  ctx->reset();
}

void 
GuiVar::synchronize()
{
  TCLInterface::synchronize();
}

void GuiVar::disable_environment_replace()
{
  ctx->dontReplaceEnv();
}

void GuiVar::enable_environment_replace()
{
  ctx->doReplaceEnv();
}

template class GuiSingle<std::string>;
template class GuiSingle<double>;
template class GuiSingle<int>;
template class GuiTriple<Point>;
template class GuiTriple<Vector>;
