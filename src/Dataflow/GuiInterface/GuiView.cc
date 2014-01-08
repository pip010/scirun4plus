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
 *  GuiView.cc  Structure that provides for easy access of view information.
 *              The view information is interactively provided by the user.
 *
 *  Written by:
 *   Steven Parker
 *   Department of Computer Science
 *   University of Utah
 *
 *   separated from the Viewer code by me (Aleksandra)
 *   in May 1996
 *
 */

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Thread/Time.h>

#include <Dataflow/GuiInterface/GuiView.h>
#include <Dataflow/GuiInterface/GuiContext.h>
#include <Dataflow/GuiInterface/TCLInterface.h>

using namespace SCIRun;

#include <iostream>

GuiView::GuiView(GuiContext* ctx)
: GuiVar(ctx), eyep_(ctx->subVar("eyep")),
  lookat_(ctx->subVar("lookat")), up_(ctx->subVar("up")),
  fov_(ctx->subVar("fov")), eyep_offset_(ctx->subVar("eyep_offset"))
{
}

GuiView::~GuiView()
{
}

View
GuiView::get()
{
  View v(eyep_.get(), lookat_.get(), up_.get(), fov_.get());
  return v;
}

View
GuiView::get_cached()
{
  View v(eyep_.get_cached(), lookat_.get_cached(), up_.get_cached(), fov_.get_cached());
  return v;
}

void
GuiView::request()
{
  eyep_.request();
  lookat_.request();
  up_.request();
  fov_.request();
}

void
GuiView::set(const View& view)
{
  eyep_.set(view.eyep());
  lookat_.set(view.lookat());
  up_.set(view.up());
  fov_.set(view.fov());}

void
GuiView::send(const View& view)
{
  eyep_.send(view.eyep());
  lookat_.send(view.lookat());
  up_.send(view.up());
  fov_.send(view.fov());
}

GuiExtendedView::GuiExtendedView(GuiContext* ctx)
: GuiVar(ctx), eyep_(ctx->subVar("eyep")),
  lookat_(ctx->subVar("lookat")), up_(ctx->subVar("up")),
  fov_(ctx->subVar("fov")), eyep_offset_(ctx->subVar("eyep_offset")),
  xres_(ctx->subVar("xres")), yres_(ctx->subVar("yres")), bg_(ctx->subVar("bg"))
{
}

GuiExtendedView::~GuiExtendedView()
{
}

ExtendedView
GuiExtendedView::get()
{
  ExtendedView v(eyep_.get(), lookat_.get(), up_.get(), fov_.get(), xres_.get(),
		 yres_.get(), bg_.get()*( 1.0 / 255.0 ) );
  return v;
}

ExtendedView
GuiExtendedView::get_cached()
{
  ExtendedView v(eyep_.get_cached(), lookat_.get_cached(), up_.get_cached(), fov_.get_cached(), xres_.get_cached(),
		 yres_.get_cached(), bg_.get_cached()*( 1.0 / 255.0 ) );
  return v;
}

void
GuiExtendedView::request()
{
  eyep_.request();
  lookat_.request();
  up_.request();
  fov_.request();
  xres_.request();
  yres_.request();
  bg_.request();
}


void
GuiExtendedView::set(const ExtendedView& view)
{
  eyep_.set(view.eyep());
  lookat_.set(view.lookat());
  up_.set(view.up());
  fov_.set(view.fov());
  xres_.set(view.xres());
  yres_.set(view.yres());
  bg_.set( view.bg()*255 );}

void
GuiExtendedView::send(const ExtendedView& view)
{
  eyep_.send(view.eyep());
  lookat_.send(view.lookat());
  up_.send(view.up());
  fov_.send(view.fov());
  xres_.send(view.xres());
  yres_.send(view.yres());
  bg_.send( view.bg()*255 );
}
