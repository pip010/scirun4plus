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
 *  GuiGeom.cc: Interface to TCL variables for Geom stuff
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   June 1995
 *
 */

#include <Dataflow/GuiInterface/GuiGeom.h>
#include <Dataflow/GuiInterface/GuiContext.h>
#include <Dataflow/GuiInterface/TCLInterface.h>
#include <Core/Datatypes/Color.h>
#include <Core/Datatypes/Material.h>

using namespace SCIRun;

GuiColor::GuiColor(GuiContext* ctx)
  : GuiVar(ctx), r_(ctx->subVar("r")), g_(ctx->subVar("g")),
    b_(ctx->subVar("b"))
{
}

GuiColor::GuiColor(GuiContext* ctx,Color c)
  : GuiVar(ctx), r_(ctx->subVar("r")), g_(ctx->subVar("g")),
    b_(ctx->subVar("b"))
{
  r_.set(c.r());
  g_.set(c.g());
  b_.set(c.b());
}


GuiColor::~GuiColor()
{
}

Color 
GuiColor::get()
{
  Color c(r_.get(), g_.get(), b_.get());
  return c;
}

Color 
GuiColor::get_cached()
{
  Color c(r_.get_cached(), g_.get_cached(), b_.get_cached());
  return c;
}

void
GuiColor::request()
{
  r_.request();
  g_.request();
  b_.request();
} 

void 
GuiColor::set(const Color& p)
{
  r_.set(p.r());
  g_.set(p.g());
  b_.set(p.b());}

void 
GuiColor::send(const Color& p)
{
  r_.send(p.r());
  g_.send(p.g());
  b_.send(p.b());
}


GuiMaterial::GuiMaterial(GuiContext* ctx)
: GuiVar(ctx), ambient_(ctx->subVar("ambient")),
  diffuse_(ctx->subVar("diffuse")), 
  specular_(ctx->subVar("specular")),
  shininess_(ctx->subVar("shininess")), 
  emission_(ctx->subVar("emission")),
  reflectivity_(ctx->subVar("reflectivity")),
  transparency_(ctx->subVar("transparency")),
  refraction_index_(ctx->subVar("refraction_index"))
{
}

GuiMaterial::~GuiMaterial()
{
}

Material 
GuiMaterial::get()
{
  Material m(ambient_.get(), diffuse_.get(), specular_.get(), shininess_.get());
  m.emission=emission_.get();
  m.reflectivity=reflectivity_.get();
  m.transparency=transparency_.get();
  m.refraction_index=refraction_index_.get();

  return m;
}

Material 
GuiMaterial::get_cached()
{
  Material m(ambient_.get_cached(), diffuse_.get_cached(), specular_.get_cached(), shininess_.get_cached());
  m.emission=emission_.get_cached();
  m.reflectivity=reflectivity_.get_cached();
  m.transparency=transparency_.get_cached();
  m.refraction_index=refraction_index_.get_cached();

  return m;
}

void
GuiMaterial::request()
{
  ambient_.request();
  diffuse_.request();
  specular_.request();
  shininess_.request();
  emission_.request();
  transparency_.request();
  refraction_index_.request();
}

void 
GuiMaterial::set(const Material& m)
{
  ambient_.set(m.ambient);
  diffuse_.set(m.diffuse);
  specular_.set(m.specular);
  shininess_.set(m.shininess);
  emission_.set(m.emission);
  reflectivity_.set(m.reflectivity);
  transparency_.set(m.transparency);
  refraction_index_.set(m.refraction_index);}

void 
GuiMaterial::send(const Material& m)
{
  ambient_.send(m.ambient);
  diffuse_.send(m.diffuse);
  specular_.send(m.specular);
  shininess_.send(m.shininess);
  emission_.send(m.emission);
  reflectivity_.send(m.reflectivity);
  transparency_.send(m.transparency);
  refraction_index_.send(m.refraction_index);
}


