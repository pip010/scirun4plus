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
 *  ColorMapKey.cc: create a key for colormap
 *
 *  Written by:
 *   Philip Sutton
 *   Department of Computer Science
 *   University of Utah
 *   May 1998
 *
 *  Updated by:
 *   Michael Callahan
 *   January 2001
 *
 */

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Core/Geom/GeomGroup.h>
#include <Core/Geom/GeomText.h>
#include <Core/Geom/GeomLine.h>
#include <Core/Geom/ColorMapTex.h>
#include <Core/Geom/GeomTransform.h>
#include <Core/Geometry/Transform.h>

#include <Core/Geom/GeomSticky.h>

#include <stdio.h>
#include <string.h>
#include <iostream>

namespace SCIRun {

class CreateScaleBar : public Module {
public:
  CreateScaleBar(GuiContext*);
  virtual ~CreateScaleBar();
  virtual void execute();
  virtual void tcl_command(GuiArgs& args, void* userdata);

private:
  GeomHandle geometry_output_handle_;

  GuiDouble gui_scale_;
  GuiString gui_label_;
  GuiDouble gui_color_r_;
  GuiDouble gui_color_g_;
  GuiDouble gui_color_b_;

  MaterialHandle text_material_handle_;

  bool color_changed_;
};

  DECLARE_MAKER(CreateScaleBar)

CreateScaleBar::CreateScaleBar(GuiContext* context)
  : Module("CreateScaleBar", context, Filter, "Visualization", "SCIRun"),
    geometry_output_handle_(0),
    gui_scale_(context->subVar("scale"), 1.0),
    gui_label_(context->subVar("label"), "1mm"),
    gui_color_r_(context->subVar("color-r"), 1.0),
    gui_color_g_(context->subVar("color-g"), 1.0),
    gui_color_b_(context->subVar("color-b"), 1.0),
    text_material_handle_(new Material(Color(1.0, 1.0, 1.0))),
    color_changed_(false)
{
}


CreateScaleBar::~CreateScaleBar()
{
}


void
CreateScaleBar::execute()
{
  if( !geometry_output_handle_.get_rep() ||
      gui_scale_.changed( true ) ||
      gui_label_.changed( true ) ||
      color_changed_ == true )
  {
    // Inform module that execution started
    update_state(Executing);

    GeomGroup *all = new GeomGroup();

    Point  ref1(0.0, 0.0, 0.0);
    Vector out(0.0, 0.0, 0.0);
    Vector along(0.0, 0.0, 0.0);

    color_changed_ = false;

    text_material_handle_->diffuse =
      Color(gui_color_r_.get(), gui_color_g_.get(), gui_color_b_.get());

    out = Vector(0.0, -0.12, 0.0);
    
    double scale = gui_scale_.get() * 7.0/10.0;
    ref1 = Point(15.0/16.0 - scale/2.0, -1.0, 0.0);
    along = Vector(scale/2.0 - 1.0/16.0, 0.0, 0.0);
    
    const Point  ref0(ref1 - out);
    std::string label = gui_label_.get();
    if (label.length() > 50) {
      error("Length of label string is too long.  Make it smaller than 50 characters please.");
      return;
    }

    const int numTicks = 11;

    // Fill in the text.

    scale *= 10./7.0;
    Point p0  = ref0 - out * 0.02; 
    GeomLines *lines = new GeomLines();
    GeomTexts *texts = new GeomTexts();
    texts->set_is_2d(true);
    texts->set_font_index(0);
    lines->add(p0+out*0.5, text_material_handle_, 
	       p0+along+out*0.5, text_material_handle_);
    for(int i = 0; i < numTicks; i++ )
      {
 	const Point loc = p0 + along * (i/(numTicks-1.0));
	if (i==numTicks-1) texts->add(label, loc-Vector(0.04,0,0)*((label.length()-1.0)/2.0),
			     text_material_handle_->diffuse);
	if (i==0 || i==numTicks-1) {
	  lines->add(loc,             text_material_handle_,
		     loc + out * 0.5, text_material_handle_);
	} else {
	  lines->add(loc + out * 0.2, text_material_handle_,
		     loc + out * 0.5, text_material_handle_);
	}
      }
    
    all->add(texts);
    all->add(lines);

    GeomSticky *sticky = new GeomSticky(all);
    geometry_output_handle_ = GeomHandle( sticky );
    send_output_handle( "Geometry",
			geometry_output_handle_,
			"CreateScaleBar Sticky" );
  }
}

void CreateScaleBar::tcl_command(GuiArgs& args, void* userdata)
{
  if(args.count() < 2) 
  {
    args.error("CreateScaleBar needs a minor command");
    return;
  }

  if (args[1] == "color_change") 
  {
    color_changed_ = true;
  } 
  else 
  {
    Module::tcl_command(args, userdata);
  }
}

} // End namespace SCIRun
