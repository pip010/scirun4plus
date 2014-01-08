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
#include <Dataflow/Network/Ports/ColorMapPort.h>
#include <Dataflow/Network/Ports/FieldPort.h>
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

class ShowColorMap : public Module {
  public:
    ShowColorMap(GuiContext*);
    virtual ~ShowColorMap();
    virtual void execute();
    virtual void tcl_command(GuiArgs& args, void* userdata);

  private:
    GeomHandle geometry_output_handle_;

    GuiString gui_length_;
    GuiString gui_side_;
    GuiInt    gui_numlabels_;
    GuiDouble gui_scale_;
    GuiInt    gui_num_sig_digits_;
    GuiString gui_units_;
    GuiDouble gui_color_r_;
    GuiDouble gui_color_g_;
    GuiDouble gui_color_b_;
    GuiInt    gui_text_fontsize_;
    GuiInt    gui_extra_padding_;

    MaterialHandle text_material_handle_;

    bool color_changed_;
};

  DECLARE_MAKER(ShowColorMap)

ShowColorMap::ShowColorMap(GuiContext* context)
  : Module("ShowColorMap", context, Filter, "Visualization", "SCIRun"),
    geometry_output_handle_(0),
    gui_length_(context->subVar("length"), "half2"),
    gui_side_(context->subVar("side"), "left"),
    gui_numlabels_(context->subVar("numlabels"), 5),
    gui_scale_(context->subVar("scale"), 1.0),
    gui_num_sig_digits_(context->subVar("numsigdigits"), 2),
    gui_units_(context->subVar("units"), ""),
    gui_color_r_(context->subVar("color-r"), 1.0),
    gui_color_g_(context->subVar("color-g"), 1.0),
    gui_color_b_(context->subVar("color-b"), 1.0),
    gui_text_fontsize_(context->subVar("text-fontsize"), 2),
    gui_extra_padding_(context->subVar("extra-padding"), 0),
    text_material_handle_(new Material(Color(1.0, 1.0, 1.0))),
    color_changed_(false)
{
}


ShowColorMap::~ShowColorMap()
{
}


void
ShowColorMap::execute()
{
  ColorMapHandle cmHandle = 0;
  if( !get_input_handle( "ColorMap", cmHandle, true ) ) return;

  if( inputs_changed_ ||
      !geometry_output_handle_.get_rep() ||
      gui_length_.changed( true ) ||
      gui_side_.changed( true ) ||
      gui_numlabels_.changed( true ) ||
      gui_scale_.changed( true ) ||
      gui_num_sig_digits_.changed( true ) ||
      gui_units_.changed( true ) ||
      gui_text_fontsize_.changed( true ) ||
      gui_extra_padding_.changed( true ) ||
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

    if (gui_side_.get() == "left")
    {
      out = Vector(-0.05, 0.0, 0.0);
      if (gui_length_.get() == "full" || gui_length_.get() == "half2") 
      {
        if (gui_extra_padding_.get()) 
        {
          ref1 = Point(-1.0, -14.0/16.0, 0.0);
        } 
        else 
        {
          ref1 = Point(-1.0, -15.0/16.0, 0.0);
        }
      } 
      else
      {
        ref1 = Point(-1.0, 1.0/16.0, 1.0);
      }
	
      if (gui_length_.get() == "full") 
      {
        if (gui_extra_padding_.get()) 
        {
          along = Vector(0.0, 29.0/16.0, 0.0);
        } 
        else 
        {
          along = Vector(0.0, 30.0/16.0, 0.0);
        }
      } 
      else 
      {
        along = Vector(0.0, 14.0/16.0, 0.0);
      }
    } 
    else if (gui_side_.get() == "bottom") 
    {
      out = Vector(0.0, -0.05, 0.0);
      if (gui_length_.get() == "full" || gui_length_.get() == "half1") 
      {
        ref1 = Point(-15.0/16.0, -1.0, 0.0);
      } 
      else 
      {
        ref1 = Point(1.0/16.0, -1.0, 0.0);
      }

      if (gui_length_.get() == "full") 
      {
        if (gui_extra_padding_.get()) 
        {
          along = Vector(29.0/16.0, 0.0, 0.0);
        } 
        else 
        {
          along = Vector(30.0/16.0, 0.0, 0.0);
        }
      } 
      else 
      {
        if (gui_extra_padding_.get()) 
        {
          along = Vector(13.0/16.0, 0.0, 0.0);
        } 
        else 
        {
          along = Vector(14.0/16.0, 0.0, 0.0);
        }
      }
    }
      
    const Point  ref0(ref1 - out);
    // Create a new colormap that we can send to ColorMapTex.  We need
    // to do this, because if the colormap min/max are too close you get
    // problems.  This is because the min and max are used to lookup
    // into the texture as floats.  Precion problems occur when the min
    // == max in floats, but not as doubles.
    ColorMapHandle cmap_rescaled;
    // Only rescale it when the min != max or min and max are too close.
    float too_close = Abs((float)(cmHandle->getMin()) - (float)(cmHandle->getMax()));
    // Replace zero compared with too_close with an epsilon if desired.
    if (too_close <= 0) 
    {
      // Make a copy of the colormap we can rescale
      cmap_rescaled = cmHandle->clone();
      cmap_rescaled->Scale(0, 1);
    } 
    else 
    {
      cmap_rescaled = cmHandle;
    }
    
    ColorMapTex *sq = new ColorMapTex(ref0, ref0 + along,
                              ref0 + along + out, ref0 + out, cmap_rescaled);
    all->add( sq );

    double scale = gui_scale_.get();
    std::string str = gui_units_.get();
    if (str.empty()) str = cmHandle->units();
    // So if the maximum number of digits the number will take up is
    // at most 25 then the length of str better be less than 80-25-1.
    // See size of value and num_sig_digits below.
    if (str.length() > 50) 
    {
      error("Length of units string is too long.  Make it smaller than 50 characters please.");
      return;
    }

    const int numlabels = gui_numlabels_.get();
    if (numlabels > 1 && numlabels < 50)
    {
      // Fill in the text.
      const double minval = cmHandle->getMin()*scale;
      const double maxval = cmHandle->getMax()*scale;

      Point p0  = ref0 - out * 0.02; 
      char value[80];
      GeomLines *lines = new GeomLines();
      GeomTexts *texts = new GeomTexts();
      texts->set_is_2d(true);
      texts->set_font_index(gui_text_fontsize_.get());
      int num_sig_digits = gui_num_sig_digits_.get();
      if (num_sig_digits < 1) 
      {
        warning("Number of significant digits needs to be at least 1.  Setting the number of significant digits to 1.");
        gui_num_sig_digits_.set(1);
        num_sig_digits = 1;
      }
      if (num_sig_digits > 20) 
      {
        warning("Number of significant digits needs to be less than or equal to 20.  Setting the number of significant digits to 20");
        gui_num_sig_digits_.set(20);
        num_sig_digits = 20;
      }
      for(int i = 0; i < numlabels; i++ )
      {
        sprintf(value, "%.*g %s", num_sig_digits,
          minval + (maxval-minval)*(i/(numlabels-1.0)),
          str.c_str());
        const Point loc = p0 + along * (i/(numlabels-1.0));
        texts->add(value, loc, text_material_handle_->diffuse);
        lines->add(loc, text_material_handle_,
                     loc + out * 0.5, text_material_handle_);
      }

      all->add(texts);
      all->add(lines);
    }

    GeomSticky *sticky = new GeomSticky(all);
    geometry_output_handle_ = GeomHandle( sticky );
    
    send_output_handle( "Geometry", geometry_output_handle_, "ShowColorMap Sticky" );
  }
}

void ShowColorMap::tcl_command(GuiArgs& args, void* userdata)
{
  if(args.count() < 2) 
  {
    args.error("ShowColorMap needs a minor command");
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
