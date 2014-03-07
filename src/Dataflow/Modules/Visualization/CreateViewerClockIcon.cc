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
 *  CreateViewerClockIcon.cc: Choose one input field to be passed downstream
 *
 *  Written by:
 *   Allen R. Sanderson
 *   SCI Institute
 *   University of Utah
 *   January 2004
 *
 */

#include <stdio.h>


#include <Dataflow/Network/Module.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Ports/NrrdPort.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Core/Geom/GeomGroup.h>
#include <Core/Geom/GeomBox.h>
#include <Core/Geom/GeomDisk.h>
#include <Core/Geom/GeomLine.h>
#include <Core/Geom/GeomText.h>
#include <Core/Geom/GeomTransform.h>
#include <Core/Geom/GeomMaterial.h>
#include <Core/Geom/GeomSticky.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

namespace SCIRun {

class CreateViewerClockIcon : public Module {
  public:
    CreateViewerClockIcon(GuiContext *context);
    virtual ~CreateViewerClockIcon() {}

    virtual void execute();

  protected:
    GeomHandle generateAnalog ();
    GeomHandle generateDigital();
    GeomHandle generateTime( int &nchars );

    virtual void tcl_command(GuiArgs& args, void* userdata);

  private:
    GeomHandle geometry_out_handle_;

    GuiInt    gui_type_;
    GuiInt    gui_bbox_;
    GuiString gui_format_;
    GuiDouble gui_min_;
    GuiDouble gui_max_;
    GuiString gui_current_;
    GuiInt    gui_size_;
    GuiDouble gui_location_x_;
    GuiDouble gui_location_y_;
    GuiDouble gui_color_r_;
    GuiDouble gui_color_g_;
    GuiDouble gui_color_b_;

    MaterialHandle material_handle_;

    bool color_changed_;
};


DECLARE_MAKER(CreateViewerClockIcon)


CreateViewerClockIcon::CreateViewerClockIcon(GuiContext *context)
  : Module("CreateViewerClockIcon", context, Source, "Visualization", "SCIRun"),
    geometry_out_handle_(0),
    gui_type_(context->subVar("type"), 0),
    gui_bbox_(context->subVar("bbox"), 1),
    gui_format_(context->subVar("format"), "%8.3f seconds"),
    gui_min_(context->subVar("min"), 0),
    gui_max_(context->subVar("max"), 1),
    gui_current_(context->subVar("current"), "0"),
    gui_size_(context->subVar("size"), 100),
    gui_location_x_(context->subVar("location-x"), -31.0/32.0),
    gui_location_y_(context->subVar("location-y"),  31.0/32.0),
    gui_color_r_(context->subVar("color-r"), 1.0),
    gui_color_g_(context->subVar("color-g"), 1.0),
    gui_color_b_(context->subVar("color-b"), 1.0),
    material_handle_(new Material(Color(1.0, 1.0, 1.0))),
    color_changed_(false)
{
}

void CreateViewerClockIcon::execute(){

  // Get the time via a matrix
  MatrixHandle input_matrix_handle;
  get_input_handle( "Time Matrix", input_matrix_handle, false );

  // Get the time via a nrrd
  NrrdDataHandle input_nrrd_handle;
  get_input_handle( "Time Nrrd", input_nrrd_handle, false );

  double time;

  if( input_matrix_handle.get_rep() ) 
  {
    if( input_matrix_handle->nrows() == 1 &&
      input_matrix_handle->ncols() == 1 )
    {
      time = input_matrix_handle->get(0, 0);
    }
    else 
    {
      error( "Input time matrix does not contain a single value." );
      return;
    }

  } 
  else if( input_nrrd_handle.get_rep() ) 
  {

    if( input_nrrd_handle->nrrd_->dim == 1 &&
	input_nrrd_handle->nrrd_->axis[0].size == 1 )
      time = get_nrrd_value( input_nrrd_handle->nrrd_, 0 );
    else {
      error( "Input time nrrd does not contain a single value." );
      return;
    }

  } else {
    error( "No input time present" );
    return;
  }

  if( gui_current_.get() != to_string( time ) ) {
    gui_current_.set( to_string( time ) ); 
    gui_current_.reset();
    
    inputs_changed_  = true;
  }

  if( inputs_changed_ ||
      !geometry_out_handle_.get_rep() ||
      gui_current_.changed( true ) ||
      gui_type_.changed( true ) ||
      gui_bbox_.changed( true ) ||
      gui_min_.changed( true ) ||
      gui_max_.changed( true ) ||
      gui_format_.changed( true ) ||
      gui_size_.changed( true ) ||
      gui_location_x_.changed( true ) ||
      gui_location_y_.changed( true ) ||
      color_changed_  ) 
  {
    // Inform module that execution started
    update_state(Executing);

    color_changed_  = false;

    material_handle_->diffuse =
      Color(gui_color_r_.get(), gui_color_g_.get(), gui_color_b_.get());

    if( gui_type_.get() % 2 == 0 )       // Analog or Analog/Digial
      geometry_out_handle_ = generateAnalog();

    else if( gui_type_.get() == 1 )      // Digital
      geometry_out_handle_ = generateDigital();
    
    send_output_handle( std::string("Clock"),
			geometry_out_handle_,
			std::string("Clock Sticky") );
  }
}

GeomHandle CreateViewerClockIcon::generateAnalog()
{
  Vector refVec = 31.0/32.0 *
    Vector( gui_location_x_.get(), gui_location_y_.get(), 0.0 );

  GeomGroup *group = new GeomGroup();

  double gsizes[5] = { 0.25, 0.50, 1.0, 1.5, 2.0 };  
  double radius = 0.075 * gsizes[gui_size_.get()];

  double dx = 0, dy = 0;

  if( gui_type_.get() == 2 ) { 
    int nchars = 0;

    group->add( generateTime( nchars ) );

    double fsizes[5] = { 60,100,140,180,240 };
    double border = 0.0125;
    double scale = .000225;
    
    dx = nchars * fsizes[gui_size_.get()] * scale;
    dy =    1.8 * fsizes[gui_size_.get()] * scale + border;
    
    refVec += Vector( dx/2-radius, dy, 0 );
  }



  GeomGroup *circle = new GeomGroup();


  refVec += Vector( radius, radius, 0 );

  Point last = refVec + Point( 0, radius, 0 );

  for (int i=1; i<=36; i++) {
    double cx = sin( 2.0 * M_PI * i / 36.0 ) * radius;
    double cy = cos( 2.0 * M_PI * i / 36.0 ) * radius;
    
    Point next = refVec + Point( cx, cy, 0 );

    circle->add( new GeomLine( last, next ) );

    last = next;
  }

  double value;
  string_to_double(gui_current_.get(), value );
 
  double t =
    2.0 * M_PI * (value-gui_min_.get()) / (gui_max_.get()-gui_min_.get());
  double cx = sin( t ) * radius;
  double cy = cos( t ) * radius;

  circle->add( new GeomLine( refVec + Point(0,0,0),
			     refVec + Point( cx, cy, 0 ) ) );
  
  group->add( new GeomMaterial(circle, material_handle_) );
  
  double border = 0.0125;

  if( dx < 2.0 * radius + border )
    dx = 2.0 * radius + border;

  dy += 2.0 * radius + border;

  if( gui_bbox_.get() ) {
    Vector refVec = 31.0/32.0 *
      Vector( gui_location_x_.get(), gui_location_y_.get(), 0.0 );

    GeomGroup *box = new GeomGroup();
    
    box->add( new GeomLine( refVec+Point(-border,-border,0),
			    refVec+Point(dx,-border,0) ) );
    
    box->add( new GeomLine( refVec+Point(dx,-border,0),
			    refVec+Point(dx,dy,0) ) );

    box->add( new GeomLine( refVec+Point(dx,dy,0),
			    refVec+Point(-border,dy,0) ) );
    
    box->add( new GeomLine( refVec+Point(-border,dy,0),
			    refVec+Point(-border,-border,0) ) );
    
    group->add( new GeomMaterial(box, material_handle_) );
  }

  return new GeomSticky( group );
}

GeomHandle CreateViewerClockIcon::generateDigital()
{
  GeomGroup *group = new GeomGroup();

  int nchars = 0;

  group->add( generateTime( nchars ) );

  if( gui_bbox_.get() ) {
    
    double fsizes[5] = { 60,100,140,180,240 };
    double border = 0.0125;
    double scale = .000225;

    double dx = nchars * fsizes[gui_size_.get()] * scale + border;
    double dy =    1.8 * fsizes[gui_size_.get()] * scale + border;

    Vector refVec = 31.0/32.0 *
      Vector( gui_location_x_.get(), gui_location_y_.get(), 0.0 );

    GeomGroup *box = new GeomGroup();
    
    box->add( new GeomLine( refVec+Point(-border,-border,0),
			    refVec+Point(dx,-border,0) ) );
    
    box->add( new GeomLine( refVec+Point(dx,-border,0),
			    refVec+Point(dx,dy,0) ) );

    box->add( new GeomLine( refVec+Point(dx,dy,0),
			    refVec+Point(-border,dy,0) ) );
    
    box->add( new GeomLine( refVec+Point(-border,dy,0),
			    refVec+Point(-border,-border,0) ) );
    
    group->add( new GeomMaterial(box, material_handle_) );
  }

  return new GeomSticky( group );
}


GeomHandle CreateViewerClockIcon::generateTime( int &nchars )
{
  std::string format = gui_format_.get();

  char timestr[64];

  if( format.find("%") == std::string::npos ||
      format.find("%") != format.find_last_of("%") ) {
    error("Bad C Style format for the clock.");
    error("The format should be of the form: '%7.4f seconds'");
    sprintf( timestr, "Bad Format" );
  } else {
    double value;
    string_to_double(gui_current_.get(), value );
    sprintf( timestr, format.c_str(), value );
  }

  std::string istr(timestr);

  nchars = istr.size();
    
  GeomTexts* texts = new GeomTexts();
  texts->set_font_index(gui_size_.get());
  texts->set_is_2d(true);

  Vector refVec = 31.0/32.0 *
    Vector( gui_location_x_.get(), gui_location_y_.get(), 0.0 );

  Point loc = (Point) refVec;

  texts->add(istr, loc, material_handle_->diffuse);
  
  GeomHandle handle = dynamic_cast<GeomObj *>(texts);

  return (handle);
}


void 
CreateViewerClockIcon::tcl_command(GuiArgs& args, void* userdata) {

  if(args.count() < 2) {
    args.error("CreateViewerClockIcon needs a minor command");
    return;
  }

  if (args[1] == "color_change") {
    color_changed_ = true;

    // The below works for only the geometry not for the text so update.
    /*
    // Get a handle to the output geom port.
    GeometryOPort *ogeom_port = (GeometryOPort *) get_oport("Clock");
    
    gui_color_r_.reset();
    gui_color_g_.reset();
    gui_color_b_.reset();

    material_->diffuse = Color(gui_color_r_.get(),
                               gui_color_g_.get(),
                               gui_color_b_.get());
    
    if (ogeom_port) ogeom_port->flushViews();
    */
  } else {
    Module::tcl_command(args, userdata);
  }
}
} // End namespace SCIRun

