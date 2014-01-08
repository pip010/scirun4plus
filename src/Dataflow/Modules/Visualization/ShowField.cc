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
 *  ShowField.cc
 *
 *  Written by:
 *   Allen R. Sanderson
 *   Marty Cole
 *   School of Computing
 *   University of Utah
 *   March 15, 2007
 *
 */ 

#include <Core/Util/StringUtil.h>
#include <Core/Datatypes/Field.h>
#include <Core/Geom/GeomMaterial.h>
#include <Core/Geom/GeomSwitch.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <Core/Algorithms/Visualization/RenderField.h>
#include <Core/Algorithms/Util/FieldInformation.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/GeometryPort.h>

#include <typeinfo>
#include <iostream>

namespace SCIRun {

class ShowField : public Module, public RenderStateBase
{
  public:
    ShowField(GuiContext* ctx);
    virtual ~ShowField();

    virtual void execute();
    virtual void tcl_command(GuiArgs& args, void* userdata);

    virtual void post_read();

  protected:
    //! input ports
    int                      mesh_generation_;
    bool                     color_map_present_;

    //! Scene graph ID's
    int                      nodes_id_;
    int                      edges_id_;
    int                      faces_id_;
    int                      text_id_;

    //! top level nodes for switching on and off..
    //! Options for rendering nodes.
    GuiInt                   nodes_on_;
    GuiInt                   nodes_transparency_;
    GuiInt                   nodes_color_type_;
    GuiString                nodes_display_type_;

    //! Options for rendering edges.
    GuiInt                   edges_on_;
    GuiInt                   edges_transparency_;
    GuiInt                   edges_color_type_;
    GuiString                edges_display_type_;

    //! Options for rendering faces.
    GuiInt                   faces_on_;
    GuiInt                   faces_transparency_;
    GuiInt                   faces_color_type_;
    GuiInt                   faces_normals_;
    GuiInt                   faces_usetexture_;

    //! Options for rendering text.
    GuiInt                   text_on_;
    GuiInt                   text_color_type_;
    GuiDouble                text_color_r_;
    GuiDouble                text_color_g_;
    GuiDouble                text_color_b_;
    GuiInt                   text_backface_cull_;
    GuiInt                   text_always_visible_;
    GuiInt                   text_fontsize_;
    GuiInt                   text_precision_;
    GuiInt                   text_render_locations_;
    GuiInt                   text_show_data_;
    GuiInt                   text_show_nodes_;
    GuiInt                   text_show_edges_;
    GuiInt                   text_show_faces_;
    GuiInt                   text_show_cells_;
    MaterialHandle           text_material_;

    std::string                   field_data_basis_type_;

    //! default color and material
    GuiDouble                def_color_r_;
    GuiDouble                def_color_g_;
    GuiDouble                def_color_b_;
    GuiDouble                def_color_a_;
    MaterialHandle           def_material_;

    GuiDouble                nodes_scale_;
    GuiDouble                nodes_scaleNV_;
    GuiDouble                edges_scale_;
    GuiDouble                edges_scaleNV_;

    GuiString                active_tab_; //! for saving nets state
    GuiString                interactive_mode_;
    GuiInt                   showProgress_;

    GuiString                gui_field_name_;
    GuiInt                   gui_field_name_override_;

    //! Refinement resolution for cylinders and spheres
    GuiInt                   gui_nodes_resolution_;
    GuiInt                   gui_edges_resolution_;

    GuiInt                   approx_div_;

    LockingHandle<RenderFieldBase>  renderer_;

    // variables related to default scale factor usage.
    GuiInt                  gui_use_default_size_;
    double                  cur_mesh_scale_factor_;

    unsigned int render_state_[4];

    enum toggle_type_e {
      NODE = 0,
      EDGE = 1,
      FACE = 2,
      TEXT = 3,
      ALL = 4
    };

    void maybe_execute(toggle_type_e dis_type);
    void set_default_display_values();
};

ShowField::ShowField(GuiContext* ctx) :
  Module("ShowField", ctx, Filter, "Visualization", "SCIRun","2.0"),
  mesh_generation_(-1),
  color_map_present_(0),
  nodes_id_(0),
  edges_id_(0),
  faces_id_(0),
  text_id_(0),

  nodes_on_(get_ctx()->subVar("nodes_on"), 1),
  nodes_transparency_(get_ctx()->subVar("nodes_transparency"), 0),
  nodes_color_type_(get_ctx()->subVar("nodes_color_type"), 1),
  nodes_display_type_(get_ctx()->subVar("nodes_display_type"), "Points"),

  edges_on_(get_ctx()->subVar("edges_on"), 1),
  edges_transparency_(get_ctx()->subVar("edges_transparency"), 0),
  edges_color_type_(get_ctx()->subVar("edges_color_type"), 1),
  edges_display_type_(get_ctx()->subVar("edges_display_type"), "Lines"),

  faces_on_(get_ctx()->subVar("faces_on"), 1),
  faces_transparency_(get_ctx()->subVar("faces_transparency"), 0),
  faces_color_type_(get_ctx()->subVar("faces_color_type"), 1),
  faces_normals_(get_ctx()->subVar("faces_normals"), 0),
  faces_usetexture_(get_ctx()->subVar("faces_usetexture"), 0),

  text_on_(get_ctx()->subVar("text_on"), 0),
  text_color_type_(get_ctx()->subVar("text_color_type"), 0),
  text_color_r_(get_ctx()->subVar("text_color-r"), 1.0),
  text_color_g_(get_ctx()->subVar("text_color-g"), 1.0),
  text_color_b_(get_ctx()->subVar("text_color-b"), 1.0),
  text_backface_cull_(get_ctx()->subVar("text_backface_cull"), 0),
  text_always_visible_(get_ctx()->subVar("text_always_visible"), 0),
  text_fontsize_(get_ctx()->subVar("text_fontsize"), 1),
  text_precision_(get_ctx()->subVar("text_precision"), 3),
  text_render_locations_(get_ctx()->subVar("text_render_locations"), 0),
  text_show_data_(get_ctx()->subVar("text_show_data"), 1),
  text_show_nodes_(get_ctx()->subVar("text_show_nodes"), 0),
  text_show_edges_(get_ctx()->subVar("text_show_edges"), 0),
  text_show_faces_(get_ctx()->subVar("text_show_faces"), 0),
  text_show_cells_(get_ctx()->subVar("text_show_cells"), 0),
  text_material_(new Material(Color(0.75, 0.75, 0.75))),

  field_data_basis_type_("none"),

  def_color_r_(get_ctx()->subVar("def_color-r"), 0.5),
  def_color_g_(get_ctx()->subVar("def_color-g"), 0.5),
  def_color_b_(get_ctx()->subVar("def_color-b"), 0.5),
  def_color_a_(get_ctx()->subVar("def_color-a"), 0.5),

  def_material_(0),

  nodes_scale_(get_ctx()->subVar("nodes_scale"), 0.03),
  nodes_scaleNV_(get_ctx()->subVar("nodes_scaleNV"), 0.03),
  edges_scale_(get_ctx()->subVar("edges_scale"), 0.15),
  edges_scaleNV_(get_ctx()->subVar("edges_scaleNV"), 0.15),

  active_tab_(get_ctx()->subVar("active_tab"), "Nodes"),
  interactive_mode_(get_ctx()->subVar("interactive_mode"), "Interactive"),
  showProgress_(get_ctx()->subVar("show_progress"), 0),

  gui_field_name_(get_ctx()->subVar("field_name"), ""),
  gui_field_name_override_(get_ctx()->subVar("field_name_override"), 0),

  gui_nodes_resolution_(get_ctx()->subVar("nodes_resolution"), 6),
  gui_edges_resolution_(get_ctx()->subVar("edges_resolution"), 6),

  approx_div_(get_ctx()->subVar("approx_div"), 1),

  renderer_(0),

  gui_use_default_size_(get_ctx()->subVar("use_default_size"), 0),
  cur_mesh_scale_factor_(1.0)
{
  Color gray(0.5, 0.5, 0.5);
  def_material_ = new Material(gray);
  
  render_state_[NODE] = 0;
  render_state_[EDGE] = 0;
  render_state_[FACE] = 0;
  render_state_[TEXT] = 0;
}


ShowField::~ShowField()
{
}


void
ShowField::execute()
{
  if( nodes_on_.get() == 0 &&
      edges_on_.get() == 0 &&
      faces_on_.get() == 0 &&
       text_on_.get() == 0 )
    return;

  bool update_algo = false;

  GeometryOPortHandle ogeom;
  get_oport_handle("Scene Graph",ogeom);

  FieldHandle fld_handle;

  if( !get_input_handle("Mesh", fld_handle, false ))
  {
    // Delete everything on output if no new input is given
    ogeom->delAll();
    return;
  }
  
  // Update the field name but only if the user does not enter an
  // overriding name.
  if (!gui_field_name_override_.get())
  {
    std::string fname("");

    if ( !fld_handle->get_property("name", fname))
      fld_handle->mesh()->get_property("name", fname);
       
    if( fname != gui_field_name_.get() )
    {
      gui_field_name_.set(fname);
      gui_field_name_.reset();
    }
  }

  // See if the algorithm needs to be recreated due the field or mesh
  // changing. This check will be made if colormap changes also.
  if( inputs_changed_ ) 
  {
    // Set inital colors here.
    def_material_->diffuse =
      Color(def_color_r_.get(), def_color_g_.get(), def_color_b_.get());
    def_material_->ambient =
      Color(def_color_r_.get(), def_color_g_.get(), def_color_b_.get());
    def_material_->transparency = def_color_a_.get();

    text_material_->diffuse =
      Color(text_color_r_.get(), text_color_g_.get(), text_color_b_.get());
    text_material_->ambient =
      Color(text_color_r_.get(), text_color_g_.get(), text_color_b_.get());
    text_material_->transparency = 1.0;

    
    const TypeDescription *data_type_description =
      fld_handle->get_type_description(Field::BASIS_TD_E);

    if( field_data_basis_type_ != data_type_description->get_name() )
    {
      field_data_basis_type_ = data_type_description->get_name();
      update_algo = true;
    }

    // If the field basis type or the mesh has changed update the algorithm.
    if( update_algo )
    {
      renderer_ = new RenderFieldV;
      BBox bbox = fld_handle->vmesh()->get_bounding_box();

      if (bbox.valid()) 
      {  
        Vector diag = bbox.diagonal(); 
        cur_mesh_scale_factor_ = diag.length();
      } 
      else 
      {
        cur_mesh_scale_factor_ = 1.0;
      }
    
      if (gui_use_default_size_.get() ||
          sci_getenv_p("SCIRUN_USE_DEFAULT_SETTINGS")) 
      {
        set_default_display_values();
      }        

      /*
      // DOTO: REPLACE TO RENDERER
      // NEED TO FIX SCALING IN CUBIC FUNCTIONS FIRST
      const TypeDescription *ftd = fld_handle->get_type_description();
      const TypeDescription *ltd = fld_handle->order_type_description();
      // description for just the data in the field

      // Get the Algorithm.
      CompileInfoHandle ci = RenderFieldBase::get_compile_info(ftd, ltd);
      if (!module_dynamic_compile(ci, renderer_))
      {
        return;
      }

      BBox bbox = fld_handle->mesh()->get_bounding_box();
      
      if (bbox.valid()) 
      {  
        Vector diag = bbox.diagonal(); 
        cur_mesh_scale_factor_ = diag.length();
      } 
      else 
      {
        cur_mesh_scale_factor_ = 1.0;
      }
    
      if (gui_use_default_size_.get() || sci_getenv_p("SCIRUN_USE_DEFAULT_SETTINGS")) 
      {
        set_default_display_values();
      }
      */
      // set new scale defaults based on input.
    }
  }

  // If the colormap gets connected or disconnected then a redraw may
  // be need.  Do this after the algorithm checks so the check does
  // not happen unless needed.

  ColorMapHandle cmap_handle;
  if( get_input_handle("ColorMap", cmap_handle, false) )
  {
    if( !color_map_present_ )
    {
      // If colormap is added and the current selection is for the
      // default assume that the user wants to use the colormap.
      if( nodes_color_type_.get() == 0 ||
          edges_color_type_.get() == 0 ||
          faces_color_type_.get() == 0 )
        remark("Detected a colormap, using it instead of the default color.");

      if( nodes_color_type_.get() == 0 )
        nodes_color_type_.set(1);

      if( edges_color_type_.get() == 0 )
        edges_color_type_.set(1);

      if( faces_color_type_.get() == 0 )
        faces_color_type_.set(1);

      nodes_color_type_.reset();
      edges_color_type_.reset();
      faces_color_type_.reset();

      color_map_present_ = true;
      inputs_changed_ = true;
    }
  }
  else if( color_map_present_ )
  {
    color_map_present_ = false;
    inputs_changed_ = true;
  }

  // Inform module that execution started
  update_state(Executing);

  if( !color_map_present_  &&
      (nodes_color_type_.get() == 1 ||
       edges_color_type_.get() == 1 ||
       faces_color_type_.get() == 1 ||
        text_color_type_.get() == 1) )
  {
    warning("No Colormap present using default color.");

    if( nodes_color_type_.get() == 1)
      nodes_color_type_.set(0);
    if( edges_color_type_.get() == 1)
      edges_color_type_.set(0);
    if( faces_color_type_.get() == 1)
      faces_color_type_.set(0);
    if( text_color_type_.get() == 1)
      text_color_type_.set(0);

    nodes_color_type_.reset();
    edges_color_type_.reset();
    faces_color_type_.reset();
    text_color_type_.reset();
  }

  // Do this after the algorithm checks so the check does not happen
  // unless needed.
  if( gui_field_name_.changed() )
    inputs_changed_ = true;

  // Major input change so everything is dirty
  if( inputs_changed_ )
  {
    set_flag( render_state_[NODE], DIRTY);
    set_flag( render_state_[EDGE], DIRTY);
    set_flag( render_state_[FACE], DIRTY);
    set_flag( render_state_[TEXT], DIRTY);
  }

  // check to see if there is something to do.
  if (!get_flag(render_state_[NODE], DIRTY) && 
      !get_flag(render_state_[EDGE], DIRTY) && 
      !get_flag(render_state_[FACE], DIRTY) && 
      !get_flag(render_state_[TEXT], DIRTY))
  {
    return;
  }

  const int dim = fld_handle->vmesh()->dimensionality();
  if (edges_on_.get() && dim < 1)
  {
    remark("Field type contains no edges, not drawing them.");
    edges_on_.set(0);
    edges_on_.reset();
  }

  if (faces_on_.get() && dim < 2)
  {
    remark("Field type contains no faces, not drawing them.");
    faces_on_.set(0);
    faces_on_.reset();
  }

  switch_flag( render_state_[NODE], IS_ON, nodes_on_.get() );
  switch_flag( render_state_[EDGE], IS_ON, edges_on_.get() );
  switch_flag( render_state_[FACE], IS_ON, faces_on_.get() );
  switch_flag( render_state_[TEXT], IS_ON, text_on_.get() );

  switch_flag( render_state_[NODE], USE_DEFAULT_COLOR, nodes_color_type_.get()  == 0);
  switch_flag( render_state_[EDGE], USE_DEFAULT_COLOR, edges_color_type_.get()  == 0);
  switch_flag( render_state_[FACE], USE_DEFAULT_COLOR, faces_color_type_.get()  == 0);
  switch_flag( render_state_[TEXT], USE_DEFAULT_COLOR,  text_color_type_.get()  == 0);

  switch_flag( render_state_[NODE], USE_COLORMAP, nodes_color_type_.get() == 1);
  switch_flag( render_state_[EDGE], USE_COLORMAP, edges_color_type_.get() == 1);
  switch_flag( render_state_[FACE], USE_COLORMAP, faces_color_type_.get() == 1);
  switch_flag( render_state_[TEXT], USE_COLORMAP,  text_color_type_.get() == 1);

  switch_flag( render_state_[NODE], USE_COLOR_CONVERT, nodes_color_type_.get() == 2);
  switch_flag( render_state_[EDGE], USE_COLOR_CONVERT, edges_color_type_.get() == 2);
  switch_flag( render_state_[FACE], USE_COLOR_CONVERT, faces_color_type_.get() == 2);
  switch_flag( render_state_[TEXT], USE_COLOR_CONVERT,  text_color_type_.get() == 2);

  switch_flag( render_state_[NODE], USE_TRANSPARENCY, nodes_transparency_.get() );
  switch_flag( render_state_[EDGE], USE_TRANSPARENCY, edges_transparency_.get() );
  switch_flag( render_state_[FACE], USE_TRANSPARENCY, faces_transparency_.get() );

  switch_flag( render_state_[FACE], USE_NORMALS, faces_normals_.get() );
  switch_flag( render_state_[FACE], USE_TEXTURE, 0 );

  if( faces_on_.get() &&
      fld_handle->vmesh()->is_imagemesh() &&
      faces_usetexture_.get() )
  {
    if( !color_map_present_ )
      warning("No Colormap present can not render images as a texture quad.");
    else if( faces_color_type_.get() == 0 )
      warning("Default color selection overrides rendering images as a texture quad.");
    else
      switch_flag( render_state_[FACE], USE_TEXTURE, faces_usetexture_.get() );
  }

  std::string fname = clean_fieldname(gui_field_name_.get());
  if (fname != "" && fname[fname.size()-1] != ' ') { fname = fname + " "; }

  //  Nodes.
  if (get_flag(render_state_[NODE], (IS_ON|DIRTY))) 
  {
    GeomHandle mesh_geometry =
      renderer_->render_mesh_nodes( fld_handle,
				    nodes_display_type_.get(),
				    nodes_scale_.get(),
				    gui_nodes_resolution_.get(),
				    render_state_[NODE] );

    GeomHandle gmat = new GeomMaterial(mesh_geometry, def_material_);
    GeomHandle geom = new GeomSwitch(new GeomColorMap(gmat, cmap_handle));
    if (nodes_id_) ogeom->delObj(nodes_id_);
    nodes_id_ = ogeom->addObj(geom, fname +
			       (nodes_transparency_.get()?
				"Transparent Nodes":"Nodes"));
  }
  else if (!(render_state_[NODE]&IS_ON))
  {
    if (nodes_id_) ogeom->delObj(nodes_id_);  
    nodes_id_ = 0;
  }

  // Edges.
  if (get_flag(render_state_[EDGE], (IS_ON|DIRTY))) 
  {
    GeomHandle mesh_geometry =
      renderer_->render_mesh_edges( fld_handle,
				    edges_display_type_.get(),
				    edges_scale_.get(),
				    gui_edges_resolution_.get(),
				    render_state_[EDGE],
				    approx_div_.get() );
    
    GeomHandle gmat = new GeomMaterial(mesh_geometry, def_material_);
    GeomHandle geom = new GeomSwitch(new GeomColorMap(gmat, cmap_handle));
    if (edges_id_) ogeom->delObj(edges_id_);
    edges_id_ = ogeom->addObj(geom, fname +
			       (edges_transparency_.get()?
				"Transparent Edges":"Edges"));
  }
  else if (!(render_state_[EDGE]&IS_ON))
  {
    if (edges_id_) ogeom->delObj(edges_id_);
    edges_id_ = 0;
  }

  //  Faces.
  if (get_flag(render_state_[FACE], (IS_ON|DIRTY))) 
  {
    if (faces_normals_.get())
      fld_handle->mesh()->synchronize(Mesh::NORMALS_E);
    
    GeomHandle mesh_geometry =
      renderer_->render_mesh_faces(fld_handle,
				   cmap_handle,
				   render_state_[FACE],
				   approx_div_.get() );
    
    GeomHandle gmat = new GeomMaterial(mesh_geometry, def_material_);
    GeomHandle geom = new GeomSwitch(new GeomColorMap(gmat, cmap_handle));
    if (faces_id_) ogeom->delObj(faces_id_);
    faces_id_ = ogeom->addObj(geom, fname +
			       (faces_transparency_.get()?
				"Transparent Faces":"Faces"));
  }
  else if (!(render_state_[FACE]&IS_ON))
  {
    if (faces_id_) ogeom->delObj(faces_id_);
    faces_id_ = 0;
  }

  // Render Text.
  if (get_flag(render_state_[TEXT], (IS_ON|DIRTY)))
  {
    GeomHandle text_geometry =
      renderer_->render_text(fld_handle,
			     cmap_handle,
			     text_color_type_.get() == 1,
			     text_color_type_.get() == 0,
			     text_backface_cull_.get(),
			     text_fontsize_.get(),
			     text_precision_.get(),
			     text_render_locations_.get(),
			     text_show_data_.get(),
			     text_show_nodes_.get(),
			     text_show_edges_.get(),
			     text_show_faces_.get(),
			     text_show_cells_.get(),
			     text_always_visible_.get());
    
    GeomHandle gmat = new GeomMaterial(text_geometry, text_material_);
    GeomHandle geom = new GeomSwitch(new GeomColorMap(gmat, cmap_handle));
    if (text_id_) ogeom->delObj(text_id_);
    text_id_ = ogeom->addObj(geom, fname +
			      (text_backface_cull_.get()?
			       "Culled Text Data":"Text Data"));
  }
  else if (!(render_state_[TEXT]&IS_ON))
  {
    if (text_id_) ogeom->delObj(text_id_);
    text_id_ = 0;
  }

  clear_flag( render_state_[NODE], DIRTY);
  clear_flag( render_state_[EDGE], DIRTY);
  clear_flag( render_state_[FACE], DIRTY);
  clear_flag( render_state_[TEXT], DIRTY);
}


void
ShowField::set_default_display_values() 
{
  double fact = cur_mesh_scale_factor_;
  nodes_scaleNV_.set(fact * 0.01);
  edges_scaleNV_.set(fact * 0.003);

  set_flag( render_state_[NODE], DIRTY);
  set_flag( render_state_[EDGE], DIRTY);
}


void
ShowField::maybe_execute(toggle_type_e dis_type)
{
  bool do_execute = false;

  interactive_mode_.reset();
  if (interactive_mode_.get() == "Interactive") {
    switch(dis_type) {
    case NODE :
      do_execute = nodes_on_.get();
	break;
    case EDGE :
      do_execute = edges_on_.get();
      break;
    case FACE :
      do_execute = faces_on_.get();
	break;
    case TEXT :
      do_execute = text_on_.get();
	break;
    case ALL :
      do_execute = true;
	break;
    }
  }

  if (do_execute) {
    want_to_execute();
  }
}


void
ShowField::tcl_command(GuiArgs& args, void* userdata)
{
  if(args.count() < 2)
  {
    args.error("ShowField needs a minor command");
    return;
  }

  if (args[1] == "nodes_scale")
  {
    nodes_display_type_.reset();
    if ( nodes_display_type_.get() == "Spheres") {
      set_flag( render_state_[NODE], DIRTY);
      maybe_execute(NODE);
    }
  } else if (args[1] == "edges_scale") {
    edges_display_type_.reset();
    if ( edges_display_type_.get() == "Cylinders") {
      set_flag( render_state_[EDGE], DIRTY);
      maybe_execute(EDGE);
    }
  } else if (args[1] == "approx") {
    set_flag( render_state_[EDGE], DIRTY);
    set_flag( render_state_[FACE], DIRTY);

  } else if (args[1] == "default_color_change") {
    def_color_r_.reset();
    def_color_g_.reset();
    def_color_b_.reset();
    def_color_a_.reset();
    def_material_->diffuse =
      Color(def_color_r_.get(), def_color_g_.get(), def_color_b_.get());
    def_material_->ambient =
      Color(def_color_r_.get(), def_color_g_.get(), def_color_b_.get());
    def_material_->transparency = def_color_a_.get();
    GeometryOPortHandle ogeom;
    get_oport_handle("Scene Graph",ogeom);
    if (ogeom.get_rep()) ogeom->flushViews();
  } else if (args[1] == "text_color_change") {
    text_color_r_.reset();
    text_color_g_.reset();
    text_color_b_.reset();
    text_material_->diffuse =
      Color(text_color_r_.get(), text_color_g_.get(), text_color_b_.get());
    text_material_->ambient =
      Color(text_color_r_.get(), text_color_g_.get(), text_color_b_.get());
    GeometryOPortHandle ogeom;
    get_oport_handle("Scene Graph",ogeom);
    if( ogeom.get_rep() ) ogeom->flushViews();

  } 
  else if (args[1] == "toggle_display_nodes") 
  {
    // Toggle the GeomSwitches.
    nodes_on_.reset();
    if (nodes_on_.get())
    {
      set_flag( render_state_[NODE], DIRTY);
      maybe_execute(NODE);
    }
    else if (nodes_id_)
    {
      GeometryOPortHandle ogeom;
      get_oport_handle("Scene Graph",ogeom);
      ogeom->delObj(nodes_id_);
      ogeom->flushViews();
      nodes_id_ = 0;
    }
  } 
  else if (args[1] == "toggle_display_edges") 
  {
    // Toggle the GeomSwitch.
    edges_on_.reset();
    if (edges_on_.get())
    {
      set_flag( render_state_[EDGE], DIRTY);
      maybe_execute(EDGE);
    }
    else if (edges_id_)
    {
      GeometryOPortHandle ogeom;
      get_oport_handle("Scene Graph",ogeom);
      ogeom->delObj(edges_id_);
      ogeom->flushViews();
      edges_id_ = 0;
    }
  } 
  else if (args[1] == "toggle_display_faces") 
  {
    // Toggle the GeomSwitch.
    faces_on_.reset();
    if (faces_on_.get())
    {
      set_flag( render_state_[FACE], DIRTY);
      maybe_execute(FACE);
    }
    else if (faces_id_)
    {
      GeometryOPortHandle ogeom;
      get_oport_handle("Scene Graph",ogeom);
      ogeom->delObj(faces_id_);
      ogeom->flushViews();
      faces_id_ = 0;
    }
  } 
  else if (args[1] == "toggle_display_text")
  {
    // Toggle the GeomSwitch.
    text_on_.reset();
    if ((text_on_.get()) && (text_id_ == 0))
    {
      set_flag( render_state_[TEXT], DIRTY);
      maybe_execute(TEXT);
    }
    else if (!text_on_.get() && text_id_)
    {
      GeometryOPortHandle ogeom;
      get_oport_handle("Scene Graph",ogeom);
      ogeom->delObj(text_id_);
      ogeom->flushViews();
      text_id_ = 0;
    }
  } 
  else if (args[1] == "rerender_all") 
  {
    set_flag( render_state_[NODE], DIRTY);
    set_flag( render_state_[EDGE], DIRTY);
    set_flag( render_state_[FACE], DIRTY);
    set_flag( render_state_[TEXT], DIRTY);
    maybe_execute(ALL);
  } 
  else if (args[1] == "rerender_nodes") 
  {
    set_flag( render_state_[NODE], DIRTY);
    maybe_execute(NODE);
  } 
  else if (args[1] == "rerender_edges") 
  {
    set_flag( render_state_[EDGE], DIRTY);
    maybe_execute(EDGE);
  } 
  else if (args[1] == "rerender_faces") 
  {
    set_flag( render_state_[FACE], DIRTY);
    maybe_execute(FACE);
  } 
  else if (args[1] == "rerender_text") 
  {
    set_flag( render_state_[TEXT], DIRTY);
    maybe_execute(TEXT);
  } 
  else if (args[1] == "execute_policy") 
  {
  } 
  else if (args[1] == "calcdefs") 
  {
    set_default_display_values();
    maybe_execute(ALL);
  } 
  else 
  {
    Module::tcl_command(args, userdata);
  }
}

void
ShowField::post_read()
{
  char names[38][2][30] = { {"nodes-on", "nodes_on"},
			    {"nodes-transparency","nodes_transparency" },
//		            {"nodes-as-disks", "" },

// Special case 	    {"nodes-usedefcolor", "nodes_color_type" },
			    {"node_display_type", "nodes_display_type" },
			    
			    {"edges-on", "edges_on" },
			    {"edges-transparency", "edges_transparency" },
// Special case		    {"edges-usedefcolor", "edges_color_type" },
			    {"edge_display_type", "edges_display_type" },
			    
			    {"faces-on", "faces_on" },
			    {"use-normals", "faces_normals" },
			    {"use-transparency", "faces_transparency" },
// Special case		    {"faces-usedefcolor", "faces_color_type" },
			    {"faces-usetexture", "faces_usetexture" },
			    
// 		            {"vectors-on",  },
// 		            {"normalize-vectors",  },
// 		            {"has_vector_data",  },
// 		            {"bidirectional",  },
// 		            {"vectors-usedefcolor",  },

// 		            {"tensors-on",  },
// 		            {"has_tensor_data",  },
// 		            {"tensors-usedefcolor",  },
// 		            {"tensors-emphasis",  },

// 		            {"scalars-on",  },
// 		            {"scalars-transparency",  },
// 		            {"scalars-usedefcolor",  },
// 		            {"has_scalar_data",  },

			    {"text-on", "text_on" },
			    {"text-use-default-color", "text_color_type" },
			    {"text-color-r", "text_color-r" },
			    {"text-color-g", "text_color-g" },
			    {"text-color-b", "text_color-b" },
			    {"text-backface-cull", "text_backface_cull" },
			    {"text-always_visible", "text_always_visible" },
			    {"text-fontsize", "text_fontsize" },
			    {"text-precision","text_precision"  },
			    {"text-render_locations", "text_render_locations" },
			    {"text-show-data",  "text_show_data" },
			    {"text-show-nodes", "text_show_nodes" },
			    {"text-show-edges", "text_show_edges" },
			    {"text-show-faces", "text_show_faces" },
			    {"text-show-cells", "text_show_cells" },
		       
			    {"def-color-r", "def_color-r" },
			    {"def-color-g", "def_color-g" },
			    {"def-color-b", "def_color-b" },
// Special case		    {"def-color-a", "def_color-a" },

// 		            {"data_display_type",  },
// 		            {"tensor_display_type",  },
// 		            {"scalar_display_type",  },

			    {"node_scale",   "nodes_scale" },
			    {"node_scaleNV", "nodes_scaleNV" },
			    {"edge_scale",   "edges_scale" },
			    {"edge_scaleNV", "edges_scaleNV"  },

// 		            {"vectors_scale",  },
// 		            {"vectors_scaleNV",  },
// 		            {"tensors_scale",  },
// 		            {"tensors_scaleNV",  },
// 		            {"scalars_scale",  },
// 		            {"scalars_scaleNV",  },

// No Change		    {"active_tab", "active_tab" },
// No Change		    {"interactive_mode", "interactive_mode" },
// No Change		    {"show_progress", "show_progress" },
			    {"field-name","field_name"  },
			    {"field-name-override","field_name_override"  },
//		            {"field-name-update", },
			    {"node-resolution", "nodes_resolution" },
			    {"edge-resolution", "edges_resolution" },
//		            {"data-resolution",  },
			    {"approx-div", "approx_div" },
			    {"use-default-size", "use_default_size" } };

  // Get the module name
  const std::string modName = get_ctx()->getfullname() + "-";

  std::string val;

  for( unsigned int i=0; i<38; i++ )
  {
    // Get the current values for the old names
    if( TCLInterface::get(modName+names[i][0], val, get_ctx()) )
    {
      // Set the current values for the new names
      TCLInterface::set(modName+names[i][1], val, get_ctx());
    }
  }

  bool has_colormap = true; //get_oport("ColorMap")->nconnections();

  // Special cases for the old default color vars
  if( TCLInterface::get(modName+"nodes-usedefcolor", val, get_ctx()) )
  {
    if( val == std::string("1") )
      TCLInterface::set(modName+"nodes_color_type", "0", get_ctx());
    else if( has_colormap )
      TCLInterface::set(modName+"nodes_color_type", "1", get_ctx());
    else
      TCLInterface::set(modName+"nodes_color_type", "2", get_ctx());
  }

  if( TCLInterface::get(modName+"edges-usedefcolor", val, get_ctx()) )
  {
    if( val == std::string("1") )
      TCLInterface::set(modName+"edges_color_type", "0", get_ctx());
    else if( has_colormap )
      TCLInterface::set(modName+"nodes_color_type", "1", get_ctx());
    else
      TCLInterface::set(modName+"nodes_color_type", "2", get_ctx());
  }

  if( TCLInterface::get(modName+"faces-usedefcolor", val, get_ctx()) )
  {
    if( val == std::string("1") )
      TCLInterface::set(modName+"faces_color_type", "0", get_ctx());
    else if( has_colormap )
      TCLInterface::set(modName+"nodes_color_type", "1", get_ctx());
    else
      TCLInterface::set(modName+"nodes_color_type", "2", get_ctx());
  }

  if( TCLInterface::get(modName+"def-color-a", val, get_ctx()) )
  {
    double dval;
    from_string(val,dval);

    dval = dval * dval * dval * dval;
    val = to_string( dval );

    TCLInterface::set(modName+"def_color-a", val, get_ctx());
  }
}

DECLARE_MAKER(ShowField)
} // End namespace SCIRun
