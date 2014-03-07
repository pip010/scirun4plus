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
 *  ShowFieldGlyphs.cc
 *
 *  Written by:
 *   Allen R. Sanderson
 *   Marty Cole
 *   School of Computing
 *   University of Utah
 *   March 15, 2007
 *
 */


#include <Core/Datatypes/Field.h>
#include <Core/Geom/GeomMaterial.h>
#include <Core/Geom/GeomSwitch.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <Core/Algorithms/Fields/CompareFields/SimilarMeshes.h>
#include <Core/Algorithms/Visualization/RenderField.h>
#include <Core/Algorithms/Visualization/RenderFieldGlyphs.h>
#include <Core/Algorithms/Util/FieldInformation.h>


#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/ColorMapPort.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Network/Ports/FieldPort.h>

#include <typeinfo>
#include <iostream>

namespace SCIRun {

using namespace SCIRunAlgo;

class ShowFieldGlyphs : public Module, public RenderStateBase
{
  public:
    ShowFieldGlyphs(GuiContext* ctx);
    virtual ~ShowFieldGlyphs();

    virtual void execute();
    virtual void tcl_command(GuiArgs& args, void* userdata);

    virtual void post_read();

  protected:
    //! input ports
    int                      pfld_mesh_generation_;
    int                      sfld_mesh_generation_;
    int                      tfld_mesh_generation_;
    bool                     secondary_field_present_;
    bool                     tertiary_field_present_;
    bool                     color_map_present_;

    //! Scene graph ID's
    int                      data_id_;
    int                      text_id_;

    //! top level nodes for switching on and off.

    //! Options for rendering data.
    GuiInt                   scalars_has_data_;
    GuiInt                   scalars_on_;
    GuiString                scalars_display_type_;
    GuiInt                   scalars_transparency_;
    GuiInt                   scalars_normalize_;
    GuiInt                   scalars_color_type_;
    GuiDouble                scalars_scale_;
    GuiDouble                scalars_scaleNV_;
    GuiDouble                scalars_threshold_;
    GuiDouble                scalars_thresholdNV_;
    GuiInt                   scalars_resolution_;
    GuiInt                   scalars_small_is_dot_;

    GuiInt                   vectors_has_data_;
    GuiInt                   vectors_on_;
    GuiString                vectors_display_type_;
    GuiInt                   vectors_transparency_;
    GuiInt                   vectors_normalize_;
    GuiInt                   vectors_color_type_;
    GuiDouble                vectors_scale_;
    GuiDouble                vectors_scaleNV_;
    GuiDouble                vectors_threshold_;
    GuiDouble                vectors_thresholdNV_;
    GuiInt                   vectors_resolution_;
    GuiInt                   vectors_bidirectional_;
    GuiInt                   vectors_small_is_dot_;

    GuiInt                   tensors_has_data_;
    GuiInt                   tensors_on_;
    GuiString                tensors_display_type_;
    GuiInt                   tensors_transparency_;
    GuiInt                   tensors_normalize_;
    GuiInt                   tensors_color_type_;
    GuiDouble                tensors_scale_;
    GuiDouble                tensors_scaleNV_;
    GuiDouble                tensors_threshold_;
    GuiDouble                tensors_thresholdNV_;
    GuiInt                   tensors_resolution_;
    GuiDouble                tensors_emphasis_;
    GuiInt                   tensors_small_is_dot_;

    GuiInt                   secondary_has_data_;
    GuiInt                   secondary_on_;
    GuiString                secondary_display_type_;
    GuiInt                   secondary_color_type_;
    GuiInt                   secondary_alpha_;
    GuiInt                   secondary_value_;
    GuiDouble                secondary_scale_;
    GuiDouble                secondary_scaleNV_;

    GuiInt                   tertiary_has_data_;
    GuiInt                   tertiary_on_;
    GuiString                tertiary_display_type_;
    GuiInt                   tertiary_color_type_;
    GuiInt                   tertiary_alpha_;
    GuiInt                   tertiary_value_;
    GuiDouble                tertiary_scale_;
    GuiDouble                tertiary_scaleNV_;

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

    GuiString                active_tab_; //! for saving nets state
    GuiString                interactive_mode_;
    GuiInt                   showProgress_;

    GuiString                gui_field_name_;
    GuiInt                   gui_field_name_override_;

    //! Refinement resolution for cylinders and spheres
    GuiInt                   approx_div_;

    LockingHandle<RenderFieldBase>   text_renderer_;
    LockingHandle<RenderScalarField> scalar_renderer_;
    LockingHandle<RenderVectorField> vector_renderer_;
    LockingHandle<RenderTensorField> tensor_renderer_;

    // variables related to default scale factor usage.
    GuiInt                  gui_use_default_size_;
    double                  cur_mesh_scale_factor_;

    enum toggle_type_e {
      SCALAR    = 0,
      VECTOR    = 1,
      TENSOR    = 2,
      SECONDARY = 3,
      TERTIARY  = 4,
      TEXT      = 5,
      ALL   = 6
    };

    unsigned int            render_state_[ALL];

    void maybe_execute(toggle_type_e dis_type);
    void set_default_display_values();
    
    SCIRunAlgo::SimilarMeshesAlgo algo_;
};

ShowFieldGlyphs::ShowFieldGlyphs(GuiContext* ctx) :
  Module("ShowFieldGlyphs", ctx, Filter, "Visualization", "SCIRun"),
  pfld_mesh_generation_(-1),
  sfld_mesh_generation_(-1),
  tfld_mesh_generation_(-1),
  secondary_field_present_(false),
  tertiary_field_present_(false),
  color_map_present_(false),
  data_id_(0),
  text_id_(0),

  scalars_has_data_(get_ctx()->subVar("scalars_has_data"), 0),
  scalars_on_(get_ctx()->subVar("scalars_on"), 0),
  scalars_display_type_(get_ctx()->subVar("scalars_display_type"), "Spheres"),
  scalars_transparency_(get_ctx()->subVar("scalars_transparency"), 0),
  scalars_normalize_(get_ctx()->subVar("scalars_normalize"), 0),
  scalars_color_type_(get_ctx()->subVar("scalars_color_type"), 1),
  scalars_scale_(get_ctx()->subVar("scalars_scale")),
  scalars_scaleNV_(get_ctx()->subVar("scalars_scaleNV")),
  scalars_threshold_(get_ctx()->subVar("scalars_threshold"),1.0e-5),
  scalars_thresholdNV_(get_ctx()->subVar("scalars_thresholdNV")),
  scalars_resolution_(get_ctx()->subVar("scalars_resolution"), 6),
  scalars_small_is_dot_(get_ctx()->subVar("scalars_small_is_dot"), 1),

  vectors_has_data_(get_ctx()->subVar("vectors_has_data"), 0),
  vectors_on_(get_ctx()->subVar("vectors_on"), 0),
  vectors_display_type_(get_ctx()->subVar("vectors_display_type"), "Arrows"),
  vectors_transparency_(get_ctx()->subVar("vectors_transparency"), 0),
  vectors_normalize_(get_ctx()->subVar("vectors_normalize"), 0),
  vectors_color_type_(get_ctx()->subVar("vectors_color_type"), 1),
  vectors_scale_(get_ctx()->subVar("vectors_scale")),
  vectors_scaleNV_(get_ctx()->subVar("vectors_scaleNV")),
  vectors_threshold_(get_ctx()->subVar("vectors_threshold"),1.0e-5),
  vectors_thresholdNV_(get_ctx()->subVar("vectors_thresholdNV")),
  vectors_resolution_(get_ctx()->subVar("vectors_resolution"), 6),
  vectors_bidirectional_(get_ctx()->subVar("vectors_bidirectional"), 0),
  vectors_small_is_dot_(get_ctx()->subVar("vectors_small_is_dot"), 1),

  tensors_has_data_(get_ctx()->subVar("tensors_has_data"), 0),
  tensors_on_(get_ctx()->subVar("tensors_on"), 0),
  tensors_display_type_(get_ctx()->subVar("tensors_display_type"), "Colored Boxes"),
  tensors_transparency_(get_ctx()->subVar("tensors_transparency"), 0),
  tensors_normalize_(get_ctx()->subVar("tensors_normalize"), 0),
  tensors_color_type_(get_ctx()->subVar("tensors_color_type"), 1),
  tensors_scale_(get_ctx()->subVar("tensors_scale")),
  tensors_scaleNV_(get_ctx()->subVar("tensors_scaleNV")),
  tensors_threshold_(get_ctx()->subVar("tensors_threshold"),1.0e-5),
  tensors_thresholdNV_(get_ctx()->subVar("tensors_thresholdNV")),
  tensors_resolution_(get_ctx()->subVar("tensors_resolution"), 6),
  tensors_emphasis_(get_ctx()->subVar("tensors_emphasis"), 0.825),
  tensors_small_is_dot_(get_ctx()->subVar("tensors_small_is_dot"), 1),
  
  secondary_has_data_(get_ctx()->subVar("secondary_has_data"), 0),
  secondary_on_(get_ctx()->subVar("secondary_on"), 0),
  secondary_display_type_(get_ctx()->subVar("secondary_display_type"), "Major Radius"),
  secondary_color_type_(get_ctx()->subVar("secondary_color_type"), 0),
  secondary_alpha_(get_ctx()->subVar("secondary_alpha"), 0),
  secondary_value_(get_ctx()->subVar("secondary_value"), 1),
  secondary_scale_(get_ctx()->subVar("secondary_scale")),
  secondary_scaleNV_(get_ctx()->subVar("secondary_scaleNV")),

  tertiary_has_data_(get_ctx()->subVar("tertiary_has_data"), 0),
  tertiary_on_(get_ctx()->subVar("tertiary_on"), 0),
  tertiary_display_type_(get_ctx()->subVar("tertiary_display_type"), "Minor Radius"),
  tertiary_color_type_(get_ctx()->subVar("tertiary_color_type"), 0),
  tertiary_alpha_(get_ctx()->subVar("tertiary_alpha"), 0),
  tertiary_value_(get_ctx()->subVar("tertiary_value"), 1),
  tertiary_scale_(get_ctx()->subVar("tertiary_scale")),
  tertiary_scaleNV_(get_ctx()->subVar("tertiary_scaleNV")),

  text_on_(get_ctx()->subVar("text_on"), 0),
  text_color_type_(get_ctx()->subVar("text_color_type"), 0),
  text_color_r_(get_ctx()->subVar("text_color-r"), 1.0),
  text_color_g_(get_ctx()->subVar("text_color-g"), 1.0),
  text_color_b_(get_ctx()->subVar("text_color-b"), 1.0),
  text_backface_cull_(get_ctx()->subVar("text_backface_cull"), 0),
  text_always_visible_(get_ctx()->subVar("text_always_visible"), 0),
  text_fontsize_(get_ctx()->subVar("text_fontsize"), 0),
  text_precision_(get_ctx()->subVar("text_precision"),3),
  text_render_locations_(get_ctx()->subVar("text_render_locations"), 0),
  text_show_data_(get_ctx()->subVar("text_show_data"),1),
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
  def_material_(new Material(Color(0.5, 0.5, 0.5))),

  active_tab_(get_ctx()->subVar("active_tab"), "Scalars"),
  interactive_mode_(get_ctx()->subVar("interactive_mode"), "Interactive"),
  showProgress_(get_ctx()->subVar("show_progress"), 0),

  gui_field_name_(get_ctx()->subVar("field_name"), ""),
  gui_field_name_override_(get_ctx()->subVar("field_name_override"), 0),

  approx_div_(get_ctx()->subVar("approx_div"), 1),

  text_renderer_(0),
  scalar_renderer_(0),
  vector_renderer_(0),
  tensor_renderer_(0),

  gui_use_default_size_(get_ctx()->subVar("use_default_size"), 0),
  cur_mesh_scale_factor_(1.0)
{
  render_state_[SCALAR] = 0;
  render_state_[VECTOR] = 0;
  render_state_[TENSOR] = 0;
  render_state_[SECONDARY] = 0;
  render_state_[TERTIARY] = 0;
  render_state_[TEXT] = 0;

  algo_.set_progress_reporter(this);
}


ShowFieldGlyphs::~ShowFieldGlyphs()
{
}


void
ShowFieldGlyphs::execute()
{
  bool update_algo = false;

  GeometryOPortHandle ogeom;
  get_oport_handle("Scene Graph",ogeom);

  FieldHandle pfld_handle;
  if( !get_input_handle("Primary Data", pfld_handle, !in_power_app() ) )
    return;

  // Get the information on the types of the primary field.
  FieldInformation pfi(pfld_handle);

  if (!pfi.is_svt() )
  {
    error("No Scalar, Vector, or Tensor data found in the data field.");
    return;
  }

  // Set has data vars so the GUI comes up correctly.
  scalars_has_data_.set(pfi.is_scalar());
  vectors_has_data_.set(pfi.is_vector());
  tensors_has_data_.set(pfi.is_tensor());

  if(scalars_on_.get() == 0 &&
     vectors_on_.get() == 0 &&
     tensors_on_.get() == 0 &&
        text_on_.get() == 0 )
    return;

  // If the field has changed the algo must be updated.
  if( inputs_changed_ == true ) update_algo = true;

  // Update the field name but only if the user does not enter an
  // overriding name.
  if (!gui_field_name_override_.get())
  {
    std::string fname("");

    if ( !pfld_handle->get_property("name", fname))
      pfld_handle->mesh()->get_property("name", fname);
       
    if( fname != gui_field_name_.get() )
    {
      gui_field_name_.set(fname);
      gui_field_name_.reset();
    }
  }

  // See if the mesh has changed because this is not check with the ports.
  if( pfld_mesh_generation_ != pfld_handle->mesh()->generation )
  {
    pfld_mesh_generation_ = pfld_handle->mesh()->generation;
    inputs_changed_ = true;
    update_algo = true;
  }

  // Get the optional secondary field that will be used for the color and
  // secondary scaling.
  FieldHandle sfld_handle;
  if( get_input_handle("Secondary Data", sfld_handle, false) )
  {
    // If the field has changed the algo must be updated.
    if( inputs_changed_ == true )
      update_algo = true;


    if( !secondary_field_present_ )
    {
      secondary_field_present_ = true;
      inputs_changed_ = true;
      update_algo = true;
    }

    // See if the mesh has changed because this is not check with the ports.
    else if( sfld_mesh_generation_ != sfld_handle->mesh()->generation )
    {
      sfld_mesh_generation_ = sfld_handle->mesh()->generation;
      inputs_changed_ = true;
      update_algo = true;
    }
  }
  else 
  {
    if( secondary_field_present_ )
    {
      secondary_field_present_ = false;
      inputs_changed_ = true;
      update_algo = true;
    }

    sfld_handle = pfld_handle;
  }

  // Get the optional tertiary field that will be used for the color and
  // tertiary scaling.
  FieldHandle tfld_handle;
  if( get_input_handle("Tertiary Data", tfld_handle, false) )
  {
    // If the field has changed the algo must be updated.
    if( inputs_changed_ == true )
      update_algo = true;

    if( !tertiary_field_present_ )
    {
      tertiary_field_present_ = true;
      inputs_changed_ = true;
      update_algo = true;
    }

    // See if the mesh has changed because this is not check with the ports.
    else if( tfld_mesh_generation_ != tfld_handle->mesh()->generation )
    {
      tfld_mesh_generation_ = tfld_handle->mesh()->generation;
      inputs_changed_ = true;
      update_algo = true;
    }
  }
  else 
  {
    if( tertiary_field_present_ )
    {
      tertiary_field_present_ = false;
      inputs_changed_ = true;
      update_algo = true;
    }

    tfld_handle = pfld_handle;
  }

  // See if the algorithm needs to be recreated due the field or mesh
  // changing. This check will be made if colormap changes also.
  if( inputs_changed_ ) 
  {
    // Set inital colors here.
    def_material_->diffuse =
      Color(def_color_r_.get(), def_color_g_.get(), def_color_b_.get());
    def_material_->transparency = def_color_a_.get();

    text_material_->diffuse =
      Color(text_color_r_.get(), text_color_g_.get(), text_color_b_.get());
    text_material_->transparency = 1.0;
    
    // Make sure both fields have similar meshes and data locations.
    if( sfld_handle != pfld_handle || tfld_handle != pfld_handle )
    {
      // Get the information on the types of the fields
      FieldInformation sfi(sfld_handle);
      FieldInformation tfi(tfld_handle);
      
      if (!sfi.is_svt() )
      {
        error("No Scalar, Vector, or Tensor data found in the secondary field.");
        return;
      }

      if (!tfi.is_svt() )
      {
        error("No Scalar, Vector, or Tensor data found in the tertiary field.");
        return;
      }

      std::vector< FieldHandle > field_input_handles;
	
      field_input_handles.push_back( pfld_handle );
      field_input_handles.push_back( sfld_handle );
      field_input_handles.push_back( tfld_handle );
      
      if( !(algo_.run(field_input_handles ))) return;
    }
    else
    {
      if( !secondary_field_present_ ) 
      {
        secondary_on_.set(0);
        secondary_on_.reset();
      }

      if( !tertiary_field_present_ ) 
      {
        tertiary_on_.set(0);
        tertiary_on_.reset();
      }
    }

    secondary_has_data_.set(secondary_field_present_);
    secondary_has_data_.reset();

    tertiary_has_data_.set(tertiary_field_present_);
    tertiary_has_data_.reset();

    scalars_has_data_.reset();
    vectors_has_data_.reset();
    tensors_has_data_.reset();
    secondary_has_data_.reset();
    tertiary_has_data_.reset();

    switch_flag( render_state_[SCALAR], HAS_DATA, pfi.is_scalar() );
    switch_flag( render_state_[VECTOR], HAS_DATA, pfi.is_vector() );
    switch_flag( render_state_[TENSOR], HAS_DATA, pfi.is_tensor() );
    switch_flag( render_state_[SECONDARY], HAS_DATA, secondary_field_present_ );
    switch_flag( render_state_[TERTIARY], HAS_DATA, tertiary_field_present_ );

    const TypeDescription *data_type_description =
      pfld_handle->get_type_description(Field::BASIS_TD_E);

    if( field_data_basis_type_ != data_type_description->get_name() )
    {
      field_data_basis_type_ = data_type_description->get_name();
      update_algo = true;
    }

    // If the field basis type or the mesh has changed update the algorithm.
    if (update_algo)
    {
      text_renderer_ = new RenderFieldV;
      scalar_renderer_ = new RenderScalarField;
      vector_renderer_ = new RenderVectorField;
      tensor_renderer_ = new RenderTensorField;
      
      BBox bbox = pfld_handle->vmesh()->get_bounding_box();
   
      if (bbox.valid()) 
      {  
        Vector diag = bbox.diagonal(); 
        cur_mesh_scale_factor_ = diag.length();
      } 
      else 
      {
        cur_mesh_scale_factor_ = 1.0;
      }
    
      gui_use_default_size_.reset();
      if (gui_use_default_size_.get() || 
          sci_getenv_p("SCIRUN_USE_DEFAULT_SETTINGS")) 
      {
        set_default_display_values();
      }        
    }
  }

  // If the colormap gets connected or disconnected then a redraw may
  // be need.  Do this after the algorithm checks so the
  // inputs_changed_ is not affected.

  bool need_primary_color_map =
    ((scalars_on_.get() == 1 && scalars_color_type_.get() == 1) ||
     (vectors_on_.get() == 1 && vectors_color_type_.get() == 1) ||
     (tensors_on_.get() == 1 && tensors_color_type_.get() == 1) ||
     (   text_on_.get() == 1 &&    text_color_type_.get() == 1));

  bool need_secondary_color_map =
    (secondary_on_.get() == 1 && secondary_color_type_.get() == 1);

  bool need_tertiary_color_map =
    (tertiary_on_.get() == 1 && tertiary_color_type_.get() == 1);

  ColorMapHandle cmap_handle;

  // Check the tertiary first as it overrides the secondary.
  if( need_tertiary_color_map )
  {
    if( get_input_handle("Tertiary ColorMap", cmap_handle, true) )
    {
      if( !color_map_present_ )
      {
        color_map_present_ = true;
        inputs_changed_ = true;
      }
    }
    else if( color_map_present_ )
    {
      color_map_present_ = false;
      inputs_changed_ = true;
    }
  }

  // Check the secondary first as it overrides the primary.
  else if( need_secondary_color_map )
  {
    if( get_input_handle("Secondary ColorMap", cmap_handle, true) )
    {
      if( !color_map_present_ )
      {
        color_map_present_ = true;
        inputs_changed_ = true;
      }
    }
    else if( color_map_present_ )
    {
      color_map_present_ = false;
      inputs_changed_ = true;
    }
  }

  else if( need_primary_color_map )
  {
    if( get_input_handle("Primary ColorMap", cmap_handle, true) )
    {
      if( !color_map_present_ )
      {
          // If colormap is added and the current selection is for the
          // default assume that the user wants to use the colormap.
          if( scalars_color_type_.get() == 0 ||
              vectors_color_type_.get() == 0 ||
              tensors_color_type_.get() == 0 )
            remark("Detected a colormap, using it instead of the default color.");

          if( scalars_color_type_.get() == 0 )
            scalars_color_type_.set(1);

          if( vectors_color_type_.get() == 0 )
            vectors_color_type_.set(1);

          if( tensors_color_type_.get() == 0 )
            tensors_color_type_.set(1);

          scalars_color_type_.reset();
          vectors_color_type_.reset();
          tensors_color_type_.reset();

          color_map_present_ = true;
          inputs_changed_ = true;
      }
    }
    else if( color_map_present_ )
    {
      color_map_present_ = false;
      inputs_changed_ = true;
    }
  }

  // Inform module that execution started
  update_state(Executing);


  if( !color_map_present_  &&
      (need_primary_color_map ||
       need_secondary_color_map ||
       need_tertiary_color_map) )
  {
    warning("No Colormap present using default color.");

    if(scalars_color_type_.get() == 1)
      scalars_color_type_.set(0);
    if(vectors_color_type_.get() == 1)
      vectors_color_type_.set(0);
    if(tensors_color_type_.get() == 1)
      tensors_color_type_.set(0);
    if(secondary_color_type_.get() == 1)
      secondary_color_type_.set(0);
    if(tertiary_color_type_.get() == 1)
      tertiary_color_type_.set(0);
    if(text_color_type_.get() == 1)
      text_color_type_.set(0);

    scalars_color_type_.reset();
    vectors_color_type_.reset();
    tensors_color_type_.reset();
    secondary_color_type_.reset();
    tertiary_color_type_.reset();
    text_color_type_.reset();
  }

  // Do this after the algorithm checks so the check does not happen
  // unless needed.
  if( gui_field_name_.changed() )
    inputs_changed_ = true;

  // Major input change so everything is dirty
  if( inputs_changed_ )
  {
    set_flag( render_state_[SCALAR], DIRTY);
    set_flag( render_state_[VECTOR], DIRTY);
    set_flag( render_state_[TENSOR], DIRTY);
    set_flag( render_state_[TEXT  ], DIRTY);
  }

  // check to see if there is something to do.
  if (!get_flag(render_state_[SCALAR], DIRTY) && 
      !get_flag(render_state_[VECTOR], DIRTY) && 
      !get_flag(render_state_[TENSOR], DIRTY) && 
      !get_flag(render_state_[TEXT], DIRTY))
  {
    return;
  }

  switch_flag( render_state_[SCALAR], IS_ON, scalars_on_.get() );
  switch_flag( render_state_[VECTOR], IS_ON, vectors_on_.get() );
  switch_flag( render_state_[TENSOR], IS_ON, tensors_on_.get() );
  switch_flag( render_state_[SECONDARY], IS_ON, secondary_on_.get() );
  switch_flag( render_state_[TERTIARY], IS_ON, tertiary_on_.get() );
  switch_flag( render_state_[TEXT  ], IS_ON, text_on_.get() );

  switch_flag( render_state_[SCALAR],    USE_COLORMAP, color_map_present_ );
  switch_flag( render_state_[VECTOR],    USE_COLORMAP, color_map_present_ );
  switch_flag( render_state_[TENSOR],    USE_COLORMAP, color_map_present_ );
  switch_flag( render_state_[SECONDARY], USE_COLORMAP, color_map_present_ );
  switch_flag( render_state_[TERTIARY], USE_COLORMAP, color_map_present_ );
  switch_flag( render_state_[TEXT  ],    USE_COLORMAP, color_map_present_ );

  switch_flag( render_state_[SCALAR], USE_TRANSPARENCY, scalars_transparency_.get() );
  switch_flag( render_state_[VECTOR], USE_TRANSPARENCY, vectors_transparency_.get() );
  switch_flag( render_state_[TENSOR], USE_TRANSPARENCY, tensors_transparency_.get() );

  switch_flag( render_state_[SCALAR], USE_DEFAULT_COLOR, scalars_color_type_.get() == 0 );
  switch_flag( render_state_[VECTOR], USE_DEFAULT_COLOR, vectors_color_type_.get() == 0 );
  switch_flag( render_state_[TENSOR], USE_DEFAULT_COLOR, tensors_color_type_.get() == 0 );
  switch_flag( render_state_[TEXT  ], USE_DEFAULT_COLOR,    text_color_type_.get() == 0 );
  switch_flag( render_state_[SECONDARY], USE_DEFAULT_COLOR, secondary_color_type_.get() == 0 );
  switch_flag( render_state_[TERTIARY], USE_DEFAULT_COLOR, tertiary_color_type_.get() == 0 );

  switch_flag( render_state_[SCALAR], USE_COLORMAP, scalars_color_type_.get() == 1 );
  switch_flag( render_state_[VECTOR], USE_COLORMAP, vectors_color_type_.get() == 1 );
  switch_flag( render_state_[TENSOR], USE_COLORMAP, tensors_color_type_.get() == 1 );
  switch_flag( render_state_[TEXT  ], USE_COLORMAP,    text_color_type_.get() == 1 );
  switch_flag( render_state_[SECONDARY], USE_COLORMAP, secondary_color_type_.get() == 1 );
  switch_flag( render_state_[TERTIARY], USE_COLORMAP, tertiary_color_type_.get() == 1 );

  switch_flag( render_state_[SCALAR], USE_COLOR_CONVERT, scalars_color_type_.get() == 2 );
  switch_flag( render_state_[VECTOR], USE_COLOR_CONVERT, vectors_color_type_.get() == 2 );
  switch_flag( render_state_[TENSOR], USE_COLOR_CONVERT, tensors_color_type_.get() == 2 );
  switch_flag( render_state_[TEXT  ], USE_COLOR_CONVERT,    text_color_type_.get() == 2 );
  switch_flag( render_state_[SECONDARY], USE_COLOR_CONVERT, secondary_color_type_.get() == 2 );
  switch_flag( render_state_[TERTIARY], USE_COLOR_CONVERT, tertiary_color_type_.get() == 2 );

  switch_flag( render_state_[SCALAR], NORMALIZE_DATA, scalars_normalize_.get() );
  switch_flag( render_state_[VECTOR], NORMALIZE_DATA, vectors_normalize_.get() );
  switch_flag( render_state_[TENSOR], NORMALIZE_DATA, tensors_normalize_.get() );

  switch_flag( render_state_[SCALAR], SMALL_IS_DOT, scalars_small_is_dot_.get() );
  switch_flag( render_state_[VECTOR], SMALL_IS_DOT, vectors_small_is_dot_.get() );
  switch_flag( render_state_[TENSOR], SMALL_IS_DOT, tensors_small_is_dot_.get() );

  switch_flag( render_state_[VECTOR], BIDIRECTIONAL ,vectors_bidirectional_.get() );

  switch_flag( render_state_[SECONDARY], USE_ALPHA, secondary_alpha_.get() );
  switch_flag( render_state_[SECONDARY], USE_VALUE, secondary_value_.get() );
  switch_flag( render_state_[SECONDARY], USE_MAJOR_RADIUS,
	       secondary_display_type_.get() == "Major Radius" );
  switch_flag( render_state_[SECONDARY], USE_MINOR_RADIUS,
	       secondary_display_type_.get() == "Minor Radius" );
  switch_flag( render_state_[SECONDARY], USE_PITCH,
	       secondary_display_type_.get() == "Pitch" );

  switch_flag( render_state_[TERTIARY], USE_ALPHA, tertiary_alpha_.get() );
  switch_flag( render_state_[TERTIARY], USE_VALUE, tertiary_value_.get() );
  switch_flag( render_state_[TERTIARY], USE_MAJOR_RADIUS,
	       tertiary_display_type_.get() == "Major Radius" );
  switch_flag( render_state_[TERTIARY], USE_MINOR_RADIUS,
	       tertiary_display_type_.get() == "Minor Radius" );
  switch_flag( render_state_[TERTIARY], USE_PITCH,
	       tertiary_display_type_.get() == "Pitch" );

  std::string fname = clean_fieldname(gui_field_name_.get());
  if (fname != "" && fname[fname.size()-1] != ' ') { fname = fname + " "; }

  if(get_flag(render_state_[SCALAR], (IS_ON|HAS_DATA|DIRTY))) 
  {
    GeomHandle data_geometry =
      scalar_renderer_->render_data(pfld_handle,
				    sfld_handle,
				    tfld_handle,
				    scalars_display_type_.get(),
				    scalars_scale_.get(),
				    secondary_scale_.get(),
				    tertiary_scale_.get(),
            scalars_threshold_.get(),
				    scalars_resolution_.get(),
				    render_state_[SCALAR],
				    render_state_[SECONDARY],
				    render_state_[TERTIARY]);

    GeomHandle gmat = new GeomMaterial(data_geometry, def_material_);
    GeomHandle geom = new GeomSwitch(new GeomColorMap(gmat, cmap_handle));
    if (data_id_) ogeom->delObj(data_id_);
    data_id_ = ogeom->addObj(geom, fname +
			      (scalars_transparency_.get()?"Transparent Scalars":"Scalars"));
  }  
  else if(get_flag(render_state_[VECTOR], (IS_ON|HAS_DATA|DIRTY))) 
  {
    GeomHandle data_geometry =
      vector_renderer_->render_data(pfld_handle,
				    sfld_handle,
				    tfld_handle,
				    def_material_,
				    vectors_display_type_.get(),
				    vectors_scale_.get(),
				    secondary_scale_.get(),
				    tertiary_scale_.get(),
            vectors_threshold_.get(),
				    vectors_resolution_.get(),
				    render_state_[VECTOR],
				    render_state_[SECONDARY],
				    render_state_[TERTIARY]);

    GeomHandle gmat = new GeomMaterial(data_geometry, def_material_);
    GeomHandle geom = new GeomSwitch(new GeomColorMap(gmat, cmap_handle));
    const std::string vdname = vectors_transparency_.get() ?
      "Transparent Vectors":
      ((vectors_display_type_.get()=="Needles" ||
        vectors_display_type_.get()=="Comets" ) ?
       "Transparent Vectors":"Vectors");

    if (data_id_) ogeom->delObj(data_id_);
    data_id_ = ogeom->addObj(geom, fname + vdname);
  }
  else if(get_flag(render_state_[TENSOR], (IS_ON|HAS_DATA|DIRTY))) 
  {
    GeomHandle data_geometry =
      tensor_renderer_->render_data(pfld_handle,
				    sfld_handle,
				    tfld_handle,
				    tensors_display_type_.get(),
				    tensors_scale_.get(),
				    secondary_scale_.get(),
            0.0,
            tensors_threshold_.get(),
				    tensors_resolution_.get(),
				    tensors_emphasis_.get(),
				    render_state_[TENSOR],
				    render_state_[SECONDARY],
				    render_state_[TERTIARY]);

    GeomHandle gmat = new GeomMaterial(data_geometry, def_material_);
    GeomHandle geom = new GeomSwitch(new GeomColorMap(gmat, cmap_handle));
    const std::string tdname = tensors_transparency_.get() ?
      "Transparent Tensors":"Tensors";

    if (data_id_) ogeom->delObj(data_id_);
    data_id_ = ogeom->addObj(geom, fname + tdname);
  }

  if (get_flag(render_state_[TEXT], (IS_ON|DIRTY)))
  {
    GeomHandle text_geometry =
      text_renderer_->render_text(pfld_handle,
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
			      (text_backface_cull_.get()?"Culled Text Data":"Text Data"));
  }

  clear_flag( render_state_[SCALAR], DIRTY);
  clear_flag( render_state_[VECTOR], DIRTY);
  clear_flag( render_state_[TENSOR], DIRTY);
  clear_flag( render_state_[TEXT], DIRTY);
}

void
ShowFieldGlyphs::set_default_display_values() 
{
  double fact = cur_mesh_scale_factor_;
  vectors_scaleNV_.set(fact * 0.0735);
  tensors_scaleNV_.set(fact * 0.0735);
  scalars_scaleNV_.set(fact * 0.0735);

  vectors_thresholdNV_.set(fact * 1e-5);
  tensors_thresholdNV_.set(fact * 1e-5);
  scalars_thresholdNV_.set(fact * 1e-5);

  set_flag( render_state_[SCALAR], DIRTY);
  set_flag( render_state_[VECTOR], DIRTY);
  set_flag( render_state_[TENSOR], DIRTY);
}

void
ShowFieldGlyphs::maybe_execute(toggle_type_e dis_type)
{
  bool do_execute = false;

  interactive_mode_.reset();
  if (interactive_mode_.get() == "Interactive") 
  {
    switch(dis_type) 
    {
    case SCALAR :
      do_execute = scalars_on_.get();
      break;
    case VECTOR :
      do_execute = vectors_on_.get();
      break;
    case TENSOR :
      do_execute = tensors_on_.get();
      break;
    case TEXT :
      do_execute = text_on_.get();
      break;
    case ALL :
    default:
      do_execute = true;
      break;
    }
  }

  if (do_execute) 
  {
    want_to_execute();
  }
}


void
ShowFieldGlyphs::tcl_command(GuiArgs& args, void* userdata) 
{
  if(args.count() < 2)
  {
    args.error("ShowFieldGlyphs needs a minor command");
    return;
  }

  if (args[1] == "scalars_threshold" || args[1] == "scalars_scale" || args[1] == "scalars_resolution")
  {
    set_flag( render_state_[SCALAR], DIRTY);
    maybe_execute(SCALAR);

  } 
  else if (args[1] == "vectors_threshold" || args[1] == "vectors_scale" || args[1] == "vectors_resolution") 
  {
    set_flag( render_state_[VECTOR], DIRTY);
    maybe_execute(VECTOR);
  } 
  else if (args[1] == "tensors_threshold" || args[1] == "tensors_scale" ) 
  {
    set_flag( render_state_[TENSOR], DIRTY);
    maybe_execute(TENSOR);
  } 
  else if (args[1] == "tensors_resolution") 
  {
    if (tensors_display_type_.get() == "Ellipsoids" ||
        tensors_display_type_.get() == "Superquadrics")
    {
      set_flag( render_state_[TENSOR], DIRTY);
      maybe_execute(TENSOR);
    }
  } 
  else if (args[1] == "default_color_change") 
  {
    def_color_r_.reset();
    def_color_g_.reset();
    def_color_b_.reset();
    def_color_a_.reset();
    def_material_->diffuse =
      Color(def_color_r_.get(), def_color_g_.get(), def_color_b_.get());
    def_material_->transparency = def_color_a_.get();
    GeometryOPortHandle ogeom;
    get_oport_handle("Scene Graph",ogeom);
    if (ogeom.get_rep()) ogeom->flushViews();
  } 
  else if (args[1] == "text_color_change") 
  {
    text_color_r_.reset();
    text_color_g_.reset();
    text_color_b_.reset();
    text_material_->diffuse =
      Color(text_color_r_.get(), text_color_g_.get(), text_color_b_.get());
    GeometryOPortHandle ogeom;
    get_oport_handle("Scene Graph",ogeom);
    if (ogeom.get_rep()) ogeom->flushViews();

  } 
  else if (args[1] == "toggle_display_scalars")
  {
    // Toggle the GeomSwitch.
    scalars_on_.reset();
    if (scalars_on_.get())
    {
      set_flag( render_state_[SCALAR], DIRTY);
      maybe_execute(SCALAR);
    }
    else if (data_id_)
    {
      GeometryOPortHandle ogeom;
      get_oport_handle("Scene Graph",ogeom);
      if (ogeom.get_rep())
      {
        ogeom->delObj(data_id_);
        ogeom->flushViews();
      }
      data_id_ = 0;
    }
  } 
  else if (args[1] == "toggle_display_vectors")
  {
    // Toggle the GeomSwitch.
    vectors_on_.reset();
    if (vectors_on_.get())
    {
      set_flag( render_state_[VECTOR], DIRTY);
      maybe_execute(VECTOR);
    }
    else if (data_id_)
    {
      GeometryOPortHandle ogeom;
      get_oport_handle("Scene Graph",ogeom);
      if (ogeom.get_rep())
      {
        ogeom->delObj(data_id_);
        ogeom->flushViews();
      }
      data_id_ = 0;
    }
  } 
  else if (args[1] == "toggle_display_tensors")
  {
    // Toggle the GeomSwitch.
    tensors_on_.reset();
    if (tensors_on_.get())
    {
      set_flag( render_state_[TENSOR], DIRTY);
      maybe_execute(TENSOR);
    }
    else if (data_id_)
    {
      GeometryOPortHandle ogeom;
      get_oport_handle("Scene Graph",ogeom);
      if (ogeom.get_rep())
      {
        ogeom->delObj(data_id_);
        ogeom->flushViews();
      }
      data_id_ = 0;
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
      if (ogeom.get_rep())
      {
        ogeom->delObj(text_id_);
        ogeom->flushViews();
      }
      text_id_ = 0;
    }
  } 
  else if (args[1] == "rerender_all" ) 
  {
    set_flag( render_state_[SCALAR], DIRTY);
    set_flag( render_state_[VECTOR], DIRTY);
    set_flag( render_state_[TENSOR], DIRTY);
    maybe_execute(ALL);

  } 
  else if (args[1] == "rerender_scalars") 
  {
    set_flag( render_state_[SCALAR], DIRTY);
    maybe_execute(SCALAR);
  } 
  else if (args[1] == "rerender_vectors") 
  {
    set_flag( render_state_[VECTOR], DIRTY);
    maybe_execute(VECTOR);
  } 
  else if (args[1] == "rerender_tensors") 
  {
    set_flag( render_state_[TENSOR], DIRTY);
    maybe_execute(TENSOR);
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
ShowFieldGlyphs::post_read()
{
  char names[36][2][30] = { {"vectors-on", "vectors_on" },
		            {"normalize-vectors", "vectors_normalize" },
		            {"has_vector_data", "vectors_has_data" },
		            {"bidirectional", "vectors_bidirectional" },
		            {"tensors-on", "tensors_on" },
		            {"has_tensor_data", "tensors_has_data" },
		            {"tensors-emphasis", "tensors_emphasis" },
		            {"scalars-on", "scalars_on" },
		            {"scalars-transparency", "scalars_transparency" },
		            {"has_scalar_data", "scalars_has_data" },
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
		            {"data_display_type", "vectors_display_type" },
		            {"tensor_display_type", "tensors_display_type" },
		            {"scalar_display_type", "scalars_display_type" },
                {"field-name","field_name"  },
                {"field-name-override","field_name_override"  },
		            {"data-resolution", "vectors_resolution" },
                {"approx-div", "approx_div" },
                {"use-default-size", "use_default_size" } };

  // Get the module name
  const std::string modName = get_ctx()->getfullname() + "-";

  std::string val;

  for( unsigned int i=0; i<36; i++ )
  {
    // Get the current values for the old names
    if( TCLInterface::get(modName+names[i][0], val, get_ctx()) )
    {
      // Set the current values for the new names
      TCLInterface::set(modName+names[i][1], val, get_ctx());
    }
  }

  // Special cases for the old default color vars
  if( TCLInterface::get(modName+"scalars-usedefcolor", val, get_ctx()) )
  {
    if( val == std::string("1") )
      TCLInterface::set(modName+"scalars_color_type", "0", get_ctx());
  }

  if( TCLInterface::get(modName+"vectors-usedefcolor", val, get_ctx()) )
  {
    if( val == std::string("1") )
      TCLInterface::set(modName+"vectors_color_type", "0", get_ctx());
  }

  if( TCLInterface::get(modName+"tensors-usedefcolor", val, get_ctx()) )
  {
    if( val == std::string("1") )
      TCLInterface::set(modName+"tensors_color_type", "0", get_ctx());
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

DECLARE_MAKER(ShowFieldGlyphs)
} // End namespace SCIRun
