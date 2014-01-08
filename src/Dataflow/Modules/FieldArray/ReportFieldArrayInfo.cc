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


#include <Core/Util/StringUtil.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldArray.h>
#include <Core/Datatypes/VField.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Geometry/BBox.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>

namespace ModelCreation {

using namespace SCIRun;

class ReportFieldArrayInfo : public Module {
  public:
    ReportFieldArrayInfo(GuiContext* ctx);
    virtual ~ReportFieldArrayInfo() {}
    
    virtual void execute();
  private:
    GuiString gui_fldname_;
    GuiString gui_generation_;
    GuiString gui_typename_;
    GuiString gui_datamin_;
    GuiString gui_datamax_;
    GuiString gui_numnodes_;
    GuiString gui_numelems_;
    GuiString gui_numdata_;
    GuiString gui_dataat_;
    GuiString gui_cx_;
    GuiString gui_cy_;
    GuiString gui_cz_;
    GuiString gui_sizex_;
    GuiString gui_sizey_;
    GuiString gui_sizez_;
    GuiString gui_nodesx_;
    GuiString gui_nodesy_;
    GuiString gui_nodesz_;
    GuiInt    gui_selectfield_;
    GuiInt    gui_selectfield_min_;
    GuiInt    gui_selectfield_max_;
    
    double min_;
    double max_;
    Point  center_;
    Vector size_;
    Vector dimensions_;
    int    numelems_;
    int    numnodes_;
    int    numdata_;
    
    void clear_vals();
    void update_input_attributes(FieldHandle);
};


DECLARE_MAKER(ReportFieldArrayInfo)

ReportFieldArrayInfo::ReportFieldArrayInfo(GuiContext* ctx)
  : Module("ReportFieldArrayInfo", ctx, Sink, "FieldArray", "SCIRun"),
    gui_fldname_(get_ctx()->subVar("fldname", false),"---"),
    gui_generation_(get_ctx()->subVar("generation", false),"---"),
    gui_typename_(get_ctx()->subVar("typename", false),"---"),
    gui_datamin_(get_ctx()->subVar("datamin", false),"---"),
    gui_datamax_(get_ctx()->subVar("datamax", false),"---"),
    gui_numnodes_(get_ctx()->subVar("numnodes", false),"---"),
    gui_numelems_(get_ctx()->subVar("numelems", false),"---"),
    gui_numdata_(get_ctx()->subVar("numdata", false),"---"),
    gui_dataat_(get_ctx()->subVar("dataat", false),"---"),
    gui_cx_(get_ctx()->subVar("cx", false),"---"),
    gui_cy_(get_ctx()->subVar("cy", false),"---"),
    gui_cz_(get_ctx()->subVar("cz", false),"---"),
    gui_sizex_(get_ctx()->subVar("sizex", false),"---"),
    gui_sizey_(get_ctx()->subVar("sizey", false),"---"),
    gui_sizez_(get_ctx()->subVar("sizez", false),"---"),
    gui_nodesx_(get_ctx()->subVar("nodesx", false),"---"),
    gui_nodesy_(get_ctx()->subVar("nodesy", false),"---"),
    gui_nodesz_(get_ctx()->subVar("nodesz", false),"---"),
    gui_selectfield_(get_ctx()->subVar("selectfield", true),0),
    gui_selectfield_min_(get_ctx()->subVar("selectfield_min", true),0),
    gui_selectfield_max_(get_ctx()->subVar("selectfield_max", true),1),
    min_(0.0),
    max_(0.0),
    numelems_(0),
    numnodes_(0),
    numdata_(0)
{
}


void
ReportFieldArrayInfo::clear_vals()
{
  gui_fldname_.set("---");
  gui_generation_.set("---");
  gui_typename_.set("---");
  gui_datamin_.set("---");
  gui_datamax_.set("---");
  gui_numnodes_.set("---");
  gui_numelems_.set("---");
  gui_numdata_.set("---");
  gui_dataat_.set("---");
  gui_cx_.set("---");
  gui_cy_.set("---");
  gui_cz_.set("---");
  gui_sizex_.set("---");
  gui_sizey_.set("---");
  gui_sizez_.set("---");
  gui_nodesx_.set("---");
  gui_nodesy_.set("---");
  gui_nodesz_.set("---");}


void
ReportFieldArrayInfo::update_input_attributes(FieldHandle f)
{
  VField* vfield = f->vfield();
  VMesh*  vmesh  = f->vmesh();
  
  if (vfield)
  {
    // Get name of field
    std::string fldname;
    if (f->get_property("name",fldname))
    {
      gui_fldname_.set(fldname);
    }
    else
    {
      gui_fldname_.set("--- Name Not Assigned ---");
    }  
  
    // Generation
    gui_generation_.set(to_string(f->generation));   

    // Typename
    const std::string &tname = f->get_type_description()->get_name();
    gui_typename_.set(tname);  
  
    // Basis
    static const char *at_table[4] = { "Nodes", "Edges", "Faces", "Cells" };
    switch(f->basis_order())
    {
    case 3:
      gui_dataat_.set("Nodes (cubic basis)");
      break;
    case 2:
      gui_dataat_.set("Nodes (quadratic basis)");
      break;
    case 1:
      gui_dataat_.set("Nodes (linear basis)");
      break;
    case 0:
      gui_dataat_.set(at_table[f->vmesh()->dimensionality()] +
                      std::string(" (constant basis)"));
      break;
    case -1:
      gui_dataat_.set("None (nodata basis)");
      break;
    }

    Point center;
    Vector size;

    const BBox bbox = vmesh->get_bounding_box();
    if (bbox.valid())
    {
      size = bbox.diagonal();
      center = bbox.center();
      gui_cx_.set(to_string(center.x()));
      gui_cy_.set(to_string(center.y()));
      gui_cz_.set(to_string(center.z()));
      gui_sizex_.set(to_string(size.x()));
      gui_sizey_.set(to_string(size.y()));
      gui_sizez_.set(to_string(size.z()));
      
      size_ = size;
      center_ = center;
    }
    else
    {
      warning("Input Field is empty.");
      gui_cx_.set("--- N/A ---");
      gui_cy_.set("--- N/A ---");
      gui_cz_.set("--- N/A ---");
      gui_sizex_.set("--- N/A ---");
      gui_sizey_.set("--- N/A ---");
      gui_sizez_.set("--- N/A ---");

      size_ = Vector(0.0,0.0,0.0);
      center_ = Point(0.0,0.0,0.0);
    }

    if ((!vfield->is_tensor()) && (!vfield->is_vector()))
    {
      vfield->minmax(min_,max_);
      gui_datamin_.set(to_string(min_));
      gui_datamax_.set(to_string(max_));
    }
    else
    {
      gui_datamin_.set("--- N/A ---");
      gui_datamax_.set("--- N/A ---");
      
      min_ = 0.0;
      max_ = 0.0;
    }

    numdata_ = vfield->num_values();
    numnodes_ = vmesh->num_nodes();
    numelems_ = vmesh->num_elems();
    gui_numnodes_.set(to_string(numnodes_));
    gui_numelems_.set(to_string(numelems_));
    gui_numdata_.set(to_string(numdata_));

    VMesh::dimension_type dim;
    vmesh->get_dimensions(dim);
    dimensions_ = Vector(1.0,1.0,1.0);
    for (size_t p=0;p<dim.size();p++) dimensions_[p] = static_cast<double>(dim[p]);

    gui_nodesx_.set(to_string(dimensions_.x()));
    if (dim.size()>1) gui_nodesy_.set(to_string(dimensions_.y())); else gui_nodesy_.set("-- N/A --");
    if (dim.size()>2) gui_nodesz_.set(to_string(dimensions_.z())); else gui_nodesz_.set("-- N/A --");      

  }
}

void
ReportFieldArrayInfo::execute()
{  
  FieldArrayHandle fh;
  
  if(!(get_input_handle("Input FieldArray",fh,false)))
  {
    clear_vals();
    return;
  }

  if (inputs_changed_ || gui_selectfield_.changed())
  {
    std::vector<FieldHandle>& array = fh->array();
    gui_selectfield_min_.set(0);
    gui_selectfield_max_.set(array.size()-1);
    size_t select = static_cast<size_t>(gui_selectfield_.get());
    if (select >= array.size())
    {
      select = array.size()-1;
      gui_selectfield_.set(select);
    }
    
    TCLInterface::execute(get_id() + " update_range");
    
    update_input_attributes(array[select]);
  }
}

} // end SCIRun namespace
