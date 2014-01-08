/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
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
 *  GetFieldZeroCrossings.cc:
 *
 *  Written by:
 *   Allen Sanderson
 *   School of Computing
 *   University of Utah
 *   April 2005
 *
 */

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/FieldPort.h>

namespace FusionPSE {


using namespace SCIRun;


class GetFieldZeroCrossingsAlgo
{
  public:
    //! virtual interface. 
    FieldHandle execute(FieldHandle& src, int axis);
};


FieldHandle
GetFieldZeroCrossingsAlgo::execute(FieldHandle& ifield_h, int axis)
{
  // Create a new field based on same mesh
  
  FieldInformation fi(ifield_h);
  FieldHandle ofield_handle = CreateField(fi,ifield_h->mesh());
  
  VField *ofield = ofield_handle->vfield();
  VMesh*  omesh =  ofield->vmesh();

  // Zero all data
  ofield->set_all_values(0.0);

  // Nodes are working. Some mesh types need this.
  omesh->synchronize(Mesh::NODES_E);
  
  // Parallel synchonization goes faster as it will be done in parallel
  if (omesh->dimensionality() == 1) omesh->synchronize(Mesh::EDGES_E);
  if (omesh->dimensionality() == 2) omesh->synchronize(Mesh::FACES_E|Mesh::EDGES_E);
  if (omesh->dimensionality() == 3) omesh->synchronize(Mesh::CELLS_E|Mesh::FACES_E|Mesh::EDGES_E);

  VMesh::Node::array_type nodes;
  VMesh::Edge::iterator iEdgeItr;
  VMesh::Edge::iterator oEdgeItr;

  omesh->begin( iEdgeItr );
  omesh->end  ( oEdgeItr );

  Point p0, p1;

  // Iterate through each edge.
  while( iEdgeItr != oEdgeItr ) 
  {
    omesh->get_nodes(nodes, *iEdgeItr);

    // Get the point and value at the locations.
    omesh->get_center(p0, nodes[0]);
    omesh->get_center(p1, nodes[1]);

    for( unsigned int i=0; i<2; i++ ) 
    {
      Vector diff;

      // Get the derivative between the two points.
      if( i == 0 ) diff = (Vector) p0 - (Vector) p1;
      else         diff = (Vector) p1 - (Vector) p0;

      // Normalize to be -1 or +1
      if( diff.x() > 0 ) diff.x( 1 ); else if( diff.x() < 0 ) diff.x( -1 );
      if( diff.y() > 0 ) diff.y( 1 ); else if( diff.y() < 0 ) diff.y( -1 );
      if( diff.z() > 0 ) diff.z( 1 ); else if( diff.z() < 0 ) diff.z( -1 );

      // Sum the derivative sign at each point. If the sum is zero there
      // is a zero crossing.

      // Get the current sum.
      double val;
      ofield->get_value(val, nodes[i] );
      
      // Add the new edge to the sum.
      if( axis == 0 )      val += diff.x();
      else if( axis == 1 ) val += diff.y();
      else if( axis == 2 ) val += diff.z();

      // Resave the sum.
      ofield->set_value(val, nodes[i]);
    }

    ++iEdgeItr;
  }

  VMesh::Node::iterator inodeItr, inodeEnd;
  omesh->begin( inodeItr );
  omesh->end( inodeEnd );

  Point currPt;

  omesh->get_center(currPt, *inodeItr);

  Point lastPt = currPt;
  Point penultimatePt = lastPt;

  bool check = true;

  while (inodeItr != inodeEnd) 
  {
    // Get the current point.
    omesh->get_center(currPt, *inodeItr);

    double val;
    ofield->get_value(val,*inodeItr);
    
      // Find the positive zero crossing.
    if( check && 0.0 < currPt.z() &&
        penultimatePt.z() <= lastPt.z() &&
        currPt.z() < lastPt.z() ) 
    {
      val = 1.001;

      check = false;
    } 
    else 
    {
      val = 0.0;
    }

    if( currPt.z() < 0.0 && check == false ) check = true;
    
    ofield->set_value(val,*inodeItr);
    penultimatePt = lastPt;
    lastPt  = currPt;

    ++inodeItr;
  }

  return (ofield_handle);
}

class GetFieldZeroCrossings : public Module {
public:
  GetFieldZeroCrossings(GuiContext *context);

  virtual ~GetFieldZeroCrossings();

  virtual void execute();

private:
  GuiInt gui_axis_;

  FieldHandle field_output_handle_;
};


DECLARE_MAKER(GetFieldZeroCrossings)


GetFieldZeroCrossings::GetFieldZeroCrossings(GuiContext *context)
  : Module("GetFieldZeroCrossings", context, Filter, "Fields", "FusionPSE"),
    
    gui_axis_(context->subVar("axis"), 2)
{
}

GetFieldZeroCrossings::~GetFieldZeroCrossings()
{
}

void
GetFieldZeroCrossings::execute()
{
  FieldHandle field_input_handle;
  if( !get_input_handle( "Input Field", field_input_handle, true  ) ) return;

  if( field_input_handle->mesh()->topology_geometry() & Mesh::REGULAR ) {

    error( field_input_handle->get_type_description(Field::MESH_TD_E)->get_name() );
    error( "Only availible for topologically irregular data." );
    return;
  }

  if( field_input_handle->basis_order() != 1 ) {
    error( field_input_handle->get_type_description(Field::MESH_TD_E)->get_name() );
    error( "Currently only availible for node data." );
    return;
  }

  // If no data or a changed input field or axis recreate the mesh.
  if( inputs_changed_ ||
      !field_output_handle_.get_rep() ||
      gui_axis_.changed( true ) ) 
  {
    update_state(Executing);
    GetFieldZeroCrossingsAlgo algo;
    field_output_handle_ = algo.execute(field_input_handle,gui_axis_.get() );
  }

  // Send the data downstream
  send_output_handle( "Output Field", field_output_handle_, true );
}

} // End namespace SCIRun
