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
 *  CreateGTCMeshConnections.cc:
 *
 *  Written by:
 *   Allen R. Sanderson
 *   SCI Institute
 *   University of Utah
 *   March 2008
 *
 */

#include <Dataflow/Network/Module.h>

#include <Core/Datatypes/DenseMatrix.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Ports/NrrdPort.h>

namespace FusionPSE {

using namespace SCIRun;

class SubsampleGTCMesh : public Module {
public:
  SubsampleGTCMesh(SCIRun::GuiContext *context);
  virtual ~SubsampleGTCMesh();
  virtual void execute();

private:
  vector<int> nPoloidalNodes;   // Number of nodes in each poloidal contour
  vector<int> poloidalIndex;    // Starting node index of each poloidal contour

  GuiInt gui_dim_i_;
  GuiInt gui_dim_j_;

  GuiInt gui_start_i_;
  GuiInt gui_start_j_;

  GuiInt gui_stop_i_;
  GuiInt gui_stop_j_;
};

DECLARE_MAKER(SubsampleGTCMesh)

SubsampleGTCMesh::SubsampleGTCMesh(SCIRun::GuiContext *context) : 
  Module("SubsampleGTCMesh", context, Source, "Fields", "FusionPSE"),
  gui_dim_i_(context->subVar("dim-i"), 2),
  gui_dim_j_(context->subVar("dim-j"), 2),
  
  gui_start_i_(context->subVar("start-i"), 0),
  gui_start_j_(context->subVar("start-j"), 0),

  gui_stop_i_(context->subVar("stop-i"), 1),
  gui_stop_j_(context->subVar("stop-j"), 1)
{
}

SubsampleGTCMesh::~SubsampleGTCMesh()
{
}

void 
SubsampleGTCMesh::execute()
{
  update_state(NeedData);

  NrrdDataHandle nin_handle;
  if (!get_input_handle("Input Nrrd", nin_handle));

   if( inputs_changed_ )
  {
    Nrrd *nin = nin_handle->nrrd_;

    int nNodes = nin->axis[1].size;  // Number of nodes in each poloidal plane
    int nPoloidalPlanes = nin->axis[2].size;  // Number of poloidal plane

    nPoloidalNodes.clear();   // Number of nodes in each poloidal contour
    poloidalIndex.clear();    // Starting node index of each poloidal contour

    float *nrrd_in_ptr = (float *) nin->data;

    float x = *nrrd_in_ptr; ++nrrd_in_ptr;
    float y = *nrrd_in_ptr; ++nrrd_in_ptr;
    float z = *nrrd_in_ptr; ++nrrd_in_ptr;

    Point basePt = Point(x, y, z);
    int cc = 1;  // Temporary counter of the nodes in each poloidal contour

    poloidalIndex.push_back(0);

    for (int i=1; i<nNodes; i++)
    {
      x = nrrd_in_ptr[0];
      y = nrrd_in_ptr[1];
      z = nrrd_in_ptr[2];
      nrrd_in_ptr += 3;
    
      Point tmpPt = Point(x, y, z);
      ++cc;

      // In each contour the first and last node point are the same.
      // Well almost - a bit of rounding.
      if( Vector( basePt-tmpPt ).length() < 1.0e-6 )
      {
	nPoloidalNodes.push_back(cc);
	
	// Last node ? if so quit.
	if( i == nNodes-1 )
	{
	  break;
	}
	// Restart for the next poloidal contour.
	else
	{
	  poloidalIndex.push_back(++i);
	
	  x = nrrd_in_ptr[0];
	  y = nrrd_in_ptr[1];
	  z = nrrd_in_ptr[2];
	  nrrd_in_ptr += 3;
	
	  basePt = Point(x, y, z);
	  cc = 1;
	}
      }
    }

    int nFluxSurfaces = nPoloidalNodes.size();  // Number of flux surfaces

    remark( string( "Found " ) +
	    to_string(nFluxSurfaces) +
	    string( " Flux Surfaces") );

    bool update_dims = false;

    //! Check to see if the gui dimensions are different than the field.
    if( gui_dim_i_.get() != static_cast<int>(nFluxSurfaces) ) 
    {
      gui_dim_i_.set( nFluxSurfaces );
      update_dims = true;
    }

    //! Check to see if the gui dimensions are different than the field.
    if( gui_dim_j_.get() != static_cast<int>(nPoloidalPlanes) ) 
    {
      gui_dim_j_.set( nPoloidalPlanes );
      update_dims = true;
    }

    //! If the gui dimensions are different than the field then update the gui.
    if( update_dims ) 
    {
      ostringstream str;
      str << get_id() << " set_size ";
      TCLInterface::execute(str.str().c_str());
      
      reset_vars();
    }
  }

  //! Get the optional matrix handle from the port. Note if a matrix is
  //! present it is sent down stream. Otherwise it will be created.
  MatrixHandle matrix_handle = 0;
  get_input_handle( "Input Matrix", matrix_handle, false );

  //! An input matrix is present so use the values in it to override
  //! the variables set in the gui.
  //! Column 0 start  index to subsample.
  //! Column 1 stop   index to subsample.
  //! Column 2 dimensions of the data.
  if( matrix_handle.get_rep() ) 
  {
    //! The matrix is optional. If present make sure it is a 2x3 matrix.
    //! The row indices is the axis index. The column is the data.
    if( (matrix_handle->nrows() != 2 || matrix_handle->ncols() != 3) ) 
    {
      error( "Input matrix is not a 2x3 matrix" );
      return;
    }

    //! Sanity check. Make sure the gui dimensions match the matrix
    //! dimensions.
    if( gui_dim_i_.get() != matrix_handle->get(0, 2) ||
        gui_dim_j_.get() != matrix_handle->get(1, 2) ) 
    {
      ostringstream str;
      str << "The dimensions of the matrix slicing do match the field. "
          << " Expected "
          << gui_dim_i_.get() << " "
          << gui_dim_j_.get()
          << " Got "
          << matrix_handle->get(0, 2) << " "
          << matrix_handle->get(1, 2);
            
      error( str.str() );
      return;
    }
    //! Check to see what index has been selected and if it matches
    //! the gui index.
    if( gui_start_i_.get() != (int) matrix_handle->get(0, 0) ||
        gui_start_j_.get() != (int) matrix_handle->get(1, 0) ||

        gui_stop_i_.get() != (int) matrix_handle->get(0, 1) ||
        gui_stop_j_.get() != (int) matrix_handle->get(1, 1) )
    {
      gui_start_i_.set( (int) matrix_handle->get(0, 0) );
      gui_start_j_.set( (int) matrix_handle->get(1, 0) );

      gui_stop_i_.set( (int) matrix_handle->get(0, 1) );
      gui_stop_j_.set( (int) matrix_handle->get(1, 1) );

      ostringstream str;
      str << get_id() << " update_index ";

      TCLInterface::execute(str.str().c_str());
	
      reset_vars();

      inputs_changed_ = true;
    }

  }

  //! If no data or an input change recreate the field. I.e Only
  //! execute when neeed.
  if( inputs_changed_ ||

      !oport_cached("Scalar") ||
      !oport_cached("Vector") ||

      gui_start_i_.changed( true ) ||
      gui_start_j_.changed( true ) ||
      
      gui_stop_i_.changed( true ) ||
      gui_stop_j_.changed( true ) )
  {

    remark( string("Range " + to_string( poloidalIndex[gui_start_i_.get()] ) +
		   " to " + to_string( poloidalIndex[gui_stop_i_.get()] +
				       nPoloidalNodes[gui_stop_i_.get()]-1 ) ) );

    DenseMatrix *scalar_matrix = new DenseMatrix(2,2);

    MatrixHandle matrix_handle = (MatrixHandle) scalar_matrix;
    
    scalar_matrix->put(0, 0, poloidalIndex[gui_start_i_.get()] );
    scalar_matrix->put(0, 1, poloidalIndex[gui_stop_i_.get()] +
		  nPoloidalNodes[gui_stop_i_.get()]-1 );
    
    scalar_matrix->put(1, 0, gui_start_j_.get() );
    scalar_matrix->put(1, 1, gui_stop_j_.get() );
    
    send_output_handle( "Scalar", matrix_handle );


    DenseMatrix *vector_matrix = new DenseMatrix(3,2);

    matrix_handle = (MatrixHandle) vector_matrix;
    
    vector_matrix->put(0, 0, 0 );
    vector_matrix->put(0, 1, 2 );

    vector_matrix->put(1, 0, poloidalIndex[gui_start_i_.get()] );
    vector_matrix->put(1, 1, poloidalIndex[gui_stop_i_.get()] +
		  nPoloidalNodes[gui_stop_i_.get()]-1 );
    
    vector_matrix->put(2, 0, gui_start_j_.get() );
    vector_matrix->put(2, 1, gui_stop_j_.get() );
    
    send_output_handle( "Vector", matrix_handle );


    DenseMatrix *base_matrix = new DenseMatrix(2,3);

    matrix_handle = (MatrixHandle) base_matrix;
    
    base_matrix->put(0, 0, gui_start_i_.get() );
    base_matrix->put(0, 1, gui_stop_i_.get() );
    base_matrix->put(0, 2, gui_dim_i_.get() );
    
    base_matrix->put(1, 0, gui_start_j_.get() );
    base_matrix->put(1, 1, gui_stop_j_.get() );
    base_matrix->put(1, 2, gui_dim_j_.get() );
    
    send_output_handle( "Base", matrix_handle );
  }
}


} // End namespace SCIFusionPSE
