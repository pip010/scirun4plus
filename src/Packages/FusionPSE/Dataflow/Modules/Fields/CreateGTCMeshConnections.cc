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
 *   Feburary 2008
 *
 */

#include <Dataflow/Network/Module.h>

#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/Network/Ports/NrrdPort.h>

namespace FusionPSE {

using namespace SCIRun;

class CreateGTCMeshConnections : public Module {
public:
  CreateGTCMeshConnections(SCIRun::GuiContext *context);
  virtual ~CreateGTCMeshConnections();
  virtual void execute();
};

DECLARE_MAKER(CreateGTCMeshConnections)

CreateGTCMeshConnections::CreateGTCMeshConnections(SCIRun::GuiContext *context) : 
  Module("CreateGTCMeshConnections", context, Source, "Fields", "FusionPSE")
{
}

CreateGTCMeshConnections::~CreateGTCMeshConnections()
{
}

void 
CreateGTCMeshConnections::execute()
{
  update_state(NeedData);

  NrrdDataHandle nin_handle;
  if (!get_input_handle("InputNrrd", nin_handle));



  if( inputs_changed_ ||

      !oport_cached("OutputNrrd") )
  {
    Nrrd *nin = nin_handle->nrrd_;

    int nNodes = nin->axis[1].size;  // Number of nodes in each poloidal plane
    int nPoloidalPlanes = nin->axis[2].size;  // Number of poloidal plane

    vector<int> nPoloidalNodes;   // Number of nodes in each poloidal contour
    vector<int> poloidalIndex;    // Starting node index of each poloidal contour

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

    if( poloidalIndex.size() != nPoloidalNodes.size() )
    {
      error( "Could not find matching start and stop points.");
      return;
    }

    int nFluxSurfaces = nPoloidalNodes.size();  // Number of flux surfaces

    // Index of the closest node on the neighboring contour.
    int *neighborIndex = new int[nNodes];

    // Work from the outside to the inside because there are more nodes
    // on the outside. As such, neighbors will get used multiple times
    // thus allowing for degenerate connections to be found.
    for (int k=nFluxSurfaces-1; k>0; k--)
    {
      for (int j=0; j<nPoloidalNodes[k]-1; j++ )
      {
	int l = poloidalIndex[k] + j;  // Index of the working node.

	nrrd_in_ptr = (float*) nin->data + l*3;

	x = nrrd_in_ptr[0];
	y = nrrd_in_ptr[1];
	z = nrrd_in_ptr[2];
	nrrd_in_ptr += 3;
    
	Point basePt = Point(x, y, z);

	double rmin = 1.0e8;

	// Find the node on the adjacent contour that is the closest to
	// the working node. Brute force search.

	// Never search the last node because it is the same as the
	// first node.
	for (int i=0; i<nPoloidalNodes[k-1]-1; i++ )
	{
	  int m = poloidalIndex[k-1] + i;

	  nrrd_in_ptr = (float*) nin->data + m*3;

	  x = nrrd_in_ptr[0];
	  y = nrrd_in_ptr[1];
	  z = nrrd_in_ptr[2];
	  nrrd_in_ptr += 3;
	
	  Point tmpPt = Point(x, y, z);
	
	  double rdiff = Vector( basePt-tmpPt ).length();

	  if (rdiff <= rmin )
	  {
	    neighborIndex[l] = m;
	    rmin = rdiff;
	  }
	}
      }
    }

    vector<int> connections;

    // Work from the outside to the inside because there are more nodes
    // on the outside. As such, neighbors will get used multiple times
    // thus allowing for degenerate connections to be found.
    for (int k=nFluxSurfaces-1; k>0; k--)
    {
      int nDegenerate = 0;  // Count for the number of degenerate connections.

      for (int j=0; j<nPoloidalNodes[k]-1; j++ )
      {
	int l = poloidalIndex[k] + j;
	int lp1 = (l + 1);
      
	// Never use the last node cause it is the same as the first
	// node.
	if (lp1 == poloidalIndex[k] + nPoloidalNodes[k] - 1)
	  lp1 = poloidalIndex[k];

	connections.push_back( lp1 );
	connections.push_back( neighborIndex[lp1] );
	connections.push_back( l );

	if( neighborIndex[l] != neighborIndex[lp1] )
	{
	  connections.push_back( neighborIndex[lp1] );
	  connections.push_back( neighborIndex[l] );
	  connections.push_back( l );
	}
	else
	{
	  ++nDegenerate;
	}
      }

      // The number of degenerate connections should always equal the
      // difference in the number of nodes between two countors.
      ASSERT( nDegenerate == nPoloidalNodes[k] - nPoloidalNodes[k-1] );
    }

    delete[] neighborIndex;
    
    // Total number of element in one poloidal slice.
    int nElements = connections.size() / 3;
    
    Nrrd *nout = nrrdNew();
    size_t size[NRRD_DIM_MAX];
    
    // If only one slice then create connections for trisurf.
    if( nPoloidalPlanes == 1 )
    {
      // Create a local array of axis sizes (3xN), so we can allocate
      // the output Nrrd
      size[0] = 3;
      size[1] = nElements;
    
      // Allocate the nrrd's data, set the size of each axis
      nrrdAlloc_nva(nout, nrrdTypeInt, 2, size);
    
      int *nrrd_out_ptr = (int *) nout->data;
    
      for (int i=0; i<nElements; i++ )
      {
	int index = i * 3;
      
	nrrd_out_ptr[0] = connections[index + 0];
	nrrd_out_ptr[1] = connections[index + 1];
	nrrd_out_ptr[2] = connections[index + 2];

	nrrd_out_ptr += 3;
      }
    }
    // Multiple slices so create connections for prisms.
    else
    {
      // Create a local array of axis sizes (6xN), so we can allocate
      // the output Nrrd
      size[0] = 6;
      size[1] = (nPoloidalPlanes-1) * nElements;
    
      // Allocate the nrrd's data, set the size of each axis
      nrrdAlloc_nva(nout, nrrdTypeInt, 2, size);
    
      int *nrrd_out_ptr = (int *) nout->data;

      //  Connect along the toriodal direction.
      for (int n=0; n<nPoloidalPlanes-1; n++ )
      {
	int offset = n * nNodes;
	int offset1 = (n+1) * nNodes;
      
	//  Connect along the poliodal plane.
	for (int i=0; i<nElements; i++ )
	{
	  int index = i * 3;
      
	  nrrd_out_ptr[0] = connections[index + 2] + offset;
	  nrrd_out_ptr[1] = connections[index + 1] + offset;
	  nrrd_out_ptr[2] = connections[index + 0] + offset;

	  nrrd_out_ptr[3] = connections[index + 2] + offset1;
	  nrrd_out_ptr[4] = connections[index + 1] + offset1;
	  nrrd_out_ptr[5] = connections[index + 0] + offset1;

	  nrrd_out_ptr += 6;
	}
      }
    }

    // Set axis info for (new) axis 0
    nout->axis[0].kind=nrrdKindUnknown;
    nout->axis[0].spacing=AIR_NAN;
    nout->axis[0].min=AIR_NAN;
    nout->axis[0].max=AIR_NAN;

    // Create SCIRun data structure wrapped around nout
    NrrdDataHandle nout_handle(new NrrdData(nout));

    // Copy the properties
    nout_handle->copy_properties(nin_handle.get_rep());

    send_output_handle("OutputNrrd", nout_handle);
  }
}

} // End namespace SCIFusionPSE
