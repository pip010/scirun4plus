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
 *  SetupBEMatrix.cc: 
 *
 *  Written by:
 *   Andrew Keely - Northeastern University
 *   Michael Callahan - Department of Computer Science - University of Utah
 *   July 2006
 *
 */

#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/ColumnMatrix.h>

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Module.h>

namespace BioPSE {

using namespace SCIRun;

class BuildBEMSrcToInfTransferMatrix : public Module {
  public:
    BuildBEMSrcToInfTransferMatrix(GuiContext* ctx);
    virtual ~BuildBEMSrcToInfTransferMatrix() {}
    virtual void execute();
};

DECLARE_MAKER(BuildBEMSrcToInfTransferMatrix)

BuildBEMSrcToInfTransferMatrix::BuildBEMSrcToInfTransferMatrix(GuiContext *context):
  Module("BuildBEMSrcToInfTransferMatrix", context, Source, "Forward", "BioPSE")
{
}

void
BuildBEMSrcToInfTransferMatrix::execute()
{
  FieldHandle surface, dipoles;
  
  get_input_handle("Surface", surface, true);
  get_input_handle("Dipoles", dipoles, true);
 
  // TODO: Check for point cloud of vectors in dipoles!
 
  FieldInformation sfi(surface);
  FieldInformation dfi(dipoles);
 
  if (!(dfi.is_vector()))
  {
    error("Dipoles needs to have vector data");
    return;
  }
  
  VField*  sfield = surface->vfield();
  VField*  dfield = dipoles->vfield();
  
  VField::size_type num_values = sfield->num_values();
    
  MatrixHandle output = new ColumnMatrix(num_values);
  output->zero();
  double* output_data = output->get_data_pointer();
  
  VField::size_type num_dipoles = dfield->num_values();  

  Point p;
  Vector v;  
  for (VField::index_type idx = 0;idx<num_dipoles;idx++)
  {
  
    dfield->get_center(p,idx);
    dfield->get_value(v, idx);
    
    for (VField::index_type sidx = 0;sidx<num_values;sidx++)
    {
      Point surface_point;
      sfield->get_center(surface_point, sidx);

      double result;
      // Run some function of p, v, surface_point, put in result;
      
      result = (v.x() * (surface_point.x()-p.x()) +
                v.y() * (surface_point.y()-p.y()) +
                v.z() * (surface_point.z()-p.z())) /
               (4 * M_PI * pow(pow(surface_point.x()-p.x(),2) +
                           pow(surface_point.y()-p.y(),2) +
                           pow(surface_point.z()-p.z(),2), 3/2));

      output_data[sidx] = output_data[sidx]+result;
    }
  }
  
  send_output_handle("Surface Potentials", output);
}

} // end namespace BioPSE
