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


#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Module.h>


namespace BioPSE {

using namespace SCIRun;

class CalculateCurrentDensity : public Module {
public:
  CalculateCurrentDensity(GuiContext *context);
  virtual ~CalculateCurrentDensity() {}
  virtual void execute();
};

DECLARE_MAKER(CalculateCurrentDensity)

CalculateCurrentDensity::CalculateCurrentDensity(GuiContext *context)
  : Module("CalculateCurrentDensity", context, Filter, "Forward", "BioPSE")
{
}

void
CalculateCurrentDensity::execute()
{
  FieldHandle efieldH, sigmasH;
  
  get_input_handle("TetMesh EField", efieldH);
  get_input_handle("TetMesh Sigmas", sigmasH);

  if (efieldH->mesh().get_rep() != sigmasH->mesh().get_rep()) 
  {
    error("EField and Sigma Field need to have the same mesh.");
    return;
  }

  VField* efield = efieldH->vfield();
  VMesh*  emesh = efieldH->vmesh();
  VField* sfield = sigmasH->vfield();

  if (sfield->basis_order() != 0) 
  {
    error("Need sigmas at Elements");
    return;
  }
  
  if (efield->basis_order() != 0) 
  {
    error("Need efield at Elements");
    return;
  }

  if (!(efield->is_vector()))
  {
    error("Electric field needs to vector data");
    return;
  }
  
  std::vector<std::pair<std::string, Tensor> > conds;
  bool index_based = sfield->get_property("conductivity_table", conds); 
  

  FieldInformation fi(efieldH);
  fi.make_vector();
  FieldHandle output_field = CreateField(fi,efieldH->mesh());
  VField* ofield = output_field->vfield();

  VMesh::size_type num_elems = emesh->num_elems();

  Vector vec;
  Vector e;
  for (VMesh::Elem::index_type idx=0; idx<num_elems; idx++)
  {
    efield->get_value(e, idx);
    
    Tensor s;
    if (index_based) 
    {
      int sigma_idx;
      sfield->get_value(sigma_idx, idx);
      s=conds[sigma_idx].second;
    } 
    else 
    {
      sfield->get_value(s, idx);
    }

    // - sign added to vector to account for E = - Del V
    vec = Vector(-(s.mat_[0][0]*e.x()+s.mat_[0][1]*e.y()+s.mat_[0][2]*e.z()),
		 -(s.mat_[1][0]*e.x()+s.mat_[1][1]*e.y()+s.mat_[1][2]*e.z()),
		 -(s.mat_[2][0]*e.x()+s.mat_[2][1]*e.y()+s.mat_[2][2]*e.z()));

    ofield->set_value(vec,idx);
  }

  send_output_handle("Currents", output_field);
}

} // End namespace BioPSE
