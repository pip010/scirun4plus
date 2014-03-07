/*
 * The MIT License
 * 
 * Copyright (c) 2009 Scientific Computing and Imaging Institute,
 * University of Utah.
 * 
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldArray.h>

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Ports/FieldArrayPort.h>
#include <Dataflow/Network/Module.h>

#include <algorithm>

namespace SCIRun {

using namespace SCIRun;

class InsertFieldIntoFieldArray : public Module {
  public:
    InsertFieldIntoFieldArray(GuiContext*);
    virtual ~InsertFieldIntoFieldArray() {}
    virtual void execute();

  private:
    GuiInt index_min_;  
    GuiInt index_max_;  
    GuiInt index_;  
};


DECLARE_MAKER(InsertFieldIntoFieldArray)

InsertFieldIntoFieldArray::InsertFieldIntoFieldArray(GuiContext* ctx) :
  Module("InsertFieldIntoFieldArray", ctx, Source, "FieldArray", "SCIRun"),
  index_min_(get_ctx()->subVar("index-min"),0),
  index_max_(get_ctx()->subVar("index-max"),1),
  index_(get_ctx()->subVar("index"),0)
{
}

void
InsertFieldIntoFieldArray::execute()
{
  FieldArrayHandle fa;
  get_input_handle("FieldArray",fa,true);
  
  MatrixHandle index;
  get_input_handle("Index",index,false);  
  
  FieldHandle field;
  get_input_handle("Field",field,true);

  if (inputs_changed_ || index_.changed() || !oport_cached("FieldArray")  ||
      !oport_cached("Index"))
  {
    if (index.get_rep())
    {
      index_.set(static_cast<int>(index->get(0,0)));
      index_.reset();
    }
    
    std::vector<FieldHandle>& array = fa->array();
    
    index_min_.set(0);
    index_max_.set(array.size()-1);  
    TCLInterface::execute(get_id() + " update_range");
    
    int idx = index_.get();
    if (idx < 0) 
    {
      error("Index is negative");
      return;
    }
    else if (idx >= static_cast<int>(array.size()))
    {
      error("Index is out of bound");
      return;
    }

    FieldArrayHandle output = new FieldArray;
    
    std::vector<FieldHandle>& oarray = output->array();
    oarray.resize(array.size());
    std::copy(array.begin(),array.end(),oarray.begin());
    
    oarray[idx] = field;
    send_output_handle("FieldArray",output,true);

    index = new DenseMatrix(static_cast<double>(idx));
    send_output_handle("Index",index,true);
  }

}


} // End namespace SCIRun


