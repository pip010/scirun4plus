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
 *  UnuSort
 *
 *  Written by:
 *   Allen R. Sanderson
 *   SCientific Computing & Imaging Institute
 *   University of Utah
 *   March 2007
 *
 */

#include <Core/Datatypes/NrrdData.h>

#include <Dataflow/Network/Ports/NrrdPort.h>
#include <Dataflow/Network/Module.h>

#include <algorithm>

namespace SCITeem {

using namespace SCIRun;

template< class NTYPE>
class UnuSortAlgoT
{
  public:
    class UnuSortAsc : public std::binary_function<size_t,size_t,bool>
    {
      public:
        UnuSortAsc(NrrdDataHandle nrrd, size_t index)
        {
          size_ = nrrd->nrrd_->axis[0].size;
          num_elems_ = nrrd->nrrd_->axis[1].size;
          index_ = index;
          data_ = reinterpret_cast<NTYPE*>(nrrd->nrrd_->data);
        }
        
        bool operator()(size_t i1, size_t i2)
        {
          NTYPE val1 = data_[size_*i1+index_];
          NTYPE val2 = data_[size_*i2+index_];
          return (val1 < val2);
        }

      private:
        size_t size_;
        size_t num_elems_;
        size_t index_;
        NTYPE* data_;
    };  
  
  
    UnuSortAlgoT(unsigned int index) : index_(index) {}
  
    bool run(NrrdDataHandle nrrd_input_handle, 
                     NrrdDataHandle& nrrd_output_handle);
    
  private:
    unsigned int index_;
};

 

template< class NTYPE>
bool UnuSortAlgoT< NTYPE>::
run(NrrdDataHandle nrrd_input_handle, NrrdDataHandle& nrrd_output_handle)
{
  nrrd_output_handle = nrrd_input_handle->clone();
  std::vector<size_t> order;

  size_t size = nrrd_input_handle->nrrd_->axis[0].size;
  size_t numelems = nrrd_input_handle->nrrd_->axis[1].size;

  order.resize(numelems);
  for (size_t p=0; p<numelems; p++) order[p] = p;

  std::sort(order.begin(),order.end(),UnuSortAsc(nrrd_input_handle,index_));

  NTYPE* data_in = reinterpret_cast<NTYPE*>(nrrd_input_handle->nrrd_->data);
  NTYPE* data_out = reinterpret_cast<NTYPE*>(nrrd_output_handle->nrrd_->data);

  for(size_t i=0; i<numelems; i++)
  {
    size_t oidx = order[i];
    for (size_t j=0; j<size;j++)
    {
      data_out[oidx*size+j] = data_in[i*size+j];
    }
  }

  return (true);
}




class UnuSort : public Module {

public:
  UnuSort(SCIRun::GuiContext *ctx);
  virtual ~UnuSort() {}
  virtual void execute();

protected:
  GuiInt gui_index_;
  GuiInt gui_min_;
  GuiInt gui_max_;
};


DECLARE_MAKER(UnuSort)

UnuSort::UnuSort(SCIRun::GuiContext *ctx) : 
  Module("UnuSort", ctx, Filter, "UnuNtoZ", "Teem"), 
  gui_index_(get_ctx()->subVar("index"), 0),
  gui_min_(get_ctx()->subVar("min"), 0),
  gui_max_(get_ctx()->subVar("max"), 1)
{
}

void 
UnuSort::execute()
{
  NrrdDataHandle nrrd_input_handle;

  get_input_handle("Nrrd", nrrd_input_handle,true);

  if( nrrd_input_handle->nrrd_->axis[0].size != (unsigned int)gui_max_.get() )
  {
    gui_min_.set( 0 );
    gui_max_.set( nrrd_input_handle->nrrd_->axis[0].size-1 );

    std::ostringstream str;
    str << get_id() << " set_min_max ";
    TCLInterface::execute(str.str().c_str());
  }

  // Must be a 2D nrrd for sorting.
  if( nrrd_input_handle->nrrd_->dim != 2 ) 
  {
    error( "Must be a 2D for sorting." );
    return;
  }

  // Check to see if the input field has changed.
  if( inputs_changed_ ||
      gui_index_.changed( true ) ||
      !oport_cached("Nrrd") )
  {
    // Force Teem to be locked befoer calling the Teem library
    NrrdGuard nrrd_guard;

    // Inform module that execution started
    update_state(Executing);
      
    NrrdDataHandle nrrd_output_handle;
    switch (nrrd_input_handle->nrrd_->type)
    {
      case nrrdTypeChar:
      {
        UnuSortAlgoT<char> algo(gui_index_.get());
        algo.run(nrrd_input_handle,nrrd_output_handle); break;
      }
      case nrrdTypeUChar:
      {
        UnuSortAlgoT<unsigned char> algo(gui_index_.get());
        algo.run(nrrd_input_handle,nrrd_output_handle); break;
      }
      case nrrdTypeShort:
      {
        UnuSortAlgoT<short> algo(gui_index_.get());
        algo.run(nrrd_input_handle,nrrd_output_handle); break;
      }
      case nrrdTypeUShort:
      {
        UnuSortAlgoT<unsigned short> algo(gui_index_.get());
        algo.run(nrrd_input_handle,nrrd_output_handle); break;
      }
      case nrrdTypeInt:
      {
        UnuSortAlgoT<int> algo(gui_index_.get());
        algo.run(nrrd_input_handle,nrrd_output_handle); break;
      }
      case nrrdTypeUInt:
      {
        UnuSortAlgoT<unsigned int> algo(gui_index_.get());
        algo.run(nrrd_input_handle,nrrd_output_handle); break;
      }
      case nrrdTypeLLong:
      {
        UnuSortAlgoT<long long> algo(gui_index_.get());
        algo.run(nrrd_input_handle,nrrd_output_handle); break;
      }
      case nrrdTypeULLong:
      {
        UnuSortAlgoT<unsigned long long> algo(gui_index_.get());
        algo.run(nrrd_input_handle,nrrd_output_handle); break;
      }
      case nrrdTypeFloat:
      {
        UnuSortAlgoT<float> algo(gui_index_.get());
        algo.run(nrrd_input_handle,nrrd_output_handle); break;
      }
      case nrrdTypeDouble:
      {
        UnuSortAlgoT<double> algo(gui_index_.get());
        algo.run(nrrd_input_handle,nrrd_output_handle); break;
      }
      default:
      {
        error("Unnown nrrd type encountered");
        return;
      }
    }

    send_output_handle("Nrrd", nrrd_output_handle);
  }
}

} // End namespace SCITeem

