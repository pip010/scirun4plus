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
 *  UnuCrop
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   January 2000
 *
 */

#include <Dataflow/Network/Module.h>

#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/Network/Ports/NrrdPort.h>
#include <Core/Util/StringUtil.h>

#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <sstream>
#include <iostream>
using std::endl;
#include <stdio.h>

namespace SCITeem {
using namespace SCIRun;

class UnuCrop : public Module {

  public:
    UnuCrop(SCIRun::GuiContext *ctx);
    virtual ~UnuCrop() {}

    virtual void execute();
    virtual void tcl_command(GuiArgs&, void*);

    int valid_data(std::string *minS, std::string *maxS, int *min, int *max, Nrrd *nrrd);
    int getint(const char *str, int *n);

    void add_axis(int i, int max_index);
    void remove_axis(int i);

    int parse(const NrrdDataHandle& handle, const std::string& val, const int axis);

  private:
    std::vector<GuiString*> gui_mins_;
    std::vector<GuiString*> gui_maxs_;
    std::vector<GuiInt*> gui_abs_maxs_;
    GuiInt          gui_num_axes_;
    GuiInt          gui_uis_;
    GuiInt          gui_reset_data_;
    GuiInt          gui_digits_only_;
    std::vector<int>     lastmin_;
    std::vector<int>     lastmax_;
    bool            first_time_;
};

} // End namespace SCITeem
using namespace SCITeem;


DECLARE_MAKER(UnuCrop)

UnuCrop::UnuCrop(SCIRun::GuiContext *ctx) : 
  Module("UnuCrop", ctx, Filter, "UnuAtoM", "Teem"), 
  gui_num_axes_(get_ctx()->subVar("num-axes")),
  gui_uis_(get_ctx()->subVar("uis")),
  gui_reset_data_(get_ctx()->subVar("reset_data")),
  gui_digits_only_(get_ctx()->subVar("digits_only")),
  first_time_(true)
{
  // this will get overwritten when tcl side initializes, but 
  // until then make sure it is initialized.
  gui_num_axes_.set(0);
 
  lastmin_.resize(4, -1);
  lastmax_.resize(4, -1);  

  for (int a = 0; a < 4; a++) 
  {
    std::ostringstream str;
    str << "minAxis" << a;
    gui_mins_.push_back(new GuiString(get_ctx()->subVar(str.str())));
    std::ostringstream str1;
    str1 << "maxAxis" << a;
    gui_maxs_.push_back(new GuiString(get_ctx()->subVar(str1.str())));
    std::ostringstream str2;
    str2 << "absmaxAxis" << a;
    gui_abs_maxs_.push_back(new GuiInt(get_ctx()->subVar(str2.str())));
  }
}

void 
UnuCrop::execute()
{
  NrrdDataHandle input_nrrd_handle;
  get_input_handle("Nrrd", input_nrrd_handle,true);

  // Create any resample axes that might have been saved
  if (first_time_)
  {
    first_time_ = false;

    for(int i=4; i<gui_uis_.get(); i++)
    {
      std::ostringstream str, str2, str3, str4;
      str << "minAxis" << i;
      str2 << "maxAxis" << i;
      str3 << "absmaxAxis" << i;
      str4 << i;
      gui_mins_.push_back(new GuiString(get_ctx()->subVar(str.str())));
      gui_maxs_.push_back(new GuiString(get_ctx()->subVar(str2.str())));
      gui_abs_maxs_.push_back(new GuiInt(get_ctx()->subVar(str3.str())));

      gui_mins_[i]->reset();
      gui_maxs_[i]->reset();
      gui_abs_maxs_[i]->reset();

      lastmin_.push_back(parse(input_nrrd_handle, gui_mins_[i]->get(),i));
      lastmax_.push_back(parse(input_nrrd_handle, gui_maxs_[i]->get(),i));  

      TCLInterface::execute(get_id().c_str() +
			 std::string(" make_min_max " + str4.str()));
    }
  }

  gui_num_axes_.set(input_nrrd_handle->nrrd_->dim);
  gui_num_axes_.reset();

  // Add uis needed
  while((unsigned int) gui_uis_.get() < input_nrrd_handle->nrrd_->dim)
  {
    add_axis(gui_uis_.get(),
	     input_nrrd_handle->nrrd_->axis[gui_uis_.get()].size-1);
  }
  
  // Remove unused uis
  while((unsigned int) gui_uis_.get() > input_nrrd_handle->nrrd_->dim)
  {
    remove_axis(gui_uis_.get() - 1);
  }

  if (inputs_changed_)
  {
    for (int a=0; a<gui_num_axes_.get(); a++)
    {
      int max = input_nrrd_handle->nrrd_->axis[a].size - 1;
      gui_abs_maxs_[a]->set(input_nrrd_handle->nrrd_->axis[a].size-1);
      gui_abs_maxs_[a]->reset();

      if (parse(input_nrrd_handle, gui_maxs_[a]->get(),a) > max)
      {
        if (gui_digits_only_.get() == 1)
        {
          warning("Out of bounds, resetting axis min/max");
          gui_mins_[a]->set(to_string(0));
          gui_mins_[a]->reset();
          gui_maxs_[a]->set(to_string(max));
          gui_maxs_[a]->reset();
          lastmin_[a] = 0;
          lastmax_[a] = max;
        }
        else
        {
          warning("Out of bounds, setting axis min/max to 0 to M");
          gui_mins_[a]->set(to_string(0));
          gui_mins_[a]->reset();
          gui_maxs_[a]->set("M");
          gui_maxs_[a]->reset();
          lastmin_[a] = 0;
          lastmax_[a] = max;
        }
      }
    }

    TCLInterface::execute(get_id() + " update_sizes ");
  }

  if (inputs_changed_ && !first_time_ && gui_reset_data_.get() == 1)
  {
    std::ostringstream str;
    str << get_id() << " reset_vals" << endl; 
    TCLInterface::execute(str.str());  
  }
  
  if (gui_num_axes_.get() == 0)
  {
    warning("Trying to crop a nrrd with no axes" );
    return;
  }

  MatrixHandle matrixH;

  // If a matrix present use those values.
  if (get_input_handle("Current Index", matrixH, false))
  {
    if( gui_num_axes_.get() != matrixH.get_rep()->nrows() ||
        matrixH.get_rep()->ncols() != 2 ) 
    {
      error("Input matrix size does not match nrrd dimensions." );
      return;
    }

    for (int a=0; a<gui_num_axes_.get(); a++) 
    {
      int min, max;

      min = (int) matrixH.get_rep()->get(a, 0);
      max = (int) matrixH.get_rep()->get(a, 1);

      gui_mins_[a]->set(to_string(min));
      gui_mins_[a]->reset();
      gui_maxs_[a]->set(to_string(max));
      gui_mins_[a]->reset();
    }
  }

  // See if any of the sizes have changed.
  for (int i=0; i<gui_num_axes_.get(); i++)
  {
    gui_mins_[i]->reset();
    int min = parse(input_nrrd_handle, gui_mins_[i]->get(),i);
    if (lastmin_[i] != min)
    {
      inputs_changed_ = true;
      lastmin_[i] = min;
    }
    gui_maxs_[i]->reset();
    int max = parse(input_nrrd_handle, gui_maxs_[i]->get(),i);
    if (lastmax_[i] != max)
    {
      inputs_changed_ = true;
    }
  }

  if( inputs_changed_ ||
      !oport_cached("Nrrd") ||
      !oport_cached("Selected Index"))
  {
    // Force Teem to be locked before calling the Teem library
    NrrdGuard nrrd_guard;

    // Inform module that execution started
    update_state(Executing);

    Nrrd *nin = input_nrrd_handle->nrrd_;
    Nrrd *nout = nrrdNew();

    std::vector<size_t> min(gui_num_axes_.get());
    std::vector<size_t> max(gui_num_axes_.get());

    for(int i=0; i< gui_num_axes_.get(); i++)
    {
      gui_mins_[i]->reset();
      gui_maxs_[i]->reset();
      min[i] = parse(input_nrrd_handle, gui_mins_[i]->get(),i);
      max[i] = parse(input_nrrd_handle, gui_maxs_[i]->get(),i);
      
      if (nrrdKindSize(nin->axis[i].kind) > 1 &&
          (min[i] != 0 || max[i] != (size_t) gui_abs_maxs_[i]->get()))
      {
        warning("Trying to crop axis " + to_string(i) +
          " which does not have a kind of nrrdKindDomain or nrrdKindUnknown");
      }
    }

    if (nrrdCrop(nout, nin, &min[0], &max[0]))
    {
      char *err = biffGetDone(NRRD);
      error(std::string("Trouble cropping: ") + err + "\nOutputting input Nrrd");
      remark("  Input Nrrd: nin->dim=" + to_string(nin->dim));
      free(err);

      for (unsigned int a=0; a<nin->dim; a++)
      {
        gui_mins_[a]->set(to_string(0));
        gui_maxs_[a]->set(to_string(nin->axis[a].size-1));
        gui_mins_[a]->reset();
        gui_maxs_[a]->reset();

        lastmin_[a] = min[a];
        lastmax_[a] = max[a];
      }

      send_output_handle("Nrrd", input_nrrd_handle);
    }
    else
    {
      for(unsigned int a=0; a<nin->dim; a++)
      {
        lastmin_[a] = min[a];
        lastmax_[a] = max[a];
      }

      for( unsigned int i=0; i<nin->dim; i++ )
      {
        nout->axis[i].kind  = nin->axis[i].kind;
        nout->axis[i].label = airStrdup(nin->axis[i].label);
      }

      if( (nout->axis[0].kind == nrrdKind3Vector     &&
            nout->axis[0].size != 3) ||
          (nout->axis[0].kind == nrrdKind3DSymMatrix &&
           nout->axis[0].size != 6) )
      {
        nout->axis[0].kind = nrrdKindDomain;
      }
      
      NrrdData *nrrd = new NrrdData(nout);

      NrrdDataHandle output_nrrd_handle = NrrdDataHandle(nrrd);

      // Copy the properties, kinds, and labels.
      output_nrrd_handle->copy_properties(input_nrrd_handle.get_rep());

      send_output_handle("Nrrd", output_nrrd_handle);
    }
    
    DenseMatrix *indexMat = new DenseMatrix( gui_num_axes_.get(), 2 );
    MatrixHandle output_matrix_handle = MatrixHandle(indexMat);

    for(int i=0; i< gui_num_axes_.get(); i++)
    {
      indexMat->put(i, 0, (double) min[i]);
      indexMat->put(i, 1, (double) max[i]);
    }

    send_output_handle("Selected Index", output_matrix_handle);
  }
}


void 
UnuCrop::tcl_command(GuiArgs& args, void* userdata)
{
  if(args.count() < 2){
    args.error("UnuCrop needs a minor command");
    return;
  }

  if( args[1] == "add_axis" ) 
  {
    add_axis(gui_uis_.get(), 0);
  }
  else if( args[1] == "remove_axis" ) 
  {
    remove_axis(gui_uis_.get() - 1);
  }
  else 
  {
    Module::tcl_command(args, userdata);
  }
}


void 
UnuCrop::add_axis(int i, int max_index)
{
  std::ostringstream str1, str2, str3, str4;
  str1 << "minAxis" << i;
  str2 << "maxAxis" << i;
  str3 << "absmaxAxis" << i;
  str4 << i;

  gui_mins_.push_back(new GuiString(get_ctx()->subVar(str1.str())));
  gui_mins_[i]->set(0);
  gui_mins_[i]->reset();
  
  gui_maxs_.push_back(new GuiString(get_ctx()->subVar(str2.str())));
  if (gui_digits_only_.get() == 1)
  {
    gui_maxs_[i]->set(to_string( max_index ));
    gui_maxs_[i]->reset();
  }
  else
  {
    gui_maxs_[i]->set("M");
    gui_maxs_[i]->reset();
  }
  
  gui_abs_maxs_.push_back(new GuiInt(get_ctx()->subVar(str3.str())));
  gui_abs_maxs_[i]->set(max_index);
  gui_abs_maxs_[i]->reset();
  
  lastmin_.push_back(0);
  lastmax_.push_back(max_index); 
  
  TCLInterface::execute(get_id() + " make_min_max " + str4.str());

  gui_uis_.set(i);
  gui_uis_.reset();
}


void 
UnuCrop::remove_axis(int i)
{
  std::vector<GuiString*>::iterator iter1 = gui_mins_.end();
  std::vector<GuiString*>::iterator iter2 = gui_maxs_.end();
  std::vector<GuiInt*>::iterator iter3 = gui_abs_maxs_.end();
  std::vector<int>::iterator iter4 = lastmin_.end();
  std::vector<int>::iterator iter5 = lastmax_.end();

  gui_mins_.erase(iter1, iter1);
  gui_maxs_.erase(iter2, iter2);
  gui_abs_maxs_.erase(iter3, iter3);
    
  lastmin_.erase(iter4, iter4);
  lastmax_.erase(iter5, iter5);
    
  TCLInterface::execute(get_id() + " clear_axis " + to_string(i));

  gui_uis_.set(i);
  gui_uis_.reset();
}


int 
UnuCrop::parse(const NrrdDataHandle& nH, const std::string& val, const int a) {
  // parse string val which must be in the form 
  // M, M+int, M-int, m+int, int
  
  // remove white spaces
  std::string new_val;
  for(size_t i = 0; i<val.length(); i++) 
  {
    if (val[i] != ' ') new_val += val[i];
  }
  
  int val_length = new_val.length();
  
  bool has_base = false;
  int base = 0;
  
  char op = '+';
  
  int int_result = 0;
  int start = 0;
  
  if (val_length == 0) 
  {
    error("Error in UnuCrop::parse String length 0.");
    return 0;
  }
  
  if (new_val[0] == 'M') 
  {
    has_base = true;
    base = nH->nrrd_->axis[a].size - 1;
  } 
  else if (new_val[0] == 'm') 
  { 
    has_base = true;
    base = parse(nH, gui_mins_[a]->get(), a);
  }
  
  if (has_base && val_length == 1) 
  {
    return base;
  }
  
  if (has_base)  
  {
    start = 2;
    if (new_val[1] == '+') 
    {
      op = '+';
    } 
    else if (new_val[1] == '-') 
    {
      op = '-';
    } 
    else 
    {
      error("Error UnuCrop::parse Must have +/- operation when using M or m with integers");
      return 0;
    }
  }

  if (!string_to_int(new_val.substr(start,val_length), int_result)) {
    error("Error UnuCrop::could not convert to integer");
    return 0;
  }

  if (has_base) 
  {
    if (op == '+') 
    {
      int_result = base + int_result;
    } 
    else 
    {
      int_result = base - int_result;
    }
  }

  return int_result;
}
