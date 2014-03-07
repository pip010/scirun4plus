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
 *  UnuJhisto.cc 
 *  the opposite of "crop".
 *
 *  Written by:
 *   Darby Van Uitert
 *   April 2004
 *
 */


#include <Core/Util/StringUtil.h>

#include <Dataflow/Network/Ports/NrrdPort.h>
#include <Dataflow/Network/Module.h>

namespace SCITeem {

using namespace SCIRun;

class UnuJhisto : public Module {
  public:
    UnuJhisto(GuiContext*);
    virtual ~UnuJhisto() {}
    virtual void execute();

    GuiString       bins_;
    GuiString       mins_;
    GuiString       maxs_;
    GuiString       type_;
};


DECLARE_MAKER(UnuJhisto)
UnuJhisto::UnuJhisto(GuiContext* ctx)
  : Module("UnuJhisto", ctx, Source, "UnuAtoM", "Teem"),
    bins_(get_ctx()->subVar("bins")), mins_(get_ctx()->subVar("mins")), 
    maxs_(get_ctx()->subVar("maxs")), type_(get_ctx()->subVar("type"))
{
}


void
UnuJhisto::execute()
{
  update_state(NeedData);

  NrrdDataHandle nrrd_handle1;
  get_input_handle("InputNrrd1", nrrd_handle1,true);

  NrrdDataHandle weight_handle;
  Nrrd* weight = 0;
  if(get_input_handle("WeightNrrd", weight_handle, false))
  {
    weight = weight_handle->nrrd_;
  }

  std::vector<NrrdDataHandle> nrrd_handle2;
  get_dynamic_input_handles("InputNrrd2",nrrd_handle2,true);
  
  std::vector<NrrdDataHandle> nrrds;
  if (nrrd_handle1->nrrd_) nrrds.push_back(nrrd_handle1);
  for (size_t j=0;j<nrrd_handle2.size();j++) if (nrrd_handle2[j]->nrrd_) nrrds.push_back(nrrd_handle2[j]);

  if (inputs_changed_ || bins_.changed() || mins_.changed() ||
      maxs_.changed() || type_.changed() || !oport_cached("OutputNrrd"))
  {
    reset_vars();
  
    int max_dim = static_cast<int>(nrrds[0]->nrrd_->dim);
    for (size_t j=1;j<nrrds.size();j++) 
    {
      if (static_cast<int>(nrrds[j]->nrrd_->dim) > max_dim) 
      {
        max_dim = nrrds[j]->nrrd_->dim;
      }
    }
    
    // Determine the number of bins given
    std::string bins = bins_.get();
    int binsLen = 0;
    char ch;
    int i=0, start=0;
    bool inword = false;
    while (i < (int)bins.length()) 
    {
      ch = bins[i];
      if(isspace(ch)) 
      {
        if (inword) 
        {
          binsLen++;
          inword = false;
        }
      } 
      else if (i == (int)bins.length()-1) 
      {
        binsLen++;
        inword = false;
      } 
      else 
      {
        if(!inword) inword = true;
      }
      i++;
    }
    
    if ((int)nrrds.size() != binsLen) 
    {
      error("Number of input nrrds is not equal to number of bin specifications.");
      return;
    }
    
    // get bins
    size_t *bin = new size_t[nrrds.size()];
    i=0, start=0;
    int which = 0, end=0, counter=0;
    inword = false;
    while (i < (int)bins.length())
    {
      ch = bins[i];
      if(isspace(ch))
      {
        if (inword)
        {
          end = i;
          from_string(bins.substr(start,end-start),bin[counter]);
          which++;
          counter++;
          inword = false;
        }
      }
      else if (i == (int)bins.length()-1)
      {
        if (!inword)
        {
          start = i;
        }
        end = i+1;
        from_string(bins.substr(start,end-start),bin[counter]);
        which++;
        counter++;
        inword = false;
      }
      else
      {
        if(!inword)
        {
          start = i;
          inword = true;
        }
      }
      i++;
    }

    Nrrd **nrrds_array = new Nrrd *[nrrds.size()];
    NrrdRange **range = new NrrdRange *[nrrds.size()];
    for (unsigned int d = 0; d< nrrds.size(); d++)
    {
      nrrds_array[d] = nrrds[d]->nrrd_;
      range[d] = nrrdRangeNew(AIR_NAN, AIR_NAN);
    }
    
    // Determine the number of mins given
    std::string mins = mins_.get();
    int minsLen = 0;
    i = 0;
    inword = false;
    while (i < (int)mins.length())
    {
      ch = mins[i];
      if(isspace(ch))
      {
        if (inword)
        {
          minsLen++;
          inword = false;
        }
      }
      else if (i == (int)mins.length()-1)
      {
        minsLen++;
        inword = false;
      }
      else
      {
        if (!inword) 
          inword = true;
      }
      i++;
    }
    
    if ((int)nrrds.size() != minsLen)
    {
      error("Number of input nrrds is not equal to number of mins specifications.");
      return;
    }
    
    // get mins
    double *min = new double[nrrds.size()];
    i=0, start=0;
    which = 0, end=0, counter=0;
    inword = false;
    while (i < (int)mins.length())
    {
      ch = mins[i];
      if(isspace(ch))
      { 
        if (inword)
        {
          end = i;
          if (mins.substr(start,end-start) == "nan")
            min[counter] = AIR_NAN;
          else
            from_string(mins.substr(start,end-start),min[counter]);
          which++;
          counter++;
          inword = false;
        }
      }
      else if (i == (int)mins.length()-1)
      {
        if (!inword)
        {
          start = i;
        }
        end = i+1;
        if (mins.substr(start,end-start) == "nan")
          min[counter] = AIR_NAN;
        else
          from_string(mins.substr(start,end-start),min[counter]);
        which++;
        counter++;
        inword = false;
      }
      else
      {
        if (!inword)
        {
          start = i;
          inword = true;
        }
      }
      i++;
    }
    
    for (int d=0; d<(int)nrrds.size(); d++) {
      range[d]->min = min[d];
    }
    
    // Determaxe the number of maxs given
    std::string maxs = maxs_.get();
    int maxsLen = 0;
    inword = false;
    i = 0;
    while (i < (int)maxs.length())
    {
      ch = maxs[i];
      if(isspace(ch))
      {
        if (inword)
        {
          maxsLen++;
          inword = false;
        }
      }
      else if (i == (int)maxs.length()-1)
      {
        maxsLen++;
        inword = false;
      }
      else
      {
        if(!inword) 
          inword = true;
      }
      i++;
    }
    
    if ((int)nrrds.size() != maxsLen)
    {
      error("Number of input nrrds is not equal to number of maxs specifications.");
      return;
    }
    
    // get maxs
    double *max = new double[nrrds.size()];
    i=0, start=0;
    which = 0, end=0, counter=0;
    inword = false;
    while (i < (int)maxs.length())
    {
      ch = maxs[i];
      if(isspace(ch))
      {
        if (inword)
        {
          end = i;
          if (maxs.substr(start,end-start) == "nan")
            max[counter] = AIR_NAN;
          else
            from_string(maxs.substr(start,end-start),max[counter]);
          which++;
          counter++;
          inword = false;
        }
      }
      else if (i == (int)maxs.length()-1)
      {
        if (!inword)
              {
          start = i;
        }
        end = i+1;
        if (maxs.substr(start,end-start) == "nan")
          max[counter] = AIR_NAN;
        else
          from_string(maxs.substr(start,end-start),max[counter]);
        which++;
        counter++;
        inword = false;
      }
      else
      {
        if(!inword)
        {
          start = i;
          inword = true;
        }
      }
      i++;
    }
    
    for (int d=0; d<(int)nrrds.size(); d++) 
    {
      range[d]->max = max[d];
    }
    
    int clamp[NRRD_DIM_MAX];
    for (int d=0; d<(int)nrrds.size(); d++) 
    {
      clamp[d] = 0;
    }

    // nrrdHistoJoint crashes if min == max
    for (int d = 0; d < (int)nrrds.size(); d++) 
    {
      if (range[d]->min == range[d]->max) 
      {
        warning("range has 0 width, not computing.");
        return;
      }
    }

    Nrrd *nout = nrrdNew();

    if (nrrdHistoJoint(nout, nrrds_array, range,
		       (unsigned int)nrrds.size(), weight, bin,
                       string_to_nrrd_type(type_.get()), 
		       clamp))
    {
      char *err = biffGetDone(NRRD);
      error(std::string("Error performing Unu Jhisto: ") +  err);
      free(err);
      return;
    }

    NrrdDataHandle nrrdH = new NrrdData(nout);

    if (nrrds.size())
    {
      if (airIsNaN(range[0]->min) || airIsNaN(range[0]->max))
      {
        NrrdRange *minmax = nrrdRangeNewSet(nrrds_array[0],
                                                  nrrdBlind8BitRangeFalse);
        if (airIsNaN(range[0]->min)) range[0]->min = minmax->min;
        if (airIsNaN(range[0]->max)) range[0]->max = minmax->max;
        nrrdRangeNix(minmax);
      }
      nrrdKeyValueAdd(nrrdH->nrrd_, "jhisto_nrrd0_min", 
		      to_string(range[0]->min).c_str());
      nrrdKeyValueAdd(nrrdH->nrrd_, "jhisto_nrrd0_max", 
		      to_string(range[0]->max).c_str());
    }
    
    for (int d=0; d < (int)nrrds.size(); ++d)
    {
      nrrdRangeNix(range[d]);
    }
    delete[] range;
    delete[] nrrds_array;

    delete bin;
    delete min;
    delete[] max;

    send_output_handle("OutputNrrd", nrrdH, true);
  }

}


} // End namespace Teem


