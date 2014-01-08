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
 *  WriteNrrd.cc: Nrrd writer
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   February 2001
 *
 */

#include <Core/Persistent/Pstreams.h>
#include <Core/ImportExport/Nrrd/NrrdIEPlugin.h>


#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/NrrdPort.h>

#include <Dataflow/GuiInterface/GuiVar.h>

#include <sstream>
#include <fstream>

using std::ifstream;
using std::ostringstream;

namespace SCITeem {

using namespace SCIRun;

class WriteNrrd : public Module {
    NrrdIPort*  inport_;
    GuiFilename filename_;
    GuiString   filetype_;
    GuiString   exporttype_;
    GuiString   types_;
public:
    WriteNrrd(GuiContext *ctx);
    virtual ~WriteNrrd();
    virtual void execute();
private:
  static Piostream* make_stream(const std::string&, const std::string&);
};

} // end namespace SCITeem

using namespace SCITeem;
DECLARE_MAKER(WriteNrrd)

WriteNrrd::WriteNrrd(GuiContext *ctx)
: Module("WriteNrrd", ctx, Filter, "DataIO", "Teem"), 
  filename_(get_ctx()->subVar("filename"), ""),
  filetype_(get_ctx()->subVar("filetype"), "Binary"),
  exporttype_(get_ctx()->subVar("exporttype"), ""),
  types_(get_ctx()->subVar("types"))
{
  NrrdIEPluginManager mgr;
  std::vector<std::string> exporters;
  mgr.get_exporter_list(exporters);
  
  std::string exporttypes = "{";
  exporttypes += "{{Nrrd}                {.nrrd} } ";
  exporttypes += "{{Nrrd Header and Raw} {.nhdr .raw} } ";
  exporttypes += "{{Nrrd File Any} {.*} } ";

  exporttypes += "}";

  types_.set(exporttypes);
}

WriteNrrd::~WriteNrrd()
{
}

Piostream* WriteNrrd::make_stream(const std::string& ft, const std::string& nrrd_data_fn)
{
  if (ft == "Binary") 
    return new BinaryPiostream(nrrd_data_fn, Piostream::Write);
  return new TextPiostream(nrrd_data_fn, Piostream::Write);
}

void
WriteNrrd::execute()
{
  // Read data from the input port
  NrrdDataHandle handle;
  if (!get_input_handle("Input Data", handle)) return;

  // If no name is provided, return
  std::string fn(filename_.get());
  if(fn == "") {
    error("Warning: no filename in WriteNrrd");
    return;
  }

  // get root of filename (no extension)
  std::string::size_type e = fn.find_last_of(".");
  std::string root = fn;
  if (e != std::string::npos) root = fn.substr(0,e);


  // determine which case we are writing out based on the
  // original filename.  
  bool writing_nrrd = false;
  bool writing_nhdr = false;
  bool writing_nd = false;

  if (fn.find(".nrrd",0) != std::string::npos) writing_nrrd = true;
  else if (fn.find(".nhdr",0) != std::string::npos) writing_nhdr = true;
  else if (fn.find(".nd",0) != std::string::npos) writing_nd = true;

  // If the filename doesn't have an extension
  // use the export type to determine what it should be
  if (!writing_nrrd && !writing_nhdr && !writing_nd)
  {
    std::string type(exporttype_.get());
    std::cerr << "exporttype = " << type << std::endl;
    if (type.find(".nrrd",0) != std::string::npos) writing_nrrd = true;
    else writing_nhdr = true;
  }

  // only write out the .nd extension if that was the filename
  // specified or if there are properties.  In the case that it
  // was an .nd filename specified, write out a .nrrd file also.
  if (handle->nproperties() > 0) 
  {
    if (!writing_nd) 
    {
      writing_nd = true;
      if (!writing_nhdr) writing_nrrd = true;
    }
  }

  std::string ft(filetype_.get());

  // writing out NrrdData - use Piostream
  if (writing_nd) 
  {
    std::string nrrd_data_fn = root + ".nd";
    
    // set NrrdData's write_nrrd variable to indicate
    // whether NrrdData's io method should write out a .nrrd or .nhdr
    if (writing_nhdr) handle.get_rep()->write_nrrd_ = false;
    else handle.get_rep()->write_nrrd_ = true;

    PiostreamPtr stream(make_stream(ft, nrrd_data_fn));
    
    if (stream->error()) 
    {
      error("Could not open file for writing" + nrrd_data_fn);
    } 
    else 
    {
      // Write the file
      Pio(*stream, handle); // will also write out a separate nrrd.
    } 
  } 
  else 
  {
    // writing out ordinary .nrrd .nhdr file
    // Restrict TEEM access: it is not thread safe
    NrrdData::lock_teem();
      
    std::string nrrd_fn = root;
    if (writing_nhdr) nrrd_fn += ".nhdr";
    else nrrd_fn += ".nrrd";

    Nrrd *nin = handle->nrrd_;
    
    NrrdIoState *nio = nrrdIoStateNew();
    if (ft == "Binary") {
      // set encoding to be raw
      nio->encoding = nrrdEncodingArray[1];
    }
    else {
      // set encoding to be ascii
      nio->encoding = nrrdEncodingArray[2];
    }
    // set format to be nrrd
    nio->format = nrrdFormatArray[1];
    // set endian to be endian of machine
    nio->endian = airMyEndian;
    
    if (AIR_ENDIAN != nio->endian) 
    {
      nrrdSwapEndian(nin);
    }
    
    if (writing_nhdr) 
    {
      if (nio->format != nrrdFormatNRRD) 
      {
        nio->format = nrrdFormatNRRD;
      }
    }
    
    if (nrrdSave(nrrd_fn.c_str(), nin, nio)) 
    {
      char *err = biffGet(NRRD);      
      error(std::string("Error writing nrrd ") +nrrd_fn+": "+std::string(err));
      
      free(err);
      biffDone(NRRD);
      NrrdData::unlock_teem();
      return;
    }
    NrrdData::unlock_teem();
  }
}


