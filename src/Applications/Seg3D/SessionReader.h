//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : SessionReader.h
//    Author : McKay Davis
//    Date   : Tue Oct 17 15:48:19 2006

#ifndef SEG3D_SessionReader_h
#define SEG3D_SessionReader_h

#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/NrrdVolume.h>
#include <libxml/xmlreader.h>

#include <wx/tarstrm.h>

#include <string>

using std::string;

namespace SCIRun {
class Painter;
class SessionReader {
  typedef std::map<string, NrrdDataHandle, less<string> > tgz_map_type;

public:
  SessionReader         (Painter *painter);
  ~SessionReader();
  bool                  load_session(string filename);

private:
  
  bool                  load_session_tgz(string filename);
  bool                  load_session_xml(string filename);
  bool                  parse_session_xml(const char *xmlbuf, size_t bufsize);


  size_t                load_tgz_xml(string tgzname, char **buffer);
  size_t                load_tgz_nrrd(string tgzname, tgz_map_type &tgz_nrrds);

  void                  eval_volume_node (const xmlNodePtr, NrrdVolumes pars);
  void                  eval_appearance_node (const xmlNodePtr node );
  void                  eval_version_node(const xmlNodePtr node);
  void                  eval_seg3d_node (const xmlNodePtr node);
  Painter *             painter_;
  string                dir_;
  NrrdVolumeHandle      first_volume_;
  NrrdVolumes           volumes_;
  unsigned int          version_major_;
  unsigned int          version_minor_;
  unsigned int          version_release_;

  tgz_map_type          tgz_nrrds_;
};


} // end namespace SCIRun

#endif


