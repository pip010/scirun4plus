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
//    File   : SessionWriter.cc
//    Author : McKay Davis
//    Date   : Tue Oct 17 21:27:22 2006

#include <Applications/Seg3D/SessionWriter.h>
#include <Applications/Seg3D/Seg3DVersion.h>
#include <Core/XMLUtil/XMLUtil.h>
#include <Core/Util/StringUtil.h>
#include <Core/Util/Environment.h>
#include <Core/Util/Assert.h>
#include <Core/Util/FileUtils.h>
#include <Core/Skinner/XMLIO.h>
#include <Core/Skinner/Variables.h>
#include <libxml/xmlreader.h>
#include <libxml/catalog.h>
#include <libxml/xinclude.h>
#include <libxml/xpathInternals.h>
#include <iostream>

#include <wx/wfstream.h>
#include <wx/tarstrm.h>
#include <wx/zstream.h>
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

namespace SCIRun {

struct ndh_less_than
{
  bool operator()(const NrrdDataHandle &s1, const NrrdDataHandle &s2) const
  {
    return s1.get_rep() < s2.get_rep();
  }
};


static string
color_to_string(const Color &c)
{
  char temp[24];
  sprintf(temp, "#%02X%02X%02X%02X",
          (unsigned char)(c.r() * 255.0 + 0.5),
          (unsigned char)(c.g() * 255.0 + 0.5),
          (unsigned char)(c.b() * 255.0 + 0.5),
          255);
  return string(temp);
}


bool
SessionWriter::write_session(const string &filename, Painter *painter)
{
  /*
   * Initialize libxml and check potential ABI mismatches
   * between the version it was compiled for and the actual shared
   * library used.
   */
  
  LIBXML_TEST_VERSION;
  
  // TODO:  We could pack the label volume bitplanes here for smaller
  // output file size.
  // painter->pack_labels();

  // Open up the session archive.
  wxFFileOutputStream out(std2wx(filename));
  wxZlibOutputStream gzout(out, wxZ_BEST_SPEED, wxZLIB_GZIP);
  wxTarOutputStream tar(gzout);

  // Create the session.xml document.
  xmlDocPtr doc = xmlNewDoc(to_xml_ch_ptr("1.0"));
  xmlNodePtr root = xmlNewNode(0,to_xml_ch_ptr("Seg3D"));
  xmlDocSetRootElement(doc, root);

  xmlNewProp(root, to_xml_ch_ptr("version"), to_xml_ch_ptr(SEG3D_VERS_STRING));

  std::map<NrrdDataHandle, string, ndh_less_than> entries;

  int dcounter = 1;
  int lcounter = 1;
  for (size_t i = 0; i < painter->volumes_.size(); i++)
  {
    NrrdVolumeHandle volume = painter->volumes_[i];

    string entry;
    std::map<NrrdDataHandle, string, ndh_less_than>::iterator entry_loc =
      entries.find(volume->nrrd_handle_);
    if (entry_loc == entries.end())
    {
      if (volume->label_)
      {
        entry = "label" + to_string(lcounter++) + ".nrrd";
      }
      else
      {
        entry = "data" + to_string(dcounter++) + ".nrrd";
      }
      entries[volume->nrrd_handle_] = entry;
    }
    else
    {
      entry = entries[volume->nrrd_handle_];
    }

    xmlNodePtr cnode = xmlNewChild(root, 0, to_xml_ch_ptr("volume"), 0);
    add_var_node(cnode, "name", volume->name_);
    add_var_node(cnode, "filename", entry);
    add_var_node(cnode, "label", to_string(volume->label_));
    add_var_node(cnode, "label_color",
                 color_to_string(volume->get_label_color()));
    add_var_node(cnode, "visible", to_string(volume->visible()));
    add_var_node(cnode, "opacity", to_string(volume->get_opacity()));
  }

  add_appearance_nodes(root, painter);
  
  xmlChar *buffer;
  int buffersize;
  xmlDocDumpFormatMemory(doc, &buffer, &buffersize, 1);

  // Add session.xml to the session archive.
  wxTarEntry *tarentry =
    new wxTarEntry(_T("session.xml"), wxDateTime::Now(), buffersize);
  tar.PutNextEntry(tarentry);
  tar.Write(buffer, buffersize);

  xmlFree(buffer);
  xmlFreeDoc(doc);

  // Add the volumes to the session archive.

  std::map<NrrdDataHandle, string, ndh_less_than>::iterator eitr = entries.begin();
  while (eitr != entries.end())
  {
    NrrdDataHandle nrrd_handle = (*eitr).first;
    const string entry = (*eitr).second;
    const string nhdr = NrrdVolume::create_nrrd_header(nrrd_handle->nrrd_);
    const size_t nsize = VolumeOps::nrrd_data_size(nrrd_handle);
    
    wxTarEntry *tarentry =
      new wxTarEntry(std2wx(entry), wxDateTime::Now(), nhdr.size() + nsize);
    tar.PutNextEntry(tarentry);

    tar.Write(nhdr.c_str(), nhdr.size());

    size_t totalwsize = nsize;
    char *wbuffer = (char *)(nrrd_handle->nrrd_->data);
    while (totalwsize)
    {
      size_t wsize = 1<<30;
      if (wsize > totalwsize) wsize = totalwsize;
      tar.Write(wbuffer, wsize);
      totalwsize -= wsize;
      wbuffer += wsize;
    }
    
    ++eitr;
  }

  tar.Close();

  return true;
}


static string
remove_extension(const string &filename)
{
  string result;
  std::vector<string> fileext = split_string(filename, '.');
  for (int i = 0; i < (int)fileext.size()-1; i++)
  {
    if (i > 0) result = result + ".";
    result = result + fileext[i];
  }
  return result;
}


void
SessionWriter::add_appearance_nodes(xmlNodePtr node, Painter *painter)
{
  painter->set_session_appearance_vars();
  Skinner::Variables *vars = painter->get_vars();
  xmlNodePtr cnode = xmlNewChild(node, 0, to_xml_ch_ptr("appearance"), 0);
  
  add_var_node(cnode, "Painter::volume_visible", vars);
  add_var_node(cnode, "Painter::active_volume_index", vars);
  // Probe is currently broken, so don't use it.
  //add_var_node(cnode, "Painter::probe::x", vars);
  //add_var_node(cnode, "Painter::probe::y", vars);
  //add_var_node(cnode, "Painter::probe::z", vars);
}


bool
SessionWriter::write_volumes(NrrdVolumes &volumes, const string &dir)
{
  NrrdVolumes::iterator iter = volumes.begin();  
  NrrdVolumes::iterator end = volumes.end();
  
  std::map<NrrdDataHandle, string, ndh_less_than> written;

  for (; iter != end; ++iter) {
    static int volume_id = 0;
    NrrdVolumeHandle volume = *iter;
    
    // Painter::volumes_ no longer contains a tree.  Write each nrrd
    // out only once.
    std::map<NrrdDataHandle, string, ndh_less_than>::iterator written_loc =
      written.find(volume->nrrd_handle_);
    if (written_loc == written.end())
    {
      string basename = remove_extension(volume->filename_);
      if (basename.empty()) {
        basename = remove_extension(volume->name_);
      }

      if (basename.empty()) {
        basename = "Volume";
      }
    
      string filename = basename + ".nrrd";
      while (validFile(dir + "/" + filename))
      {
        filename = basename + to_string(volume_id++) + ".nrrd";
      }

      volume->filename_ = filename;
      if (!volume->write(dir + "/" + volume->filename_)) return false;

      written[volume->nrrd_handle_] = filename;
    }
    else
    {
      volume->filename_ = written_loc->second;
    }
  }
  return true;
}
    

void
SessionWriter::add_volume_nodes(xmlNodePtr node, NrrdVolumes &volumes)
{
  NrrdVolumes::iterator iter = volumes.begin();  
  NrrdVolumes::iterator end = volumes.end();

  for (; iter != end; ++iter) {
    NrrdVolumeHandle volume = *iter;
    xmlNodePtr cnode = xmlNewChild(node, 0, to_xml_ch_ptr("volume"),0);
    add_var_node(cnode, "name", volume->name_);
    if (volume->filename_ != "")
    {
      add_var_node(cnode, "filename", volume->filename_);
    }
    add_var_node(cnode, "label", to_string(volume->label_));
    add_var_node(cnode, "label_color",
                 color_to_string(volume->get_label_color()));
    add_var_node(cnode, "visible", to_string(volume->visible()));
    add_var_node(cnode, "opacity", to_string(volume->get_opacity()));
  }
}


void
SessionWriter::add_var_node(xmlNodePtr node, 
                            const string &name,
                            const string &value)
{
  xmlNodePtr cnode = xmlNewChild(node, 0, to_xml_ch_ptr("var"), 
                                 to_xml_ch_ptr(value.c_str()));
  xmlNewProp(cnode, to_xml_ch_ptr("name"), to_xml_ch_ptr(name.c_str()));
}


void
SessionWriter::add_var_node(xmlNodePtr node, 
                            const string &name,
                            Skinner::Variables *vars)
{
  if (vars->exists(name))
  {
    add_var_node(node, name, vars->get_string(name));
  }
}


} // end namespace SCIRun
