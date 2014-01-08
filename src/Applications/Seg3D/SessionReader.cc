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
//    File   : SessionReader.cc
//    Author : McKay Davis
//    Date   : Tue Oct 17 15:55:03 2006

#include <Applications/Seg3D/SessionReader.h>
#include <Applications/Seg3D/Seg3DVersion.h>
#include <Core/XMLUtil/XMLUtil.h>
#include <Core/Util/StringUtil.h>
#include <Core/Util/Environment.h>
#include <Core/Util/Assert.h>
#include <Core/Events/EventManager.h>
#include <Core/Skinner/XMLIO.h>
#include <Core/Skinner/Variables.h>
#include <libxml/xmlreader.h>
#include <libxml/catalog.h>
#include <libxml/xinclude.h>
#include <libxml/xpathInternals.h>
#include <iostream>
#include <stdio.h>
#include <string.h>

#include <wx/wfstream.h>
#include <wx/zstream.h>
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

namespace SCIRun {

// Any value over 1<<8 will be invalid.  Just pick one that's easy to debug.
#define INVALID_LABEL_MARKER 31773

SessionReader::SessionReader(Painter *painter) :
  painter_(painter),
  dir_(),
  version_major_(0),
  version_minor_(0),
  version_release_(0)
{
}

SessionReader::~SessionReader()
{
}


bool
SessionReader::load_session(string filename)
{
  bool result = false;

  try {
    result = load_session_tgz(filename);
  }
  catch (...)
  {
    result = load_session_xml(filename);
  }

  return result;
}


bool
SessionReader::load_session_xml(string filename)
{
  filename = substituteTilde(filename);

  // Split off the directory and filename and save it so that it can
  // be used when grabbing the file listed in the session file.
  pair<string, string> dir_file = split_filename(filename);
  dir_ = dir_file.first;

  // Read the xml file into a buffer.
  wxFFileInputStream in(std2wx(filename));

  string xmlbuffer;
  char tmpbuffer[2];

  tmpbuffer[1] = '\0';

  while( !in.Read(&tmpbuffer, 1).Eof() )
  {
    xmlbuffer.append( tmpbuffer );
  }

  // Parse the xml buffer.
  const bool result = parse_session_xml( xmlbuffer.c_str(), xmlbuffer.size() );

  return result;
}


bool
SessionReader::load_session_tgz(string filename)
{
  // Creating this suppresses error message popups, such as that we
  // are trying to open and uncompressed file.
  // TODO: This turns off all errors.  Really only the compressed file
  // check should be turned off.
  wxLogNull suppress_error_popups;

  filename = substituteTilde(filename);

  char * xmlbuffer;
  size_t xmlsize;

  // Find the session xml file in the archive. The other files will be
  // found when the session file is parsed.

  // NOTE: The session file must be read first so that the versioning
  // information can be found. The versioning information is needed to
  // process the nrrd data properly.
  if( (xmlsize = load_tgz_xml(filename, &xmlbuffer)) )
  {
    // Get all of the nrrds in the archive.
    if( !load_tgz_nrrd(filename, tgz_nrrds_) )
    {
      throw "Error could not find any nrrd files in the archive.";
    }

    const bool result = parse_session_xml(xmlbuffer, xmlsize);

    delete [] xmlbuffer;

    return result;
  }
  else
  {
    throw "Error could not find the xml session file in the archive.";

    return false;
  }
}


size_t
SessionReader::load_tgz_xml(string tgzname, char **buffer)
{
  // NOTE: Because the tar ball is compressed it is a non-seekable
  // stream as such each entry in the archive must be read in.
  wxFFileInputStream in(std2wx(tgzname));
  wxZlibInputStream gzin(in);
  wxTarInputStream tarin(gzin);
  wxTarEntry *tar_entry;

  // Loop through all of the archive entries to find the session.xml
  // file.
  while( (tar_entry = tarin.GetNextEntry()) != NULL)
  {
    string name = wx2std(tar_entry->GetName());

    // Get the entry data no matter what.
    size_t size = tar_entry->GetSize();

    // and read it into a buffer.
    *buffer = new char[size];
    tarin.Read(*buffer, size);

    // If the session file was found return.
    if (name == "session.xml")
    {
      return size;
    }

    // Otherwise delete what was read and try the next entry.
    else
    {
      delete [] *buffer;
    }
  }

  return 0;
}


size_t
SessionReader::load_tgz_nrrd(string tgzname, tgz_map_type &tgz_nrrds)
{
  // NOTE: Because the tar ball is compressed it is a non-seekable
  // stream as such each entry in the archive must be read in.
  wxFFileInputStream in(std2wx(tgzname));
  wxZlibInputStream gzin(in);
  wxTarInputStream tarin(gzin);
  wxTarEntry *tar_entry;

  size_t nentries = 0;

  // Loop through all of the archive entries to find the nrrd files.
  while( (tar_entry = tarin.GetNextEntry()) != NULL)
  {
    string name = wx2std(tar_entry->GetName());

    if (name.find(".nrrd") != string::npos)
    {
      size_t size = tar_entry->GetSize();

      const size_t maxhbufsize = 1<<18; // 256k is enough for anyone?
      char *hbuffer = new char[maxhbufsize];

      // First read in the nrrd header. Read character by character so
      // that it s possible to find the linebreak which signifies the
      // end of the header.
      size_t hsize = 0;
      while (hsize < size && hsize < maxhbufsize-1 && tarin.LastRead())
      {
	hbuffer[hsize++] = tarin.GetC();
	    
	// Check for both Unix and DOS linebreak formats in the nrrd header.
	if (hsize >= 2 &&
	    hbuffer[hsize-1] == '\n' && hbuffer[hsize-2] == '\n' ||
	    hsize >= 4 &&
	    hbuffer[hsize-2] == '\r' && hbuffer[hsize-1] == '\n' &&
	    hbuffer[hsize-4] == '\r' && hbuffer[hsize-3] == '\n')
	{
	  hbuffer[hsize] = '\0';
	  break;
	}
      }

      if (hsize == size || hsize == maxhbufsize-1 || !tarin.LastRead())
      {
	// Throwing errors when reading xml files screws up the xml
	// reader so throw an error but also do a cerr first.
	cerr << "Can not find the nrrd header for " << name << "." << endl;
	throw "Can not find the nrrd header for " + name + ".";
      }
	  
      NrrdDataHandle nrrdh = new NrrdData();
      NrrdIoState *nio = nrrdIoStateNew();
      nrrdIoStateSet(nio, nrrdIoStateSkipData, AIR_TRUE);
      if (nrrdStringRead(nrrdh->nrrd_, hbuffer, nio))
      {
	char *err = biffGetDone(NRRD);
	// Throwing errors when reading xml files screws up the xml
	// reader so throw an error but also do a cerr first.
	cerr << err << endl;
	throw err;

      }

      delete [] hbuffer;

      // Verify that the format was RAW.
      if( nio->encoding != nrrdEncodingArray[1]) // raw
      {
	// Throwing errors when reading xml files screws up the xml
	// reader so throw an error but also do a cerr first.
	cerr << "Nrrd data is not raw " << name << "." << endl;
	throw "Nrrd data is not raw for " + name + ".";
      }

      // Verify uchar and float formats.
      if( nrrdh->nrrd_->type != nrrdTypeUChar &&
	  nrrdh->nrrd_->type != nrrdTypeFloat)
      {
	// Throwing errors when reading xml files screws up the xml
	// reader so throw an error but also do a cerr first.
	cerr << "Nrrd data is not uchar or float for " << name << "." << endl;
	throw "Nrrd data is not uchar or float for " + name + ".";
      }

      const size_t dsize = size - hsize;
      ASSERT(dsize == VolumeOps::nrrd_data_size(nrrdh));
      char *dbuffer = new char[dsize];
      nrrdh->nrrd_->data = dbuffer;
	  
      // Read in 1G blocks to work around 4G filesystem read size limits.
      size_t totalrsize = dsize;
      char *rbuffer = dbuffer;
      while (totalrsize)
      {
	size_t rsize = 1<<30;
	if (rsize > totalrsize) rsize = totalrsize;
	tarin.Read(rbuffer, rsize);
	totalrsize -= rsize;
	rbuffer += rsize;
      }

      // Get endianness correct.
      if (nrrdh->nrrd_->type == nrrdTypeFloat && nio->endian != AIR_ENDIAN)
      {
	nrrdSwapEndian(nrrdh->nrrd_);
      }

      nrrdIoStateNix(nio);

      tgz_nrrds_[name] = nrrdh;

      ++nentries;
    }

    // Was not the correct entry but read and discard it. This is done
    // because the archieve is not seekable.
    else
    {
      size_t size = tar_entry->GetSize();
      
      char * buffer = new char[size];
      tarin.Read(buffer, size);

      delete [] buffer;
    }
  }

  return nentries;
}


bool
SessionReader::parse_session_xml(const char *xmlbuffer, size_t xmlsize)
{
  /*
   * This initialize the library and check potential ABI mismatches
   * between the version it was compiled for and the actual shared
   * library used.
   */
  LIBXML_TEST_VERSION;

  // The parser context
  xmlParserCtxtPtr ctxt = xmlNewParserCtxt();
  if (!ctxt) {
    cerr << "SessionReader failed xmlNewParserCtx()" << endl;
    return false;
  }

  // Parse the file, activating the DTD validation option.
  xmlDocPtr doc = 
    xmlCtxtReadMemory(ctxt, xmlbuffer, xmlsize, "session.xml",
                      0, XML_PARSE_PEDANTIC);
  
  if (!doc) {
    cerr << "Skinner::XMLIO::load failed to parse session.xml" << endl;
    return false;
  } 

  // Parse the doc at network node.
  for (xmlNode *cnode=doc->children; cnode!=0; cnode=cnode->next) {
    if (XMLUtil::node_is_element(cnode, "Seg3D")) {

      // Get the Seg3D version
      eval_version_node(cnode);

      // Evaluate the Volume and Apperance nodes.
      eval_seg3d_node(cnode);
    }

    else if (!XMLUtil::node_is_comment(cnode))
      throw "Unknown node type.";
  }               
  
  xmlFreeDoc(doc);
  xmlFreeParserCtxt(ctxt);  
#ifndef _WIN32
  xmlCleanupParser();
#endif

  return true;
}


void
SessionReader::eval_seg3d_node(const xmlNodePtr node)
{
  volumes_.clear();
  painter_->clear_all_volumes();

  for (xmlNode *cnode=node->children; cnode!=0; cnode=cnode->next) {
    if (XMLUtil::node_is_element(cnode, "volume")) {
      NrrdVolumes empty;
      eval_volume_node(cnode, empty);
    }
  }

  if (volumes_.size())
  {
    painter_->volumes_ = volumes_;
    painter_->current_volume_ = volumes_.back();
  }

  for (xmlNode *cnode=node->children; cnode!=0; cnode=cnode->next) {
    if (XMLUtil::node_is_element(cnode, "appearance")) {
      eval_appearance_node(cnode);
    }
  }

  painter_->extract_all_window_slices();
  painter_->rebuild_layer_buttons();

  if (first_volume_.get_rep())
  {
    painter_->current_vrender_target_ = first_volume_;
    painter_->current_vrender_target_deferred_ = true;
    event_handle_t unused;
    painter_->ShowVolumeRendering2(unused);
      
    first_volume_->get_geom_group();
    painter_->Autoview(unused);
  }

  painter_->redraw_all();
}


void
SessionReader::eval_volume_node(const xmlNodePtr node,
				NrrdVolumes parents)
{
  Skinner::Variables *vars = new Skinner::Variables("");
  for (xmlNode *cnode=node->children; cnode!=0; cnode=cnode->next) {
    if (XMLUtil::node_is_element(cnode, "var")) {
      Skinner::XMLIO::eval_var_node(cnode, vars);
    } 
  }

  NrrdVolumeHandle volume = 0;

  unsigned int label = 0;
  if (vars->exists("label")) {
    label = vars->get_int("label");
  }

  int label_index = 0;
  if (label > 0xffffff) { label_index = 3; label = label >> 24; }
  else if (label > 0xffff) { label_index = 2; label = label >> 16; }
  else if (label > 0xff) { label_index = 1; label = label >> 8; }

  if (!parents.size())
  {
    string filename = vars->get_string("filename");

    // In this case the xml file is part of a tgz file so pull the
    // nrrd data from the tar archive.
    if (tgz_nrrds_.size())
    {
      tgz_map_type::iterator loc = tgz_nrrds_.find(filename);

      if (loc != tgz_nrrds_.end())
      {
	NrrdDataHandle nrrdh = (*loc).second;

	if( nrrdh.get_rep() )
	{
	  // TODO: The nrrd spacings fields is being used poorly
	  // poorly. As such, move from nice orientation here back to
	  // min/max and spacings fields.
	  if( version_major_ == 1 && version_minor_ == 10 )
	  {
	    // This appears to be need because someone made the dumb
	    // decision to save nrrd using ITK which assumes a
	    // particular format of the nrrd. So when read back in it
	    // is necessary to "patch" them.
	    
	    // However, with 1.11 the nrrd writer is used so this is
	    // not needed.
	    for (size_t i = 0; i < nrrdh->nrrd_->dim; i++)
	    {
	      nrrdh->nrrd_->axis[i].spacing =
		nrrdh->nrrd_->axis[i].spaceDirection[i];
	      nrrdh->nrrd_->axis[i].min = nrrdh->nrrd_->spaceOrigin[i];
	      nrrdh->nrrd_->axis[i].max = nrrdh->nrrd_->axis[i].min +
		(nrrdh->nrrd_->axis[i].size-1) * nrrdh->nrrd_->axis[i].spacing;
	    }
	    nrrdh->nrrd_->spaceDim = 0;
	  }
 
	  if (label)
	  {
	    bool found = false;
	    
	    for (size_t i = 0; i < volumes_.size(); i++)
	    {
	      if ( nrrdh.get_rep() == volumes_[i]->nrrd_handle_.get_rep() )
	      {
		volume = volumes_[i]->create_child_label_volume(label);
		found = true;
		break;
	      }
	    }
	    if (!found)
	    {
	      volume = new NrrdVolume(painter_, "", nrrdh, label);
	    }
	  }
	  else
	  {
	    volume = new NrrdVolume(painter_, "", nrrdh, 0);
	    
	    if (!first_volume_.get_rep())
	    {
	      // TODO: Emulate the bad behavior of Seg3D and volume
	      // selection here.  There is not currently a way to
	      // select which volume is volume-rendered so we just
	      // pick an arbitrary one (the first).  Change this to a
	      // variable when we can select which are visible.
	      first_volume_ = volume;
	    }
	  }
	}
      }
      else
      {
	// Throwing errors when reading xml files screws up the xml
	// reader so throw an error but also do a cerr first.
	cerr << "Error could find " + filename + " in the archive." << endl;
	
	throw "Error could find " + filename + " in the archive.";
      }
    }

    // Parsing the xml session file as a normal file so data should be loaded.
    else
    {
      pair<string, string> dir_file = split_filename(filename);

      if (dir_file.first.empty())
      {
	dir_file.first = dir_;
      }

      if (!label)
      {
        volume = painter_->load_volume<float>(dir_file.first+filename);
        if (!volume.get_rep())
        {
	  // Throwing errors when reading xml files screws up the xml
	  // reader so throw an error but also do a cerr first.
          cerr << "Error reading auxiliary data volume "
	       << filename << "." << endl;
          throw "Error reading auxiliary data volume " + filename + ".";
        }
        if (!first_volume_.get_rep())
        {
          // TODO:  Emulate the bad behavior of Seg3D and volume selection
          // here.  There is not currently a way to select which volume
          // is volume-rendered so we just pick an arbitrary one (the first).
          // Change this to a variable when we can select which are visible.
          first_volume_ = volume;
        }
      }
      else if ( version_major_ == 1 && version_minor_ <= 9)
      {
        NrrdVolumeHandle tmp =
          painter_->load_volume<unsigned int>(dir_file.first+filename);

        if (!tmp.get_rep())
        {
	  // Throwing errors when reading xml files screws up the xml
	  // reader so throw an error but also do a cerr first.
          cerr << "Error reading auxiliary label volume "
	       << filename << "." << endl;
          throw "Error reading auxiliary label volume " + filename + ".";
        }
        vector<NrrdDataHandle> nrrds;
        VolumeOps::dice_32_into_8(tmp->nrrd_handle_, nrrds);
        for (size_t i = 0; i < nrrds.size(); i++)
        {
          NrrdVolumeHandle vol =
            new NrrdVolume(painter_, "", nrrds[i], INVALID_LABEL_MARKER);
          parents.push_back(vol);
        }
        volume = parents[label_index];
      }
      else //if ( version_major_ == 1 && version_minor_ >= 10)
      {
        bool found = false;
        for (size_t i = 0; i < volumes_.size(); i++)
        {
          if (filename == volumes_[i]->filename_)
          {
            volume = volumes_[i]->create_child_label_volume(label);
            found = true;
            break;
          }
        }

        if (!found)
        {
          volume =
	    painter_->load_volume<unsigned char>(dir_file.first+filename);

          if (!volume.get_rep())
          {
	    // Throwing errors when reading xml files screws up the xml
	    // reader so throw an error but also do a cerr first.
	    cerr << "Error reading auxiliary label volume "
		 << filename << "." << endl;
            throw "Error reading auxiliary label volume " + filename + ".";
          }
        }
      }
    }

    if (!volume.get_rep())
    {
      // Throwing errors when reading xml files screws up the xml
      // reader so throw an error but also do a cerr first.
      cerr << "Bad session volume pointer." << endl;
      throw "Bad session volume pointer.";
    }

    volume->filename_ = filename;
  }
  else
  {
    if (parents[label_index]->label_ == INVALID_LABEL_MARKER)
    {
      volume = parents[label_index];
    }
    else
    {
      volume = parents[label_index]->create_child_label_volume(label);
    }
  }

  volume->set_label(label);

  if (label)
  {
    if (vars->exists("label_color"))
    {
      Skinner::Color sc = vars->get_color("label_color");
      Color c(sc.r, sc.g, sc.b);
      volume->set_label_color(c);
    }
    else
    {
      unsigned int colormap_index = 0;
      if (vars->exists("colormap_index"))
      {
        colormap_index = vars->get_int("colormap_index");
      }
      volume->set_label_color_legacy(colormap_index);
    }
  }

  if (vars->exists("name")) {
    volume->name_ = vars->get_string("name");
  }
  
  if (vars->exists("visible")) {
    volume->tmp_visible_ = vars->get_bool("visible");
  }

  if (vars->exists("opacity")) {
    volume->set_opacity( vars->get_double("opacity") );
  }

  volumes_.push_back(volume);

  for (xmlNode *cnode=node->children; cnode!=0; cnode=cnode->next) {
    if (XMLUtil::node_is_element(cnode, "volume")) {
      eval_volume_node(cnode, parents);
    } 
  }

  delete vars;
}


void
SessionReader::eval_appearance_node(const xmlNodePtr node)
{
  Skinner::Variables *vars = new Skinner::Variables("");
  for (xmlNode *cnode = node->children; cnode != 0; cnode = cnode->next) {
    if (XMLUtil::node_is_element(cnode, "var")) {
      Skinner::XMLIO::eval_var_node(cnode, vars);
    }
  }
  
  painter_->get_session_appearance_vars(vars);
  delete vars;
}


void
SessionReader::eval_version_node(const xmlNodePtr cnode)
{
  const char *vstr = XMLUtil::node_att_as_const_char_str(cnode, "version");

  string version = string(vstr);
      
  // Previous developers did not know how to do versioning
  // properly so do some clean up.
  if( version == string("1.0") ) {
    version_major_ = 1;
    version_minor_ = 9;
    version_release_ = 0;
  }

  else if( version == string("2.0") ) {
    version_major_ = 1;
    version_minor_ = 10;
    version_release_ = 0;
  }

  else {
    version_major_ = GET_MAJOR_VERSION( version );
    version_minor_ = GET_MINOR_VERSION( version );
    version_release_ = GET_RELEASE_VERSION( version );
  }

  if(check_version( version_major_, version_minor_, version_release_) > 0)
  {
    throw "Unable to read session file with new version.";
  }
} 



} // end namespace SCIRun
