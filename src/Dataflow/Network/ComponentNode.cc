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


/* ComponentNode.cc
 * 
 * written by 
 *   Chris Moulding
 *   Sept 2000
 *   University of Utah
 */

#include <Dataflow/Network/ComponentNode.h>
#include <Dataflow/GuiInterface/TCLInterface.h>
#include <Core/XMLUtil/XMLUtil.h>
#include <libxml/xmlreader.h>
#include <Core/Util/Environment.h>

#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <sstream>
#include <fstream>

namespace SCIRun {
using std::map;
using std::cout;
using std::endl;
using std::ostringstream;

template <class PInfo>
void 
set_port_info(std::vector<PInfo*> &ports, xmlNodePtr snode)
{
  xmlNodePtr ipnode = snode->children;
  for (; ipnode != 0; ipnode = ipnode->next) {
    if (std::string(to_char_ptr(ipnode->name)) == std::string("port"))
    {
      PInfo *pinfo = new PInfo;
      xmlNodePtr pnode = ipnode->children;
      for (; pnode != 0; pnode = pnode->next) {
	if (std::string(to_char_ptr(pnode->name)) == std::string("name")) {
	  pinfo->name = std::string(to_char_ptr(pnode->children->content));
	}
	if (std::string(to_char_ptr(pnode->name)) == std::string("datatype")) 
	{
	  pinfo->datatype = std::string(to_char_ptr(pnode->children->content));
	}
      }
      ports.push_back(pinfo);
    }
  } 
}


bool
set_port_info(ModuleInfo &mi, const xmlNodePtr cnode)
{
  xmlNodePtr node = cnode->children;
  for (; node != 0; node = node->next) {
    if (std::string(to_char_ptr(node->name)) == std::string("io")) {
      xmlNodePtr ionode = node->children;
      for (; ionode != 0; ionode = ionode->next) {
	if (ionode && (std::string(to_char_ptr(ionode->name)) == 
		       std::string("inputs"))) 
	{
	  xmlAttrPtr lpd_att = get_attribute_by_name(ionode, 
						     "lastportdynamic");
	  if (lpd_att) {
	    mi.last_port_dynamic_ = false;
	    if (std::string(to_char_ptr(lpd_att->children->content)) == 
		std::string("yes")) 
	    {
	      mi.last_port_dynamic_ = true;
	    }
	  } else {
	    std::cerr << "Missing attribute lastportdynamic for module: " 
		      << mi.module_name_ << std::endl;
	    return false;
	  }
	  // set input port info.
	  set_port_info(mi.iports_, ionode);
	}
	else if (ionode && (std::string(to_char_ptr(ionode->name)) == 
			    std::string("outputs"))) 
	{
	  // set input port info.
	  set_port_info(mi.oports_, ionode);
	} 
      }
    } 
  }
  return false;
}


bool
set_description(ModuleInfo &mi, const xmlNodePtr cnode)
{
  xmlNodePtr onode = cnode->children;
  for (; onode != 0; onode = onode->next) {
    if (std::string(to_char_ptr(onode->name)) == std::string("overview")) {
      
      xmlNodePtr node = onode->children;
      for (; node != 0; node = node->next) {
        if (std::string(to_char_ptr(node->name)) == std::string("summary")) 
	{
	  mi.summary_ = get_serialized_children(node);
	} 
	else if (std::string(to_char_ptr(node->name)) == std::string("authors")) 
	{
	  xmlNodePtr anode = node->children;
	  for (; anode != 0; anode = anode->next) {
	    if (std::string(to_char_ptr(anode->name)) == std::string("author")) 
	    {
	      std::string s = std::string(to_char_ptr(anode->children->content));
	      mi.authors_.push_back(s);
	    }
	  }
	} 
      }
    }
  }
  return false;
}


bool
set_gui_info(ModuleInfo &mi, const xmlNodePtr cnode)
{
  mi.has_gui_node_ = false;
  xmlNodePtr onode = cnode->children;
  for (; onode != 0; onode = onode->next) {
    if (std::string(to_char_ptr(onode->name)) == std::string("gui")) {
      mi.has_gui_node_ = true;
      return true;
    }
  }
  return false;
}

void
write_component_file(const ModuleInfo &mi, const char* filename) 
{
  xmlDocPtr doc = 0;        /* document pointer */
  xmlNodePtr root_node = 0; /* node pointers */
  xmlDtdPtr dtd = 0;        /* DTD pointer */
  
  LIBXML_TEST_VERSION;
  
  /* 
   * Creates a new document, a node and set it as a root node
   */
  doc = xmlNewDoc(BAD_CAST "1.0");
  root_node = xmlNewNode(0, BAD_CAST "component");
  xmlDocSetRootElement(doc, root_node);
  
  /*
   * Creates a DTD declaration.
   */
  std::string fname(filename);
  std::string dtdstr;
  if (fname.find("Package") == std::string::npos) {
    dtdstr = std::string("component.dtd");
  } else {
    dtdstr = std::string("../../../../Dataflow/XML/component.dtd");
  }
  dtd = xmlCreateIntSubset(doc, BAD_CAST "component", 0, 
			   BAD_CAST dtdstr.c_str());
 
  xmlChar license[] = \
  "\n"
  "The MIT License\n"
  "Copyright (c) 2009 Scientific Computing and Imaging Institute,\n"
  "University of Utah.\n"
  "\n"
  "Permission is hereby granted, free of charge, to any person obtaining a\n"
  "copy of this software and associated documentation files (the \"Software\"),\n"
  "to deal in the Software without restriction, including without limitation\n"
  "the rights to use, copy, modify, merge, publish, distribute, sublicense,\n"
  "and/or sell copies of the Software, and to permit persons to whom the\n"
  "Software is furnished to do so, subject to the following conditions:\n"
  "\n"
  "The above copyright notice and this permission notice shall be included\n"
  "in all copies or substantial portions of the Software.\n"
  "\n"
  "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS\n"
  "OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
  "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL\n"
  "THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
  "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING\n"
  "FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER\n"
  "DEALINGS IN THE SOFTWARE.\n";
 
  xmlNodePtr licenseComment = xmlNewComment(license);
  xmlAddChild(root_node, licenseComment);
    
  /* 
   * xmlNewChild() creates a new node, which is "attached" as child node
   * of root_node node. 
   */
  xmlNewProp(root_node, BAD_CAST "name", BAD_CAST mi.module_name_.c_str());
  xmlNewProp(root_node, BAD_CAST "category", BAD_CAST mi.category_name_.c_str());
  const char* opt = mi.optional_ ? "true" : "false";
  xmlNewProp(root_node, BAD_CAST "optional", BAD_CAST opt);

  // overview
  xmlNodePtr node = xmlNewChild(root_node, 0, BAD_CAST "overview", 0); 
  // authors
  xmlNodePtr tmp = xmlNewChild(node, 0, BAD_CAST "authors", 0); 
  std::vector<std::string>::const_iterator aiter = mi.authors_.begin();
  while (aiter != mi.authors_.end()) {
    const std::string &a = *aiter++;
    xmlNewChild(tmp, 0, BAD_CAST "author", BAD_CAST a.c_str()); 
  }
  // summary
  xmlNewChild(node, 0, BAD_CAST "summary", BAD_CAST mi.summary_.c_str());

  // io
  node = xmlNewChild(root_node, 0, BAD_CAST "io", 0); 

  tmp = xmlNewChild(node, 0, BAD_CAST "inputs", 0); 
  const char* lpd = mi.last_port_dynamic_ ? "yes" : "no";
  xmlNewProp(tmp, BAD_CAST "lastportdynamic", BAD_CAST lpd);
  std::vector<IPortInfo*>::const_iterator iter = mi.iports_.begin();
  while(iter != mi.iports_.end()) 
  {
    IPortInfo *p = *iter++;
    xmlNodePtr port = xmlNewChild(tmp, 0, BAD_CAST "port", 0); 
    xmlNewChild(port, 0, BAD_CAST "name", BAD_CAST p->name.c_str()); 
    xmlNewChild(port, 0, BAD_CAST "datatype", 
		BAD_CAST p->datatype.c_str()); 
  }

  if (mi.oports_.size()) {
    tmp = xmlNewChild(node, 0, BAD_CAST "outputs", 0); 
  
    std::vector<OPortInfo*>::const_iterator iter = mi.oports_.begin();
    while(iter != mi.oports_.end()) 
    {
      OPortInfo *p = *iter++;
      xmlNodePtr port = xmlNewChild(tmp, 0, BAD_CAST "port", 0); 
      xmlNewChild(port, 0, BAD_CAST "name", BAD_CAST p->name.c_str()); 
      xmlNewChild(port, 0, BAD_CAST "datatype", 
		  BAD_CAST p->datatype.c_str()); 
    }
  }

  // write the file
  xmlSaveFormatFileEnc(filename, doc, "UTF-8", 1);
  
  // free the document
  xmlFreeDoc(doc);
}

bool
read_component_file(ModuleInfo &mi, const char* filename) 
{
  /*
   * this initialize the library and check potential ABI mismatches
   * between the version it was compiled for and the actual shared
   * library used.
   */
  LIBXML_TEST_VERSION;
  
  xmlParserCtxtPtr ctxt; /* the parser context */
  xmlDocPtr doc; /* the resulting document tree */

  /* create a parser context */
  ctxt = xmlNewParserCtxt();
  if (ctxt == 0) {
    std::cerr << "ComponentNode.cc: Failed to allocate parser context\n";
    return false;
  }
  /* parse the file, activating the DTD validation option */
  doc = xmlCtxtReadFile(ctxt, filename, 0, (XML_PARSE_DTDATTR | 
					    XML_PARSE_NOERROR));
  /* check if parsing suceeded */
  if (doc == 0 || ctxt->valid == 0) {

    xmlError* error = xmlCtxtGetLastError(ctxt);
    std::string mtype = "Parse ";
    if (doc) {
      mtype = "Validation ";
    }
    ostringstream msg;
    msg << "createSciDialog -error -title \"Component XML " 
        << mtype << "Error\" -message \"" 
        << endl << mtype << "Failure for: " << endl << filename << endl
        << endl << "Error Message: " << endl << error->message << endl << "\"";
    TCLInterface::eval(msg.str());
    return false;
  } 
  
  xmlNode* node = doc->children;
  for (; node != 0; node = node->next) 
  {
    // skip all but the component node.
    if (node->type == XML_ELEMENT_NODE && 
        std::string(to_char_ptr(node->name)) == std::string("component")) 
    {
      //! set all the ModuleInfo
      xmlAttrPtr name_att = get_attribute_by_name(node, "name");
      xmlAttrPtr cat_att = get_attribute_by_name(node, "category");
      xmlAttrPtr version_att = get_attribute_by_name(node, "version");
            
      // set the component attributes.
      if (name_att == 0 || cat_att == 0) 
      {
        std::cerr << "Attibutes missing from component node in : " 
            << filename << std::endl;
        return false;
      }
      mi.module_name_ = std::string(to_char_ptr(name_att->children->content));
      mi.category_name_ = std::string(to_char_ptr(cat_att->children->content));
      xmlAttrPtr opt_att = get_attribute_by_name(node, "optional");
 
      mi.module_version_ = "1.0";
      if (version_att != 0) 
      {
        mi.module_version_ = std::string(to_char_ptr(version_att->children->content));
      }
      
      mi.optional_ = false;
      if (opt_att != 0) 
      {
        if (std::string(to_char_ptr(opt_att->children->content)) == std::string("true")) 
        {
          mi.optional_ = true;
        }
      }
			
      if (!sci_getenv_p("SCIRUN_SHOW_HIDDEN_MODULES"))
      {
        opt_att = get_attribute_by_name(node, "hide");			
        mi.hide_ = false;
        if (opt_att != 0) 
        {
          if (std::string(to_char_ptr(opt_att->children->content)) == 
              std::string("true")) 
          {
            mi.hide_ = true;
          }
        }
      }

      mi.dynamic_ = false;

      //set the description std::string.
      set_description(mi, node);
      set_port_info(mi, node);
      set_gui_info(mi, node);
    }
  }

  xmlFreeDoc(doc);
  /* free up the parser context */
  xmlFreeParserCtxt(ctxt);  
#ifndef _WIN32
  // there is a problem on windows when using Uintah 
  // which is either caused or exploited by this
  xmlCleanupParser();
#endif
  return true;
}

} // End namespace SCIRun

