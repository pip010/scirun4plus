//
//  For more information, please see: http://software.sci.utah.edu
//
//  The MIT License
//
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//
//
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
//
//    File   : XMLIO.cc
//    Author : McKay Davis
//    Date   : Tue Jun 27 13:01:28 2006

#include <Core/Skinner/Skinner.h>
#include <Core/Skinner/Variables.h>
#include <Core/Skinner/XMLIO.h>

#include <Core/XMLUtil/XMLUtil.h>
#include <Core/Util/StringUtil.h>
#include <Core/Util/Environment.h>
#include <Core/Util/Assert.h>

#include <libxml/xmlreader.h>
#include <libxml/catalog.h>
#include <libxml/xinclude.h>
#include <libxml/xpathInternals.h>

#include <iostream>

#include <string.h>
#include <stdio.h>

#include <sci_debug.h>

namespace SCIRun {
namespace Skinner {


static bool debug_xmlio = false;
XMLIO::DrawableMakerMap_t XMLIO::makers_;


XMLIO::XMLIO()
{
}


XMLIO::~XMLIO()
{
}


Skinner *
XMLIO::load(const string &filename)
{
  debug_xmlio = sci_getenv_p("SKINNER_XMLIO_DEBUG");
  /*
   * this initialize the library and check potential ABI mismatches
   * between the version it was compiled for and the actual shared
   * library used.
   */

  LIBXML_TEST_VERSION;

  xmlParserCtxtPtr ctxt; /* the parser context */

  std::string dtd = std::string(sci_getenv("SCIRUN_SRCDIR")) +
    std::string("/Core/Skinner/skinner.dtd");

  xmlInitializeCatalog();
  xmlCatalogAdd(XMLUtil::char_to_xmlChar("public"),
                XMLUtil::char_to_xmlChar("-//Skinner/Drawable DTD"),
                XMLUtil::char_to_xmlChar(dtd.c_str()));

  /* create a parser context */
  ctxt = xmlNewParserCtxt();
  if (!ctxt) {
    std::cerr << "XMLIO::load failed xmlNewParserCtx()\n";
    return 0;
  }

  int flags = XML_PARSE_DTDATTR | XML_PARSE_DTDVALID | XML_PARSE_PEDANTIC;

  /* parse the file, activating the DTD validation option */
  xmlDocPtr doc = xmlCtxtReadFile(ctxt, filename.c_str(), 0, flags);

#ifdef DEBUG
  std::cerr << "Reading SKINNER xml file : " << filename.c_str() << std::endl;
#endif

  /* apply the XInclude process, this should trigger the I/O just
   * registered. */
  int inc = xmlXIncludeProcess(doc);
  if (inc < 0) 
  {
    std::cerr << "XInclude processing failed\n";
    return 0;
  }

  xmlCtxtUseOptions(ctxt, flags);

  /* check if parsing suceeded */
  if (!doc) 
  {
    std::cerr << "Skinner::XMLIO::load failed to parse "
              << filename << std::endl;
    return 0;
  }


  xmlValidCtxtPtr valid_ctx = xmlNewValidCtxt();
  xmlValidateDtdFinal(valid_ctx, doc);

  // parse the doc at network node.
  Skinner *skinner = 0;
  for (xmlNode *cnode=doc->children; cnode!=0; cnode=cnode->next) 
  {
    if (XMLUtil::node_is_dtd(cnode, "skinner")) 
    {
      continue;
    } 
    else if (XMLUtil::node_is_element(cnode, "skinner")) 
    {
      skinner = eval_skinner_node(cnode);
    } 
    else if (!XMLUtil::node_is_comment(cnode)) 
    {
      throw "Unknown node type";
    }
  }

  const char *output = sci_getenv("SKINNER_OUTPUT_FILE");
  if (output) 
  {
    xmlKeepBlanksDefault(0);
    xmlIndentTreeOutput = 1;
    xmlSaveFormatFile(output, doc, 1);
  }

  xmlFreeDoc(doc);

  /* free up the parser context */
  xmlFreeParserCtxt(ctxt);
#ifndef _WIN32
  // there is some problem in the windows version which is
  // either caused or exploited by this - LIBXML_TEST_VERSION
  // will fail when called from a different module
  xmlCleanupParser();
#endif

  return skinner;
}


xmlNodePtr
XMLIO::find_definition(definition_nodes_t &definitions,
                       const char *classname)
{
  // Go backwards through the vector becasue we want to search
  // UP the current object context tree defined in the skinner file
  // looking for definitions in the nearest ancestor before looking
  // at their parent node
  definition_nodes_t::reverse_iterator def_map_iter = definitions.rbegin();
  definition_nodes_t::reverse_iterator def_map_end = definitions.rend();
  for (;def_map_iter != def_map_end; ++def_map_iter) {
    string_node_map_t::iterator found_def = def_map_iter->find(classname);
    if (found_def != def_map_iter->end()) {
      return found_def->second;
    }
  }
  return 0;
}


void
XMLIO::eval_merged_object_nodes_and_push_definitions
(merged_nodes_t &merged_nodes, definition_nodes_t &defs)
{
  defs.push_back(string_node_map_t());
  string_node_map_t &node_def_map = defs.back();
  // Iterate through the xml var nodes and set their values
  for (merged_nodes_t::iterator mnode = merged_nodes.begin();
       mnode != merged_nodes.end(); ++mnode) {
    for (xmlNode *cnode = (*mnode)->children; cnode; cnode = cnode->next) {
      if (XMLUtil::node_is_element(cnode, "definition")) {
        eval_definition_node(cnode, node_def_map);
      }
    }
  }
}


Drawable *
XMLIO::eval_object_node(const xmlNodePtr node,
                        Drawable *parent,
                        definition_nodes_t &definitions)
{
  const char *classname = 0;
  std::string unique_id = "";
  try 
  {
    classname = XMLUtil::maybe_get_att_as_const_char_str(node, "name");
    if (!classname) 
    {
      classname = "Skinner::Skinner";
    }

    // Object id is not required, create unique one if not found
    // The first merged node contains the id

    if (!XMLUtil::maybe_get_att_as_string(node, "id",  unique_id)) 
    {
      // Note: This isnt absolutely guaranteed to make unique ids
      // It can be broken by hardcoding the id in the XML
      static int objcount = 0;
      unique_id = std::string(classname)+"-"+to_string(objcount++);
    }

    // First, before we can even construct the object, we need to
    // create the Variables that determine this instance's unique properties
    // So, create a new Variables context with our Unique ID
    // Note: this creates new memory that needs to be freed by the object
    Variables *parent_vars = parent ? parent->get_vars() : 0;
    Variables *variables = new Variables(unique_id, parent_vars);

    // The classname could refer to a previously parsed <definition> node
    // in that case, we instance the contained single object node as it were
    // created in the current context.  <var>s and <signals> are tags
    // are allowed in both this <object> tag and the <definition>
    // encapsulated <object> tag.  The Variables, signals/catchers of the
    // encapsulated tag are merged last.
    // Is this object tag, just an alias to another object node?
    // Definitions can contain nested already decalred definitions as
    // their object classname, unroll the definitions  until we find
    // an non previously defined classname
    merged_nodes_t merged_nodes(1,node);
    
    xmlNodePtr dnode = find_definition(definitions, classname);
    while (dnode)
    {
      // When searching for vars, signals, and children
      // we need to merge the encapsulated nodes
      merged_nodes.push_back(dnode);
      // If we are just a reference to a definition, we need to switch
      // classnames to the encapsulated object
      classname = XMLUtil::node_att_as_const_char_str(dnode, "name");

#ifdef DEBUG
      std::cerr << variables->get_id()
           << " is actually classname: "
           << classname << "\n";
#endif

      // Iteratre to the next nested definition, if there is one...
      dnode = find_definition(definitions, classname);
    }

    // Iterate through the xml <var> nodes and record the values
    for (merged_nodes_t::iterator mnode = merged_nodes.begin();
         mnode != merged_nodes.end(); ++mnode) 
    {
      xmlNodePtr node = *mnode;
      eval_vars(node, variables);
    }


    // Set the class variable before constructing the class
    std::vector<std::string> classnames = split_string(classname, ':');
    variables->insert("name", classnames.back(), "string", false);

    // Now that we have the Variables set up doen, time to create the Object
    Drawable * object = 0;

    if (!parent) 
    {
      object = new Skinner(variables);
    } 
    else 
    {
      // First, see if the current catchers can create and throw back
      // an object of type "classname"
      std::string makerstr = std::string(classname)+"_Maker";

      // The special Signal to ask for a maker is created
      MakerSignal *makersig = new MakerSignal(makerstr, variables);
      event_handle_t makerevent = makersig;

      // And thrown to the Catcher...
      // Parent is guaranteed to be non-null here.
      parent->throw_signal_extended(makerevent);

      ASSERT(dynamic_cast<MakerSignal*>(makerevent.get_rep()));

      // And we see what the cather returned
      makersig = (MakerSignal *)(makerevent.get_rep());
      if (makersig->get_signal_name() == (makerstr+"_Done")) 
      {
        // MakerSignal returned with success, so
        // The MakerSignal thrower is the newly constructed class
        object = dynamic_cast<Drawable *>(makersig->get_signal_thrower());
      } 
      else 
      {
        // Search the static_makers table and see if it contains
        // the wanted classname
        if (makers_.find(classname) != makers_.end()) 
        {
          object = (*makers_[classname])(variables);
        }
      }
    }

    // At this point, the if the object is uninstantiatable, return
    if (!object) 
    {
      std::cerr << variables->get_id()
           << " - Skinner::XMLIO::eval_object_node - object class=\""
           << classname << "\" cannot find maker\n";
      delete variables;
      return 0;
    }


    object->set_parent(parent);


    // Now the Catchers On Deck are ready, look if the xml file has
    // created any signals to hookup
    for (merged_nodes_t::iterator mnode = merged_nodes.begin();
         mnode != merged_nodes.end(); ++mnode) 
    {
      for (xmlNode *cnode = (*mnode)->children; cnode; cnode = cnode->next) 
      {
        if (XMLUtil::node_is_element(cnode, "signal")) 
        {
          eval_signal_node(cnode, object);
        }

        if (XMLUtil::node_is_element(cnode, "bridge")) 
        {
          eval_bridge_node(cnode, object);
        }

      }
    }

    eval_merged_object_nodes_and_push_definitions(merged_nodes, definitions);

    // Time to look for children object nodes
    Drawables_t children(0);
    bool have_unwanted_children = false;
    Parent *isa_parent = dynamic_cast<Parent *>(object);
    // Search the merged nodes for object nodes.
    for (merged_nodes_t::iterator mnode = merged_nodes.begin();
         mnode != merged_nodes.end(); ++mnode) 
    {
      for (xmlNode *cnode = (*mnode)->children; cnode; cnode = cnode->next) 
      {
        if (XMLUtil::node_is_element(cnode, "object")) 
        {
          if (isa_parent) 
          {
            Drawable *child =
              eval_object_node(cnode, object,definitions);
            if (child) 
            {
              children.push_back(child);
            }
          } 
          else 
          {
            have_unwanted_children = true;
          }
        }
      }
    }

    if (isa_parent && children.size()) 
    {
      isa_parent->set_children(children);
    }

    if (have_unwanted_children) 
    {
      std::cerr << "class : " << classname << " does not allow <object>\n";
    }

    definitions.pop_back();

    return object;

  } 
  catch (...) 
  {
    std::cerr << "Error evaluating object: "
         << (classname ? classname : "UNKNOWN") << std::endl;
    std::cerr << "id: " << unique_id << std::endl;
    return 0;
  }
}


Skinner *
XMLIO::eval_skinner_node(const xmlNodePtr node)
{
  ASSERT(XMLUtil::node_is_element(node, "skinner"));
  definition_nodes_t definitions;
  Drawable * object = eval_object_node(node, 0, definitions);
  Skinner *skinner = 0;
  if (object) 
  {
    ASSERT(dynamic_cast<Skinner*>(object));
    skinner = (Skinner *)object;
  }

  return skinner;
}


void
XMLIO::eval_definition_node(const xmlNodePtr node,
                            string_node_map_t &definitions)
{
  ASSERT(XMLUtil::node_is_element(node, "definition"));
  const char *name = XMLUtil::node_att_as_const_char_str(node, "name");
  for (xmlNode *cnode=node->children; cnode!=0; cnode=cnode->next) 
  {
    if (XMLUtil::node_is_element(cnode, "object")) 
    {
      definitions[name] = cnode;
    }
  }
}


void
XMLIO::eval_vars(const xmlNodePtr &node, Variables *variables)
{
  for (xmlNode *cnode = node->children; cnode; cnode = cnode->next) {
    if (XMLUtil::node_is_element(cnode, "var")) {
      eval_var_node(cnode, variables);
    }
  }
}


void
XMLIO::eval_var_node(const xmlNodePtr node,
                     Variables *variables, bool force)
{
  ASSERT(XMLUtil::node_is_element(node, "var"));
  ASSERT(variables);
  const char *varname = XMLUtil::node_att_as_const_char_str(node, "name");

  const char *overwrite_str =
    XMLUtil::maybe_get_att_as_const_char_str(node, "overwrite");

  if (overwrite_str && !strcmp(overwrite_str,"no") &&
      variables->exists(varname)) {
    return;
  }

  bool propagate = force;
  if (!propagate) {
    const char *propagate_str =
      XMLUtil::maybe_get_att_as_const_char_str(node, "propagate");
    propagate = (propagate_str && !strcmp(propagate_str,"yes"));
  }

  const char *value_str = 0;
  if (node->children) {
    value_str = XMLUtil::xmlChar_to_char(node->children->content);
  }
  const string value(value_str ? value_str : "");

  const char *type_str =
    XMLUtil::maybe_get_att_as_const_char_str(node, "type");
  string typ(type_str ? type_str : "unknown");

  variables->insert(varname, value, typ, propagate);
}


void
XMLIO::eval_signal_node(const xmlNodePtr node,
                        Drawable *object)
{
  ASSERT(XMLUtil::node_is_element(node, "signal"));
  string signalname = XMLUtil::node_att_as_string(node, "name");
  string signaltarget = XMLUtil::node_att_as_string(node, "target");

  Variables *vars = new Variables();
  for (xmlNode *cnode = node->children; cnode; cnode = cnode->next) {
    if (XMLUtil::node_is_element(cnode, "var")) {
      eval_var_node(cnode, vars, 1);
    }
  }

  object->add_callback(signalname, signaltarget, vars);
}


void
XMLIO::eval_bridge_node(const xmlNodePtr node,
                        Drawable *object)
{
  ASSERT(XMLUtil::node_is_element(node, "bridge"));
  string name = XMLUtil::node_att_as_string(node, "name");
  const char *target =
    XMLUtil::maybe_get_att_as_const_char_str(node, "target");

  if (!target) {
    object->make_bridge(name);
  } else {
    object->use_bridge(name, object, target);
  }
}


void
XMLIO::register_maker(const string &name, DrawableMakerFunc_t *maker)
{
  makers_[name] = maker;
}


}
}
