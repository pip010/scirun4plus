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
 *  IndexedGroup.cc: ?
 *
 *  Written by:
 *   Author: ?
 *   Department of Computer Science
 *   University of Utah
 *   Date: ?
 *
 */

#include <Core/Util/Debug.h>

#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/IndexedGroup.h>

namespace SCIRun {

//----------------------------------------------------------------------
GeomIndexedGroup::GeomIndexedGroup()
{
  DEBUG_CONSTRUCTOR("GeomIndexedGroup")
}

//----------------------------------------------------------------------
GeomIndexedGroup::~GeomIndexedGroup()
{
  delAll();  // just nuke everything for now...
  DEBUG_DESTRUCTOR("GeomIndexedGroup")
}

//----------------------------------------------------------------------
GeomObj* 
GeomIndexedGroup::clone()
{
  return new GeomIndexedGroup(*this);
}

//----------------------------------------------------------------------
void 
GeomIndexedGroup::reset_bbox()
{
  MapIntGeomObj::iterator itr = objs.begin();
  while (itr != objs.end())
  {
    (*itr).second->reset_bbox();
    ++itr;
  }
}

//----------------------------------------------------------------------
void 
GeomIndexedGroup::get_bounds(BBox& bbox)
{
  MapIntGeomObj::iterator iter;
  for (iter = objs.begin(); iter != objs.end(); iter++) 
  {
    (*iter).second->get_bounds(bbox);
  }  
}

void
GeomIndexedGroup::draw(DrawInfoOpenGL* di, Material* m, double time)
{
  MapIntGeomObj::iterator iter;
  for (iter = objs.begin(); iter != objs.end(); iter++)
  {
    (*iter).second->draw(di, m, time);
  }
}

//----------------------------------------------------------------------
void 
GeomIndexedGroup::addObj(GeomHandle obj, int id)
{
  objs[id] = obj;
}

//----------------------------------------------------------------------
GeomHandle 
GeomIndexedGroup::getObj(int id)
{
  MapIntGeomObj::iterator iter = objs.find(id);
  if (iter != objs.end()) 
  {
    return (*iter).second;
  }
  else 
  {
    //    cerr << "couldn't find object in GeomIndexedGroup::getObj!\n";
  }
  return 0;
}

//----------------------------------------------------------------------
void 
GeomIndexedGroup::delObj(int id)
{
  MapIntGeomObj::iterator iter = objs.find(id);
  if (iter != objs.end()) 
  {
    objs.erase(iter);
  }
  else 
  {
    //   cerr << "invalid id in GeomIndexedGroup::delObj()!\n";
  }
}

//----------------------------------------------------------------------
void 
GeomIndexedGroup::delAll()
{
  objs.clear();
}

//----------------------------------------------------------------------
GeomIndexedGroup::IterIntGeomObj 
GeomIndexedGroup::getIter(void)
{
  return IterIntGeomObj(objs.begin(), objs.end());
}

//----------------------------------------------------------------------
GeomIndexedGroup::MapIntGeomObj* 
GeomIndexedGroup::getHash(void)
{
  return &objs;
}

} // End namespace SCIRun


