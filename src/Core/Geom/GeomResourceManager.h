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
 *  GeomResourceManager.h.h
 *
 *  Written by:
 *   Michael Callahan
 *   SCI Institute
 *   University of Utah
 *   2008
 *
 */

#ifndef SCIRun_Core_Geom_GeomResourceManager_h
#define SCIRun_Core_Geom_GeomResourceManager_h

#include <Core/Thread/Mutex.h>
#include <vector>
#include <Core/Geom/share.h>

namespace SCIRun {

class SCISHARE GeomResourceManager {
public:
  
  // There is no way to delete textures outside of the context that
  // they are created in.  This means that if you delete them from a
  // GeomObj destructor the program will crash unless a context
  // happens to be active.  Instead just put the texture IDs on a list
  // and check to see if any need to be deleted at the beginning of
  // each draw loop.
  //
  // Note that for this to work delete_pending_objects needs to be added
  // at the top of the draw loop.

  // TODO:  Shader objects should also be deleted this way.
  // TODO:  SLIVR/Volume Rendering probably needs this also.
  // TODO:  delete_pending_objects needs to be added to all the draw loops.

  GeomResourceManager();
  ~GeomResourceManager();

  static void delete_pending_objects();
  static void delete_texture_id(unsigned int id);

private:

  void delete_pending_objects_nonstatic();
  void delete_texture_id_nonstatic(unsigned int id);

  Mutex                     lock_;
  std::vector<unsigned int> texture_ids_;

  static GeomResourceManager manager_;
};

}

#endif
