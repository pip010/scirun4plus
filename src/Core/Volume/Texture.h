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
//    File   : Texture.h
//    Author : Milan Ikits
//    Date   : Thu Jul 15 01:00:36 2004

#ifndef Volume_Texture_h
#define Volume_Texture_h

#include <slivr/Texture.h>
#include <Core/Datatypes/Datatype.h>
#include <Core/Datatypes/NrrdData.h>
#include <Core/Persistent/Persistent.h>
#include <Core/Containers/LockingHandle.h>
#include <Core/Volume/share.h>

namespace SCIRun {

class SCISHARE Texture : public Datatype, public SLIVR::Texture
{
public:
  Texture();
  virtual ~Texture();
 
  virtual void io(Piostream&);
  static PersistentTypeID type_id;
  virtual std::string dynamic_type_name() const { return type_id.type; }
  
  // As SLIVR does not do memory management
  // we need to keep an handle so memory will not
  // be destroyed before the object is removed
  NrrdDataHandle value_nrrd_;
  NrrdDataHandle gradient_magnitude_nrrd_;
  NrrdDataHandle mask_nrrd_;
};

typedef LockingHandle<Texture> TextureHandle;

} // namespace SCIRun

#endif // Volume_Texture_h
