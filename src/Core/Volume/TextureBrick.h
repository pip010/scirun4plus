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
//    File   : TextureBrick.h
//    Author : Milan Ikits
//    Date   : Wed Jul 14 15:55:55 2004

#ifndef Volume_TextureBrick_h
#define Volume_TextureBrick_h


#include <slivr/TextureBrick.h>
#include <Core/Geometry/BBox.h>
#include <Core/Thread/Mutex.h>
#include <Core/Containers/LockingHandle.h>
#include <Core/Util/TypeDescription.h>
#include <Core/Datatypes/NrrdData.h>
#include <vector>

#include <Core/Volume/share.h>

namespace SCIRun {

class SCISHARE TextureBrick : public SLIVR::TextureBrick
{
public:
  TextureBrick(NrrdDataHandle n0, NrrdDataHandle n1,
	       int nx, int ny, int nz, int nc, int *nb,
	       int ox, int oy, int oz, int mx, int my, int mz,
	       const BBox& bbox, const BBox& tbox);

  virtual ~TextureBrick();

  enum tb_td_info_e {
    FULL_TD_E,
    TB_NAME_ONLY_E,
    DATA_TD_E
  };

  static const std::string type_name(int n = -1);

  //! needed for our smart pointers -- LockingHandle<T>
  int ref_cnt;
  Mutex lock;
};

typedef LockingHandle<SCIRun::TextureBrick> TextureBrickHandle;

} // namespace SCIRun

#endif // Volume_TextureBrick_h
