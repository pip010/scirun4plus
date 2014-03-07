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
//    File   : Undo.h
//    Author : Michael Callahan
//    Date   : June 2008

#ifndef SEG3D_Undo_H
#define SEG3D_Undo_H

#include <Applications/Seg3D/NrrdVolume.h>

namespace SCIRun {

using std::list;

class Painter;

class Undo
{
public:
  Undo(Painter *painter, const string &label);
  virtual ~Undo();

  virtual void undo() = 0;
  virtual int mem_est() = 0;

  string get_label() { return label_; }

  // For LockingHandle<NrrdVolume> NrrdVolumeHandle typedef after this class
  ThreadLock          lock;
  unsigned int        ref_cnt;

protected:
  Painter *painter_;
  string label_;
};

typedef LockingHandle<Undo> UndoHandle;


class UndoManager
{
public:
  UndoManager();
  ~UndoManager();

  void push_undo(const UndoHandle &undo_obj);
  void undo_last();
  void clear_undo();

protected:
  list<UndoHandle> undo_stack_;

  void update_gui();
};


class UndoReplaceLayer : public Undo
{
public:
  UndoReplaceLayer(Painter *painter,
                   const string &label,
                   const NrrdVolumeHandle &oldvol,
                   const NrrdVolumeHandle &newvol,
                   size_t loc);

  virtual ~UndoReplaceLayer();

  virtual void undo();
  virtual int mem_est();

private:
  NrrdVolumeHandle oldvol_;
  NrrdVolumeHandle newvol_;
  size_t loc_;
};


class UndoLabelInvertFilter : public Undo
{
public:
  UndoLabelInvertFilter(Painter *painter,
                        const NrrdVolumeHandle &volume);

  virtual ~UndoLabelInvertFilter();

  virtual void undo();
  virtual int mem_est();

private:
  NrrdVolumeHandle volume_;
};


class UndoReplaceSlice : public Undo
{
public:
  UndoReplaceSlice(Painter *painter,
                   const string &label,
                   NrrdVolumeHandle &vol,
                   NrrdDataHandle &slice,
                   int axis, int coord);

  virtual ~UndoReplaceSlice();
  
  // Takes in stuff required to fill in the nrrdsplice.

  virtual void undo();
  virtual int mem_est();

private:
  NrrdVolumeHandle volume_;
  NrrdDataHandle slice_;
  int axis_;
  int coord_;
};


class UndoReplaceVolume : public Undo
{
public:
  UndoReplaceVolume(Painter *painter,
                    const string &label,
                    NrrdVolumeHandle &volume,
                    NrrdDataHandle &data);

  virtual ~UndoReplaceVolume();

  virtual void undo();
  virtual int mem_est();

private:
  NrrdVolumeHandle volume_;
  NrrdDataHandle data_;
};


class UndoFlipAxis : public Undo
{
public:
  UndoFlipAxis(Painter *painter,
	       const string &label,
	       NrrdVolumeHandle &volume,
	       int axis);

  virtual ~UndoFlipAxis();

  virtual void undo();
  virtual int mem_est();

private:
  NrrdVolumeHandle volume_;
  int axis_;
};


class UndoInvert : public Undo
{
public:
  UndoInvert(Painter *painter,
	     const string &label,
	     NrrdVolumeHandle &volume);

  virtual ~UndoInvert();

  virtual void undo();
  virtual int mem_est();

private:
  NrrdVolumeHandle volume_;
};


} // namespace SCIRun

#endif
