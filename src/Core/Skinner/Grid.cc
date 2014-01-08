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
//    File   : Grid.cc
//    Author : McKay Davis
//    Date   : Tue Jun 27 13:03:07 2006

#include <Core/Skinner/Variables.h>
#include <Core/Skinner/Grid.h>
#include <Core/Math/MiscMath.h>
#include <Core/Math/MinMax.h>
#include <Core/Util/Assert.h>
#include <Core/Util/StringUtil.h>
#include <Core/Events/EventManager.h>
#include <iostream>


namespace SCIRun {
namespace Skinner {


Grid::Grid(Variables *variables) :
  Parent(variables),
  rows_(variables,"rows", 1),
  cols_(variables,"cols", 1),
  cell_info_(),
  cell_frac_width_(cols_, AIR_NEG_INF),
  cell_frac_height_(rows_, AIR_NEG_INF),
  resize_mx_(0.0),
  resize_my_(0.0),
  resize_row_(0),
  resize_col_(0)
{
  REGISTER_CATCHER_TARGET(Grid::do_PointerEvent);
}


Grid::~Grid() {
  // Parent class deletes children for us
}


BaseTool::propagation_state_e
Grid::process_event(event_handle_t &event)
{
  if (!visible_()) return STOP_E;

  if (event->is_redraw_event()) {
    ReLayoutCells(event);

    const RectRegion &region = get_region();
    const int rows = rows_;
    const int cols = cols_();
    const double wid = region.width();
    const double hei = region.height();

    vector<double> posx(cols+1, 0);
    for (int i = 1; i <= cols; ++i) {
      double dx = Max(cell_frac_width_[i-1] * wid, 0.0);
      posx[i] = Round(Min(posx[i-1] + dx, wid));
    }

    vector<double> posy(rows+1, 0);
    for (int i = 1; i <= rows; ++i) {
      double dy = Max(cell_frac_height_[i-1] * hei, 0.0);
      posy[i] = Round(Min(posy[i-1] + dy, hei));
    }

    for (unsigned int i = 0; i < cell_info_.size(); ++i) {
      CellInfo_t &cell = cell_info_[i];
      if (cell.row_begin_.exists() && cell.row_end_.exists() &&
          cell.col_begin_.exists() && cell.col_end_.exists()) {
        children_[i]->set_region
          (RectRegion(region.x1() + posx[cell.col_begin_()-1],
                      region.y2() - posy[cell.row_end_()],
                      region.x1() + posx[cell.col_end_()],
                      region.y2() - posy[cell.row_begin_()-1]));
      } else {
        int r = cell.row_()-1;
        int c = cell.col_()-1;

        children_[i]->set_region(RectRegion(region.x1() + posx[c],
                                            region.y2() - posy[r + 1],
                                            region.x1() + posx[c + 1],
                                            region.y2() - posy[r]));
      }
    }
  }

  for (unsigned int i = 0; i < cell_info_.size(); ++i) {
    children_[i]->process_event(event);
  }

  return Drawable::process_event(event);
}


void
Grid::set_children(const Drawables_t &children)
{
  children_ = children;
  cell_info_.resize(children_.size());
  for (unsigned int i = 0; i < children_.size(); ++i) {
    Variables *cvars = children_[i]->get_vars();
    CellInfo_t &cell = cell_info_[i];
    cell.row_ = Var<int>(cvars, "row", 1);
    cell.col_ = Var<int>(cvars, "col", 1);
    cell.width_ = Var<double>(cvars, "cell-width", -1.0);
    cell.height_ = Var<double>(cvars, "cell-height", -1.0);
    cell.col_begin_ = Var<int>(cvars, "col-begin");
    cell.col_end_ = Var<int>(cvars, "col-end");
    cell.row_begin_ = Var<int>(cvars, "row-begin");
    cell.row_end_ = Var<int>(cvars, "row-end");
    cell.resize_horizontal_ = Var<bool>(cvars, "resize-horizontal", 0);
    cell.resize_vertical_ = Var<bool>(cvars, "resize-vertical", 0);
  }
}


BaseTool::propagation_state_e
Grid::ReLayoutCells(event_handle_t &)
{
  vector<double>cell_frac_width(cols_, AIR_NEG_INF);
  vector<double>cell_frac_height(rows_, AIR_NEG_INF);

  const double wid = get_region().width();
  const double hei = get_region().height();

  for (unsigned int i = 0; i < cell_info_.size(); ++i) {
    CellInfo_t &cell = cell_info_[i];
    //        if (!cell.row_.exists() || !cell.col_.exists()) continue;
    const int row = Clamp(cell.row_ - 1, 0, rows_-1);
    const int col = Clamp(cell.col_ - 1, 0, cols_-1);

    if (cell.width_ >= 0.0) {
      double frac_width = (cell.width_<1.0)?cell.width_:(cell.width_/wid);
      cell_frac_width[col] = Max(cell_frac_width[col], frac_width);
    }

    if (cell.height_ >= 0.0) {
      double frac_height=(cell.height_<1.0)?cell.height_:(cell.height_/hei);
      cell_frac_height[row] = Max(cell_frac_height[row], frac_height);
    }
  }

  double total = 0.0;
  int count = 0;
  for (int i = 0; i < cols_; i++) {
    if (cell_frac_width[i] >= 0.0) {
      total += cell_frac_width[i];
      count++;
    }
  }

  if (count < cols_) {
    double left_over = Max(0.0, (1.0-total)/(cols_ - double(count)));
    for (int i = 0; i < cols_; i++) {
      if (cell_frac_width[i] < 0.0) {
        cell_frac_width[i] = left_over;
      }
    }
  }

  total = 0.0;
  count = 0;
  for (int i = 0; i < rows_; i++) {
    if (cell_frac_height[i] >= 0.0) {
      total += cell_frac_height[i];
      count++;
    }
  }

  if (count < rows_) {
    double left_over = Max(0.0, (1.0-total)/(rows_ - double(count)));
    for (int i = 0; i < rows_; i++) {
      if (cell_frac_height[i] < 0.0) {
        cell_frac_height[i] = left_over;
      }
    }
  }

  cell_frac_width_ = cell_frac_width;
  cell_frac_height_ = cell_frac_height;

  return STOP_E;
}


BaseTool::propagation_state_e
Grid::do_PointerEvent(event_handle_t &event)
{
  ASSERT(dynamic_cast<PointerSignal *>(event.get_rep()));
  PointerSignal *signal = (PointerSignal *)(event.get_rep());
  PointerEvent *pointer = signal->get_pointer_event();

  bool b1 = pointer->get_pointer_state() & PointerEvent::BUTTON_1_E;
  if (!b1) return CONTINUE_E;


  bool rel = pointer->get_pointer_state() & PointerEvent::BUTTON_RELEASE_E;

  if (rel) {
    resize_row_ = 0;
    resize_col_ = 0;
    return CONTINUE_E;
  }

  bool press = pointer->get_pointer_state() & PointerEvent::BUTTON_PRESS_E;

  const RectRegion &reg = get_region();

  const double mx = (pointer->get_x()-reg.x1())/reg.width();
  const double my = (pointer->get_y()-reg.y1())/reg.height();

  if (press) {
    resize_mx_ = mx;
    resize_my_ = my;
    for (unsigned int i = 0; i < cell_info_.size(); ++i) {
      RectRegion creg = children_[i]->get_region();
      creg = RectRegion(creg.x1()-2,creg.y1()-2, creg.x2()+2, creg.y2()+2);
      if (!creg.inside(pointer->get_x(), pointer->get_y())) continue;
      CellInfo_t &cell = cell_info_[i];
      const int row = cell.row_ -1;
      if (cell.resize_vertical_ && row > 0 && row < rows_-1)
      {
        resize_row_ = row;
        resize_vertical_bounds_ =
          make_pair(cell_frac_height_[row-1], cell_frac_height_[row+1]);
      }

      const int col = cell.col_ -1;
      if (cell.resize_horizontal_ && col > 0 && col < cols_-1)
      {
        resize_col_ = col;
        resize_horizontal_bounds_ =
          make_pair(cell_frac_width_[col-1], cell_frac_width_[col+1]);
      }
    }
  }
  else
  {
    bool redraw = false;

    double dx = -mx + resize_mx_;
    double dy = my - resize_my_;

    for (unsigned int i = 0; i < cell_info_.size(); ++i) {
      CellInfo_t &cell = cell_info_[i];

      if (resize_row_) {
        const int row = cell.row_ -1;
        const double max = (resize_vertical_bounds_.first +
                            resize_vertical_bounds_.second);
        if (row == resize_row_ - 1) {
          cell.height_ = Clamp(resize_vertical_bounds_.first - dy,0., max);
          redraw = true;
        } else if (row == resize_row_ + 1) {
          cell.height_ = Clamp(resize_vertical_bounds_.second + dy,0.,max);
          redraw = true;
        }
      }

      if (resize_col_) {
        const int col = cell.col_ -1;
        const double max = (resize_horizontal_bounds_.first +
                            resize_horizontal_bounds_.second);

        if (col == resize_col_ - 1) {
          cell.width_ = Clamp(resize_horizontal_bounds_.first-dx,0.,max);
          redraw = true;
        } else if (col == resize_col_ + 1) {
          cell.width_ = Clamp(resize_horizontal_bounds_.second+dx,0.,max);
          redraw = true;
        }
      }
    }

    if (redraw) {
      throw_signal("Grid::resized");
    }
  }

  return CONTINUE_E;
}


int
Grid::get_signal_id(const string &signalname) const
{
  if (signalname == "Grid::resized") return 1;
  return 0;
}


}
}
