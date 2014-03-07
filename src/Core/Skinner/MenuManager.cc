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
//    File   : MenuManager.cc
//    Author : McKay Davis
//    Date   : Fri Aug 11 20:45:24 2006

#include <Core/Skinner/MenuList.h>
#include <Core/Skinner/MenuButton.h>
#include <Core/Skinner/MenuManager.h>
#include <Core/Skinner/Variables.h>
#include <Core/Events/EventManager.h>
#include <Core/Geom/FontManager.h>
#include <Core/Util/StringUtil.h>

namespace SCIRun {
namespace Skinner {

MenuManager::MenuManager(Variables *vars) :
  Parent(vars),
  mutex_("MenuManager"),
  menu_lists_(),
  menu_buttons_(),
  visible_menulists_()
{
  REGISTER_CATCHER_TARGET(MenuManager::MenuList_Maker);
  REGISTER_CATCHER_TARGET(MenuManager::MenuButton_Maker);
  REGISTER_CATCHER_TARGET(MenuManager::show_MenuList);
  REGISTER_CATCHER_TARGET(MenuManager::hide_MenuList);
  REGISTER_CATCHER_TARGET(MenuManager::do_PointerEvent);
}


MenuManager::~MenuManager()
{
}


BaseTool::propagation_state_e
MenuManager::MenuList_Maker(event_handle_t &maker_signal)
{
  menu_lists_.push_back
    (construct_child_from_maker_signal<MenuList>(maker_signal));
  menu_lists_.back()->visible_ = false;
  return STOP_E;
}


BaseTool::propagation_state_e
MenuManager::MenuButton_Maker(event_handle_t &maker_signal)
{
  menu_buttons_.push_back
    (construct_child_from_maker_signal<MenuButton>(maker_signal));
  return STOP_E;
}


BaseTool::propagation_state_e
MenuManager::show_MenuList(event_handle_t &event)
{
  Signal *signal = dynamic_cast<Signal *>(event.get_rep());
  Variables *vars = signal->get_vars();
  Var<string> id(vars, "id");
  Var<double> button_x1(vars, "MenuButton::x1");
  Var<double> button_y1(vars, "MenuButton::y1");
  Var<int> window_width(get_vars(), "GLWindow::width");
  if (!id.exists() || id().empty() ||
      !button_x1.exists() || !button_y1.exists()) {
    return STOP_E;
  }
  MenuLists_t::iterator citer = menu_lists_.begin();
  for (; citer != menu_lists_.end(); ++citer) {
    MenuList *menulist = *citer;
    const string menu_id = menulist->get_id();
    if (ends_with(menu_id,id)) {
      double x1 = button_x1;
      if (window_width.exists() && ((x1+menulist->width_)>window_width()))
      {
        x1 = window_width - menulist->width_;
      }

      menulist->x1_ = x1;
      menulist->y1_ = button_y1 - menulist->height_;
      menulist->x2_ = x1 + menulist->width_;
      menulist->y2_ = button_y1;

      menulist->visible_ = true;
      mutex_.lock();
      visible_menulists_.insert(make_pair(menu_id,*citer));
      mutex_.unlock();
      break;
    }
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e
MenuManager::do_PointerEvent(event_handle_t &event)
{
  PointerSignal *signal = dynamic_cast<PointerSignal *>(event.get_rep());
  PointerEvent *pointer = signal->get_pointer_event();
  x_ = pointer->get_x();
  y_ = pointer->get_y();
  return CONTINUE_E;
}


BaseTool::propagation_state_e
MenuManager::hide_MenuList(event_handle_t &event)
{
  Signal *signal = dynamic_cast<Signal *>(event.get_rep());
  if (!signal->get_vars()->exists("id"))
    return CONTINUE_E;

  const string id = signal->get_vars()->get_id();
  Var<bool> ignoreinside(signal->get_vars(),"ignoreinside",0);
  mutex_.lock();
  string erase_id = "";
  VisibleMenuLists_t::iterator citer = visible_menulists_.begin();
  for (; citer != visible_menulists_.end(); ++citer) {
    MenuList *menulist = citer->second;
    string menu_id = menulist->get_id();
    if (ends_with(menu_id,id)) {
      if (ignoreinside && menulist->get_region().inside(x_,y_)) {
        continue;
      }
      menulist->visible_ = false;
      erase_id = menu_id;
      break;
    }
  }

  if (!erase_id.empty()) {
    visible_menulists_.erase(erase_id);
  }

  mutex_.unlock();
  return CONTINUE_E;
}


}
}
