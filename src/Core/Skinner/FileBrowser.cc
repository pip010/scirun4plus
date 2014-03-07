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
//    File   : FileBrowser.cc
//    Author : McKay Davis
//    Date   : Sun Feb 11 01:55:02 2007

#include <Core/Skinner/FileBrowser.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/MiscMath.h>
#include <Core/Util/FileUtils.h>
#include <Core/Util/StringUtil.h>
#include <Core/Events/keysyms.h>

#include <algorithm> // for sort()

namespace SCIRun {
namespace Skinner {


FileBrowser::FileBrowser(Variables *vars) :
  Parent(vars),
  filename_(vars,"FileBrowser::filename"),
  directory_(vars,"FileBrowser::directory"),
  cached_filename_("INVALIDCACHE")
{
  REGISTER_CATCHER_TARGET(FileBrowser::rescan);
  REGISTER_CATCHER_TARGET(FileBrowser::do_KeyEvent);
  REGISTER_CATCHER_TARGET(FileBrowser::file_selected);
  for (int i = 0; i < 100; ++i) {
    get_vars()->insert("FileBrowser::filename_"+to_string(i),
                       "BLAH"+to_string(i), "string", 1);
  }
  event_handle_t nul = 0;
  rescan(nul);
}


BaseTool::propagation_state_e
FileBrowser::do_KeyEvent(event_handle_t &event)
{
  KeySignal *keysig = dynamic_cast<KeySignal *>(event.get_rep());
  ASSERT(keysig);
  KeyEvent *key = keysig->get_key_event();

  if (key->get_key_state() & KeyEvent::KEY_PRESS_E)
  {
    if (key->get_keyval() == SCIRun_Tab) {
      string candidate = directory_()+"/"+filename_();
      if (filename_().length() && filename_()[0] == '/') {
        candidate = filename_;
      }
      candidate = canonicalize(candidate);
      string after = autocomplete(candidate);
      int diff = after.length() - candidate.length();
      if (diff > 0) {
        filename_ = filename_() + after.substr(candidate.length(), diff);
      }
    }
  }
  return CONTINUE_E;
}


BaseTool::propagation_state_e
FileBrowser::file_selected(event_handle_t &event)
{
  Skinner::Signal *signal = dynamic_cast<Skinner::Signal *>(event.get_rep());
  ASSERT(signal);
  const string &filename = signal->get_vars()->get_string("filename");
  filename_ = filename;
  rescan(event);
  return CONTINUE_E;

}



BaseTool::propagation_state_e
FileBrowser::rescan(event_handle_t& /*event*/)
{
  string candidate = directory_()+"/"+filename_();
  if (filename_().length() && filename_()[0] == '/') {
    candidate = filename_;
  }

  if (validFile(candidate)) {
    throw_signal("FileBrowser::rescan_redundant");
  }

  if (filename_() == cached_filename_) {
    return CONTINUE_E;
  }

  pair<string, string> dirfile = split_filename(candidate);
  if (validDir(dirfile.first)) {
    directory_ = dirfile.first;
    directory_ = canonicalize(directory_);
    filename_ = dirfile.second;
  }
  cached_filename_ = filename_;
  vector<string> filenames =
    GetFilenamesStartingWith(directory_(), "");
  vector<string> nondirfilenames;

  sort(filenames.begin(), filenames.end());
  int count = 0;

  for (size_t i = 0; i < filenames.size(); ++i) {
    if (filenames[i].length() && filenames[i][0] != '.' || filenames[i] == "..") {
      if (validDir(directory_()+"/"+filenames[i])) {
        Var<string> file(get_vars(), "FileBrowser::filename_"+to_string(count+1));
        filenames[i] = filenames[i]+"/";
        file = filenames[i];
        count++;
      } else {
        nondirfilenames.push_back(filenames[i]);
      }
    }
  }

  for (size_t i = 0; i < nondirfilenames.size(); ++i) {
    Var<string> file(get_vars(), "FileBrowser::filename_"+to_string(count+1));
    file = nondirfilenames[i];
    count++;
  }

  for (size_t i = count; i < 100; ++i) {
    Var<string> file(get_vars(), "FileBrowser::filename_"+to_string(i+1));
    file = "";
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e
FileBrowser::process_event(event_handle_t &event)
{
  return Parent::process_event(event);
}


int
FileBrowser::get_signal_id(const string &id) const
{
  if (id == "FileBrowser::rescan_redundant") return 1;
  return 0;
}


}
}
