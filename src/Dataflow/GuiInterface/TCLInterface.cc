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

#include <Core/Thread/Thread.h>
#include <Core/Thread/Mailbox.h>
#include <Core/Thread/Semaphore.h>

#include <Core/Util/StringUtil.h>
#include <Core/Util/Assert.h>
#include <Core/Util/Environment.h>
#include <Core/Util/Debug.h>

#include <Dataflow/GuiInterface/GuiContext.h>
#include <Dataflow/GuiInterface/TCLInterface.h>

#include <boost/math/special_functions/fpclassify.hpp>
#include <tcl.h>
#include <tk.h>

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>
#include <iomanip>

//! We need this one to access TCL
extern "C" 
{
  SCISHARE Tcl_Interp* the_interp;
}

namespace SCIRun 
{

struct TCLCommandData {
  GuiCallback* object;
  void* userdata;
};

class SCISHARE TCLEventMessage {
  public:
    //! The constructor creates the semaphore 
    TCLEventMessage() :
      return_code_(0),
      delivery_semaphore_(new Semaphore("TCL EventMessage Delivery",0))
    {}

    //! delete the semaphore if it was created
    virtual ~TCLEventMessage()
    {
      if (delivery_semaphore_) delete delivery_semaphore_;
    }
          
    void wait_for_message_delivery()
    {
      //! Wait until TCL task has processed this event
      delivery_semaphore_->down();
    }
    
    void mark_message_delivered()
    {
      //! Tell semaphore the processing is done
      delivery_semaphore_->up();
    }
    
    virtual bool is_pause()
    {
      return (false);
    }
    
    //! Result from TCL
    std::string& result() { return result_; }

    //! Resulting return code from TCL function call
    int&  code() { return return_code_; } 

    //! Function to call by TCL thread
    virtual bool execute();
    
  private:
    std::string         result_;
    int                 return_code_;
    Semaphore *         delivery_semaphore_;
};

//! Event for getting value of a variable

class SCISHARE GetTCLEventMessage : public TCLEventMessage {
  public:
    GetTCLEventMessage(const std::string &var,
                    GuiContext* ctx = 0) : 
      var_(var), 
      ctx_(ctx) 
    {}
    virtual ~GetTCLEventMessage() {}
  
    //! function to call by event handler
   virtual bool execute();
  
    std::string         var_;
    std::string         value_;
    GuiContext*         ctx_;
};

class SCISHARE AsyncGetStringTCLEventMessage : public TCLEventMessage {
  public:
    AsyncGetStringTCLEventMessage(const std::string &var,
                       std::string* value,
                       GuiContext* ctx = 0) : 
      var_(var), 
      value_(value),
      ctx_(ctx) 
    {}
    virtual ~AsyncGetStringTCLEventMessage() {}
  
    //! function to call by event handler
   virtual bool execute();
  
    std::string         var_;
    std::string*        value_;
    GuiContext*         ctx_;
};

class SCISHARE AsyncGetDoubleTCLEventMessage : public TCLEventMessage {
  public:
    AsyncGetDoubleTCLEventMessage(const std::string &var,
                       double* value,
                       GuiContext* ctx = 0) : 
      var_(var), 
      value_(value),
      ctx_(ctx) 
    {}
    virtual ~AsyncGetDoubleTCLEventMessage() {}
  
    //! function to call by event handler
   virtual bool execute();
  
    std::string         var_;
    double*             value_;
    GuiContext*         ctx_;
};

class SCISHARE AsyncGetIntTCLEventMessage : public TCLEventMessage {
  public:
    AsyncGetIntTCLEventMessage(const std::string &var,
                       int* value,
                       GuiContext* ctx = 0) : 
      var_(var), 
      value_(value),
      ctx_(ctx) 
    {}
    virtual ~AsyncGetIntTCLEventMessage() {}
  
    //! function to call by event handler
   virtual bool execute();
  
    std::string         var_;
    int*                value_;
    GuiContext*         ctx_;
};

//! Event for setting a variable
class SCISHARE SetStringTCLEventMessage : public TCLEventMessage {
  public:
    SetStringTCLEventMessage(const std::string &var, 
                    const std::string& val, 
                    GuiContext* ctx = 0) : 
      var_(var), 
      value_(val), 
      ctx_(ctx)
    {}
    virtual ~SetStringTCLEventMessage() {}
    
    //! function to call by event handler
    virtual bool execute();
    
    std::string         var_;
    std::string         value_;
    GuiContext*         ctx_;
};

//! Event for setting a variable
class SCISHARE SetDoubleTCLEventMessage : public TCLEventMessage {
  public:
    SetDoubleTCLEventMessage(const std::string &var, 
                    const double& val, 
                    GuiContext* ctx = 0) : 
      var_(var), 
      value_(val), 
      ctx_(ctx)
    {}
    virtual ~SetDoubleTCLEventMessage() {}
    
    //! function to call by event handler
    virtual bool execute();
    
    std::string         var_;
    double              value_;
    GuiContext*         ctx_;
};


class SCISHARE AsyncSetStringTCLEventMessage : public TCLEventMessage {
  public:
    AsyncSetStringTCLEventMessage(const std::string &var, 
                    const std::string& val, 
                    GuiContext* ctx = 0) : 
      var_(var), 
      value_(val), 
      ctx_(ctx)
    {}
    virtual ~AsyncSetStringTCLEventMessage() {}
    
    //! function to call by event handler
    virtual bool execute();
    
    std::string         var_;
    std::string         value_;
    GuiContext*         ctx_;
};

class SCISHARE AsyncSetDoubleTCLEventMessage : public TCLEventMessage {
  public:
    AsyncSetDoubleTCLEventMessage(const std::string &var, 
                    const double& val, 
                    GuiContext* ctx = 0) : 
      var_(var), 
      value_(val), 
      ctx_(ctx)
    {}
    virtual ~AsyncSetDoubleTCLEventMessage() {}
    
    //! function to call by event handler
    virtual bool execute();
    
    std::string         var_;
    double              value_;
    GuiContext*         ctx_;
};

class SCISHARE AddCommandTCLEventMessage : public TCLEventMessage {
  public:
    AddCommandTCLEventMessage(const std::string &command, 
                           GuiCallback* callback, 
                           void* userdata) : 
      command_(command), 
      callback_(callback), 
      userdata_(userdata) 
    {} 
    virtual ~AddCommandTCLEventMessage() {}

    //! function to call by event handler
    virtual bool execute();
    
    std::string         command_;
    GuiCallback*        callback_;
    void*               userdata_;
};


class SCISHARE DeleteCommandTCLEventMessage : public TCLEventMessage {
  public:
    DeleteCommandTCLEventMessage(const std::string &command) : 
      command_(command) 
    {} 
    virtual ~DeleteCommandTCLEventMessage() {}

    //! function to call by event handler
    virtual bool execute();

    std::string         command_;
  };


class SCISHARE CommandTCLEventMessage : public TCLEventMessage {
  public:
    CommandTCLEventMessage(const std::string &command,
                        GuiContext* ctx = 0) : 
      command_(command),
      ctx_(ctx) 
    {} 
    virtual ~CommandTCLEventMessage() {}

    //! function to call by event handler
    virtual bool execute();

    std::string         command_;
    GuiContext*         ctx_;
  };


class SCISHARE AsyncCommandTCLEventMessage : public TCLEventMessage {
  public:
    AsyncCommandTCLEventMessage(const std::string &command,
                        GuiContext* ctx = 0) : 
      command_(command),
      ctx_(ctx) 
    {} 
    virtual ~AsyncCommandTCLEventMessage() {}

    //! function to call by event handler
    virtual bool execute();

    std::string         command_;
    GuiContext*         ctx_;
  };


class SCISHARE PauseTCLEventMessage : public TCLEventMessage {
  public:
    PauseTCLEventMessage() {} 
    virtual ~PauseTCLEventMessage() {}


    //! function to call by event handler
    virtual bool execute();
    virtual bool is_pause() { return (true); }
    
  };

class SCISHARE SynchronizeTCLEventMessage : public TCLEventMessage {
  public:
    SynchronizeTCLEventMessage() {} 
    virtual ~SynchronizeTCLEventMessage() {}

    //! function to call by event handler
    virtual bool execute();
  };


class SCISHARE StatsTCLEventMessage : public TCLEventMessage {
  public:
    StatsTCLEventMessage() {} 
    virtual ~StatsTCLEventMessage() {}

    //! function to call by event handler
    virtual bool execute();
  };
  
///////////////////////////////////////////////////////////////////////////
// All of these are static as there can only be one TCL Thread with one
// interface to TCL.

//! The thread accessing TCL needs to have its own locking mechanism
//! Hence we need to be able to check whether we are inside the TCL thread
//! or not. Hence we store the pointer to the task having access to TCL
static Thread* tcl_thread = 0;

//! This counter counts the number of internal locks and unlocks
//! As this one is only accessed from the TCL Thread itself it does not need a
//! mutex.
static int tcl_use_cnt = 0;

//! Set up a signaling loop back to the TCL interface for signaling when 
//! TCL is free for receiving instructions
static bool signal_tcl_when_free = false;
static Tcl_ThreadId tcl_thread_id = 0;

//! Lock guarding the access to TCL
static RecursiveMutex tcl_lock("TCL Thread Lock");

static RecursiveMutex event_lock("TCL Event");
static bool event_pending = false;

//! Lock for pausing TCL/TK thread
static Mutex tcl_pause("TCL Pause lock");
static ConditionVariable tcl_wait_pause("TCL Pause Condition");

// Mail that the TCL threads checks when idle to see if any of the other
// threads posted instructions that need to be executed
static Mailbox<TCLEventMessage*> tclQueue("TCL command mailbox", 5000);

// Function call that TCL eventually calls back up on
int eventCallback(Tcl_Event* event, int flags);

static int stat_num_tcl_callbacks = 0;
static int stat_num_tcl_sync_callbacks = 0;
static int stat_num_tcl_set_callbacks = 0;
static int stat_num_tcl_get_callbacks = 0;
static int stat_num_tcl_command_callbacks = 0;
static int stat_num_tcl_pause_callbacks = 0;

// Functions passed to TCL
void addEventCallback();

int eventCallback(Tcl_Event* event, int flags)
{
  event_lock.lock();
  event_pending =  false;
  event_lock.unlock();

  //! Only try to receive a message if TCL is not in use
  if (!TCLInterface::tcl_in_use())
  {
    TCLEventMessage* em;
    while (tclQueue.tryReceive(em)) 
    {
      bool sync = em->execute();
      bool pause = em->is_pause();
      stat_num_tcl_callbacks++;

      if (sync)
      {
        em->mark_message_delivered();
        stat_num_tcl_sync_callbacks++;
      }
      else
      {
        delete em;
      }
      
      if (pause)
      {
        // atomically free tcl_pause_ mutex
        // wait for the condition variable to be
        // signaled and reobtain the mutex
        tcl_wait_pause.wait(tcl_pause);
        stat_num_tcl_pause_callbacks++;
      }
    
    }
  }

  return 1;
}

void addEventCallback()
{
  event_lock.lock();
  if ( event_pending )
  {
    event_lock.unlock();
    return;
  }
  
  event_pending = true;
  event_lock.unlock();
  
  Tcl_Event* event = reinterpret_cast<Tcl_Event*>(Tcl_Alloc(sizeof(Tcl_Event)));
  event->proc = eventCallback;
  event->nextPtr = 0;
  Tcl_ThreadQueueEvent(tcl_thread_id, event, TCL_QUEUE_HEAD);
  Tcl_ThreadAlert(tcl_thread_id);
}

int TCLInterface::OK = TCL_OK;

class TCLInterfaceHelper
{
public:
  //! Store the pointer to the TCL Thread
  static void set_tcl_thread_id(Tcl_ThreadId thread_id);

  //! Get the pointer
  static Tcl_ThreadId get_tcl_thread_id();
};

void
TCLInterface::start_processing_events()
{
  //! Forward this one to the TCLInterface as well so it can signal when internal
  //! evaluations are finished
  TCLInterfaceHelper::set_tcl_thread_id(Tcl_GetCurrentThread());
  TCLInterface::unlock();
}

void
TCLInterface::execute(const std::string& str, GuiContext* ctx)
{
  std::string result;
  eval(str, result,ctx);
}

std::string 
TCLInterface::eval(const std::string& str, GuiContext* ctx)
{
  std::string result;
  eval(str, result,ctx);
  return result;
}


void 
TCLInterface::source_once(const std::string& filename)
{
  std::string result;
  if(!eval("source \"" + filename+"\"", result)) 
  {
    char* msg = ccast_unsafe("Couldn't source file '" + filename + "'");
    Tcl_AddErrorInfo(the_interp,msg);
    Tk_BackgroundError(the_interp);
  }
}


int do_command(ClientData cd, Tcl_Interp*, int argc, CONST84 char* argv[])
{
  TCLCommandData* td= reinterpret_cast<TCLCommandData*>(cd);
  GuiArgs args(argc, argv);
  try 
  {
    td->object->tcl_command(args, td->userdata);
  } 
  catch (const std::string &message) 
  {
    args.string_ = message;
    args.have_error_ = true;
    args.have_result_ = true;
  } 
  catch (const char *message) 
  {
    args.string_ = message;
    args.have_error_ = true;
    args.have_result_ = true;
  }

  if(args.have_result_) 
  {
    Tcl_SetResult(the_interp, strdup(args.string_.c_str()),
		  (Tcl_FreeProc*)free);
  }
  return args.have_error_?TCL_ERROR:TCL_OK;
}


GuiContext* 
TCLInterface::createContext(const std::string& name)
{
  return new GuiContext(name);
}


void 
TCLInterface::post_message(const std::string& errmsg, bool err)
{
  // "Status" could also be "warning", but this function only takes a
  // bool.  In the future, perhas we should update the functions that
  // call post_message to be more expressive.
  // displayErrorWarningOrInfo() is a function in NetworkEditor.tcl.

  std::string status = "info";
  if( err ) status = "error";

  // Replace any double quotes (") with single quote (') as they break the
  // tcl parser.
  std::string fixed_errmsg = errmsg;
  for( size_t cnt = 0; cnt < fixed_errmsg.size(); cnt++ )
  {
    char ch = errmsg[cnt];
    if( ch == '"' )	ch = '\'';
    fixed_errmsg[cnt] = ch;
  }
  
  std::string command = "displayErrorWarningOrInfo \"" + fixed_errmsg + "\" " + status;
  execute( command );
}


bool
TCLInterface::get(const std::string& name, std::string& value,GuiContext* ctx)
{
  if (is_tcl_thread()) 
  {
    if (ctx && !(ctx->is_active())) return (false);

    if (string_is_list_element(name))
    {
      std::vector<int> indexes;
      std::string listname = name;
      int idx;
      while((idx = pop_last_list_index_from_string(listname)) != -1)
      {
        indexes.push_back(idx);
      }
      
      std::string list_contents;
      if (!get(listname, list_contents,ctx)) 
      {
        return false;
      }
      
      return(extract_element_from_list (list_contents, indexes, value));
    } 
    else if (string_is_map(name)) 
    {
     const char* str = Tcl_GetVar2(the_interp, 
                         ccast_unsafe(get_map_name_from_string(name)),
                         ccast_unsafe(get_map_key_from_string(name)),
                         TCL_GLOBAL_ONLY); 
     value =  std::string(str?str:"");       
      return(value != "");
    }

    // else? just do a standard gui get
    const char* str = Tcl_GetVar(the_interp, ccast_unsafe(name),TCL_GLOBAL_ONLY);
    
    value =  std::string(str?str:""); 
    if (ctx && ctx->replace_env())
    {
      replace_environment_variables(value);    
    }
      
    return (value != "");      
  }
  else 
  {
    lock();
    if (ctx && !ctx->is_active())
    {
      unlock();
      return (false);
    }   
    
    GetTCLEventMessage em(name,ctx);
    tclQueue.send(&em);
    addEventCallback();
    em.wait_for_message_delivery();
    value = em.value_;
    unlock();
    return (value != "");      
  }
}

bool
TCLInterface::get(const std::string& name, double& value,GuiContext* ctx)
{
  std::string newval;
  bool result = get(name,newval,ctx);
  from_string(newval,value);
  return (result);
}

bool
TCLInterface::get(const std::string& name, int& value,GuiContext* ctx)
{
  std::string newval;
  bool result = get(name,newval,ctx);
  from_string(newval,value);
  return (result);
}


void
TCLInterface::async_get(const std::string& name, std::string* value,GuiContext* ctx)
{
  if (is_tcl_thread()) 
  {
    if (ctx && !(ctx->is_active())) return;

    if (string_is_list_element(name))
    {
      std::vector<int> indexes;
      std::string listname = name;
      int idx;
      while((idx = pop_last_list_index_from_string(listname)) != -1)
      {
        indexes.push_back(idx);
      }
      
      std::string list_contents;
      if (!get(listname, list_contents,ctx)) return;
      
      std::string listval;
      extract_element_from_list (list_contents, indexes, listval);
      *value = listval;
      return;
    } 
    else if (string_is_map(name)) 
    {
      std::string newval;
      const char* str = Tcl_GetVar2(the_interp, 
                         ccast_unsafe(get_map_name_from_string(name)),
                         ccast_unsafe(get_map_key_from_string(name)),
                         TCL_GLOBAL_ONLY); 
      newval =  std::string(str?str:"");       
      *value = newval;
    }

    // else? just do a standard gui get
    const char* str = Tcl_GetVar(the_interp, ccast_unsafe(name),TCL_GLOBAL_ONLY);
    
    std::string newval;
    newval =  std::string(str?str:""); 
    if (ctx && ctx->replace_env())
    {
      replace_environment_variables(newval);    
    }
    
    *value = newval;
    return;      
  }
  else 
  {
    lock();
    if (ctx && !ctx->is_active())
    {
      unlock();
      return;
    }   
    
    AsyncGetStringTCLEventMessage* em =  new AsyncGetStringTCLEventMessage(name,value,ctx);
    tclQueue.send(em);
    unlock();
  }
}


void
TCLInterface::async_get(const std::string& name, double* value,GuiContext* ctx)
{
  if (is_tcl_thread()) 
  {
    if (ctx && !(ctx->is_active())) return;

    if (string_is_list_element(name))
    {
      std::vector<int> indexes;
      std::string listname = name;
      int idx;
      while((idx = pop_last_list_index_from_string(listname)) != -1)
      {
        indexes.push_back(idx);
      }
      
      std::string list_contents;
      if (!get(listname, list_contents,ctx)) return;
      
      std::string listval;
      extract_element_from_list (list_contents, indexes, listval);
      from_string(listval,*value);
      
      return;
    } 
    else if (string_is_map(name)) 
    {
      std::string newval;
      const char* str = Tcl_GetVar2(the_interp, 
                         ccast_unsafe(get_map_name_from_string(name)),
                         ccast_unsafe(get_map_key_from_string(name)),
                         TCL_GLOBAL_ONLY); 
      newval =  std::string(str?str:"");       
      from_string(newval,*value);      
    }

    // else? just do a standard gui get
    const char* str = Tcl_GetVar(the_interp, ccast_unsafe(name),TCL_GLOBAL_ONLY);
    
    std::string newval;
    newval =  std::string(str?str:""); 
    if (ctx && ctx->replace_env())
    {
      replace_environment_variables(newval);    
    }
    
    from_string(newval,*value);
    return;      
  }
  else 
  {
    lock();
    if (ctx && !ctx->is_active())
    {
      unlock();
      return;
    }   
    
    AsyncGetDoubleTCLEventMessage* em =  new AsyncGetDoubleTCLEventMessage(name,value,ctx);
    tclQueue.send(em);
    unlock();
  }
}


void
TCLInterface::async_get(const std::string& name, int* value,GuiContext* ctx)
{
  if (is_tcl_thread()) 
  {
    if (ctx && !(ctx->is_active())) return;

    if (string_is_list_element(name))
    {
      std::vector<int> indexes;
      std::string listname = name;
      int idx;
      while((idx = pop_last_list_index_from_string(listname)) != -1)
      {
        indexes.push_back(idx);
      }
      
      std::string list_contents;
      if (!get(listname, list_contents,ctx)) return;
      
      std::string listval;
      extract_element_from_list (list_contents, indexes, listval);
      from_string(listval,*value);
      
      return;
    } 
    else if (string_is_map(name)) 
    {
      std::string newval;
      const char* str = Tcl_GetVar2(the_interp, 
                         ccast_unsafe(get_map_name_from_string(name)),
                         ccast_unsafe(get_map_key_from_string(name)),
                         TCL_GLOBAL_ONLY); 
      newval =  std::string(str?str:"");       
      from_string(newval,*value);      
    }

    // else? just do a standard gui get
    const char* str = Tcl_GetVar(the_interp, ccast_unsafe(name),TCL_GLOBAL_ONLY);
    
    std::string newval;
    newval =  std::string(str?str:""); 
    if (ctx && ctx->replace_env())
    {
      replace_environment_variables(newval);    
    }
    
    from_string(newval,*value);
    return;      
  }
  else 
  {
    lock();
    if (ctx && !ctx->is_active())
    {
      unlock();
      return;
    }   
    
    AsyncGetIntTCLEventMessage* em =  new AsyncGetIntTCLEventMessage(name,value,ctx);
    tclQueue.send(em);
    unlock();
  }
}

bool
TCLInterface::set(const std::string& name, const std::string& value,GuiContext* ctx)
{
  if (is_tcl_thread()) 
  {
    if (ctx && !(ctx->is_active())) return (false);
    
    if (string_is_list_element(name))
    {
      std::vector<int> indexes;
      std::string listname = name;
      int idx;
      while((idx = pop_last_list_index_from_string(listname)) != -1)
      {
        indexes.push_back(idx);
      }
      
      std::string list_contents;
      if (!get(listname,list_contents)) 
      {
        return false;
      }
      
      if (!(set_element_in_list(list_contents, indexes, value))) return (false);
      return(set(listname, list_contents));
    } 
    else if (string_is_map(name)) 
    {
      const char* str = Tcl_GetVar2(the_interp, 
                                    ccast_unsafe(get_map_name_from_string(name)),
                                    ccast_unsafe(get_map_key_from_string(name)),
                                    TCL_GLOBAL_ONLY);
      std::string oldvalue =  std::string(str?str:""); 
    
      if (oldvalue != value || value == "")
      {
        Tcl_SetVar2(the_interp, 
                    ccast_unsafe(get_map_name_from_string(name)),
                    ccast_unsafe(get_map_key_from_string(name)),
                    ccast_unsafe(value), 
                    TCL_GLOBAL_ONLY);
      }
      return (true);
    }
    
    const char* str = Tcl_GetVar(the_interp, ccast_unsafe(name),TCL_GLOBAL_ONLY);
    std::string oldvalue =  std::string(str?str:""); 

    // Prevent write traces from getting out of control
    if (oldvalue != value || value == "")
    {
      Tcl_SetVar(the_interp, ccast_unsafe(name), ccast_unsafe(value), TCL_GLOBAL_ONLY);
    }
    return (true);
  }
  else 
  {
    lock();
    if (ctx && !ctx->is_active())
    {
      unlock();
      return (false);
    }   
    
    SetStringTCLEventMessage em(name, value ,ctx);
    tclQueue.send(&em);
    addEventCallback();
    em.wait_for_message_delivery();
    unlock();
    return (true);
  }
}

std::string 
TCLInterface::convert_double_for_display(double value, GuiContext* ctx)
{
  if (!boost::math::isinf(value))
  {
    std::ostringstream stream;
    // Print the number 17 digits wide with decimal
    stream << std::setiosflags(std::ios::showpoint) << std::setprecision(17) << value;
    // Evaluate it in TCL to pare down extra 0's at the end
    return TCLInterface::eval("expr " + stream.str(), ctx);    
  }
  else
    return (value > 0 ? "+Infinity" : "-Infinity");
}

bool
TCLInterface::set(const std::string& name, const double& value,GuiContext* ctx)
{
  if (is_tcl_thread()) 
  {
    if (ctx && !(ctx->is_active())) return (false);
    
    const std::string svalue = convert_double_for_display(value, ctx);
    
    if (string_is_list_element(name))
    {
      std::vector<int> indexes;
      std::string listname = name;
      int idx;
      while((idx = pop_last_list_index_from_string(listname)) != -1)
      {
        indexes.push_back(idx);
      }
      
      std::string list_contents;
      if (!get(listname,list_contents)) 
      {
        return false;
      }
      
      if (!(set_element_in_list(list_contents, indexes, svalue))) return (false);
      return(set(listname, list_contents));
    } 
    else if (string_is_map(name)) 
    {
      const char* str = Tcl_GetVar2(the_interp, 
                                    ccast_unsafe(get_map_name_from_string(name)),
                                    ccast_unsafe(get_map_key_from_string(name)),
                                    TCL_GLOBAL_ONLY);
      std::string oldvalue =  std::string(str?str:""); 
    
      if (oldvalue != svalue)
      {
        Tcl_SetVar2(the_interp, 
                    ccast_unsafe(get_map_name_from_string(name)),
                    ccast_unsafe(get_map_key_from_string(name)),
                    ccast_unsafe(svalue), 
                    TCL_GLOBAL_ONLY);
      }
      return (true);
    }
    
    const char* str = Tcl_GetVar(the_interp, ccast_unsafe(name),TCL_GLOBAL_ONLY);
    std::string oldvalue =  std::string(str?str:""); 

    // Prevent write traces from getting out of control
    if (oldvalue != svalue)
    {
      Tcl_SetVar(the_interp, ccast_unsafe(name), ccast_unsafe(svalue), TCL_GLOBAL_ONLY);
    }
    return (true);
  }
  else 
  {
    lock();
    if (ctx && !ctx->is_active())
    {
      unlock();
      return (false);
    }   
    
    SetDoubleTCLEventMessage em(name, value ,ctx);
    tclQueue.send(&em);
    addEventCallback();
    em.wait_for_message_delivery();
    unlock();
    return (true);
  }
}

bool
TCLInterface::set(const std::string& name, const int& value,GuiContext* ctx)
{
  return(set(name,to_string(value),ctx));
}

void
TCLInterface::async_set(const std::string& name, const std::string& value,GuiContext* ctx)
{
  if (is_tcl_thread()) 
  {
    if (ctx && !(ctx->is_active())) return;
    
    if (string_is_list_element(name))
    {
      std::vector<int> indexes;
      std::string listname = name;
      int idx;
      while((idx = pop_last_list_index_from_string(listname)) != -1)
      {
        indexes.push_back(idx);
      }
      
      std::string list_contents;
      get(listname,list_contents);
      
      if (!(set_element_in_list(list_contents, indexes, value))) return;
      set(listname, list_contents);
    } 
    else if (string_is_map(name)) 
    {
      const char* str = Tcl_GetVar2(the_interp, 
                                    ccast_unsafe(get_map_name_from_string(name)),
                                    ccast_unsafe(get_map_key_from_string(name)),
                                    TCL_GLOBAL_ONLY);
      std::string oldvalue =  std::string(str?str:""); 
    
      if (oldvalue != value || value == "")
      {
        Tcl_SetVar2(the_interp, 
                    ccast_unsafe(get_map_name_from_string(name)),
                    ccast_unsafe(get_map_key_from_string(name)),
                    ccast_unsafe(value), 
                    TCL_GLOBAL_ONLY);
      }
      return;
    }

    const char* str = Tcl_GetVar(the_interp, ccast_unsafe(name),TCL_GLOBAL_ONLY);
    std::string oldvalue =  std::string(str?str:""); 

    // Prevent write traces from getting out of control
    if (oldvalue != value || value == "")
    {        
      Tcl_SetVar(the_interp, ccast_unsafe(name), ccast_unsafe(value), TCL_GLOBAL_ONLY);
    }
    return;
  }
  else 
  {
    lock();
    if (ctx && !ctx->is_active())
    {
      unlock();
      return;
    }   
    
    AsyncSetStringTCLEventMessage* em = new AsyncSetStringTCLEventMessage(name, value ,ctx);
    tclQueue.send(em);
    unlock();
    return;
  }
}

void
TCLInterface::async_set(const std::string& name, const double& value,GuiContext* ctx)
{
  if (is_tcl_thread()) 
  {
    if (ctx && !(ctx->is_active())) return;

    std::ostringstream stream;
    // Print the number 17 digits wide with decimal
    
    stream << std::setiosflags(std::ios::showpoint) << std::setprecision(17) << value;
    // Evaluate it in TCL to pare down extra 0's at the end
    const std::string svalue = TCLInterface::eval("expr "+stream.str(),ctx); 
    
    if (string_is_list_element(name))
    {
      std::vector<int> indexes;
      std::string listname = name;
      int idx;
      while((idx = pop_last_list_index_from_string(listname)) != -1)
      {
        indexes.push_back(idx);
      }
      
      std::string list_contents;
      get(listname,list_contents);
      
      if (!(set_element_in_list(list_contents, indexes, svalue))) return;
      set(listname, list_contents);
    } 
    else if (string_is_map(name)) 
    {
      const char* str = Tcl_GetVar2(the_interp, 
                                    ccast_unsafe(get_map_name_from_string(name)),
                                    ccast_unsafe(get_map_key_from_string(name)),
                                    TCL_GLOBAL_ONLY);
      std::string oldvalue =  std::string(str?str:""); 
    
      if (oldvalue != svalue)
      {
        Tcl_SetVar2(the_interp, 
                    ccast_unsafe(get_map_name_from_string(name)),
                    ccast_unsafe(get_map_key_from_string(name)),
                    ccast_unsafe(svalue), 
                    TCL_GLOBAL_ONLY);
      }
      return;
    }

    const char* str = Tcl_GetVar(the_interp, ccast_unsafe(name),TCL_GLOBAL_ONLY);
    std::string oldvalue =  std::string(str?str:""); 

    // Prevent write traces from getting out of control
    if (oldvalue != svalue)
    {        
      Tcl_SetVar(the_interp, ccast_unsafe(name), ccast_unsafe(svalue), TCL_GLOBAL_ONLY);
    }
    return;
  }
  else 
  {
    lock();
    if (ctx && !ctx->is_active())
    {
      unlock();
      return;
    }   
    
    AsyncSetDoubleTCLEventMessage* em = new AsyncSetDoubleTCLEventMessage(name, value ,ctx);
    tclQueue.send(em);
    unlock();
    return;
  }
}

void
TCLInterface::async_set(const std::string& name, const int& value,GuiContext* ctx)
{
  async_set(name,to_string(value),ctx);
}


void 
TCLInterface::delete_command( const std::string& command )
{
  if (is_tcl_thread())
  {
    Tcl_DeleteCommand(the_interp, ccast_unsafe(command));
  }
  else
  {
    lock();
    DeleteCommandTCLEventMessage em(command);
    tclQueue.send(&em);
    addEventCallback();
    em.wait_for_message_delivery();
    unlock();
  }
}


void 
TCLInterface::add_command(const std::string&command , GuiCallback* callback,
			       void* userdata)
{
  if (is_tcl_thread())
  {
    TCLCommandData* command_data=new TCLCommandData;
    command_data->object=callback;
    command_data->userdata=userdata;
    Tcl_CreateCommand(the_interp, ccast_unsafe(command),
          &do_command, command_data, 0);
  }
  else
  {
    lock();
    AddCommandTCLEventMessage em(command,callback,userdata);
    tclQueue.send(&em);
    addEventCallback();
    em.wait_for_message_delivery();
    unlock();  
  }
}



int 
TCLInterface::eval(const std::string& command, std::string& result, GuiContext* ctx)
{
  //! if we are the TCL Thread, go ahead and execute the command, otherwise
  //! add it to the queue so the tcl thread can get it later
  if (is_tcl_thread())
  {
    //! Lock the TCL task so we are not interrupted by a callback from another
    //! thread, this will block receiving events from other threads, the event 
    //! loop will ignore any calls until done.

    if (ctx && !ctx->is_active()) return (!TCL_OK);
        
    int code = Tcl_Eval(the_interp, ccast_unsafe(command));
    if(code != TCL_OK)
    {
      Tk_BackgroundError(the_interp);
    } 
    else 
    {
      result = std::string(Tcl_GetStringResult (the_interp));
    }

    return (code == TCL_OK);
  }
  else 
  {

    //! Send message to TCL
    lock();

    if (ctx && !ctx->is_active())
    {
      unlock();
      return (!TCL_OK);
    }
    
    CommandTCLEventMessage em(command,ctx);

    tclQueue.send(&em);
    addEventCallback();
    em.wait_for_message_delivery();
    unlock();

    result = em.result();
    return (em.code() == TCL_OK);
  }
}


void
TCLInterface::async_execute(const std::string& command, GuiContext* ctx)
{
  //! if we are the TCL Thread, go ahead and execute the command, otherwise
  //! add it to the queue so the tcl thread can get it later
  if (is_tcl_thread())
  {
    //! Lock the TCL task so we are not interupted by a callback from another
    //! thread, this will block receiving events from other threads, the event 
    //! loop will ignore any calls until done.

    if (ctx && !ctx->is_active()) return;
        
    Tcl_Eval(the_interp, ccast_unsafe(command));
    
    return;
  }
  else 
  {
    //! Send message to TCL
    
    AsyncCommandTCLEventMessage* em = new AsyncCommandTCLEventMessage(command,ctx);
    tclQueue.send(em);

    return;
  }
}


bool
TCLInterface::extract_element_from_list(const std::string& contents, 
                                        const std::vector <int>& indexes,
                                        std::string& value)
{
  std::string command = "lsubindex {"+contents+"}";
  const unsigned int last = indexes.size()-1;
  for (unsigned int i = 0; i <= last; ++i)
  {
    command += " "+to_string(indexes[last-i]);
  }
  return eval(command, value);
}


bool
TCLInterface::set_element_in_list(std::string& contents, 
                                  const std::vector <int>& indexes,
                                  const std::string& value)
{
  std::string command = "lreplacesubindex {"+contents+"} "+value;
  const unsigned int last = indexes.size()-1;
  for (unsigned int i = 0; i <= last; ++i)
  {
    command += " "+to_string(indexes[last - i]);
  }
  return eval(command, contents);
}

void
TCLInterface::synchronize()
{
  if (tcl_thread != Thread::self())
  {
    // Signal TCL that we want a pause and that
    // it should give up ownership of the pause mutex
    SynchronizeTCLEventMessage sem;
    tclQueue.send(&sem);
    addEventCallback();
    sem.wait_for_message_delivery();
  }
}

static Thread* tcl_pause_thread = 0;
static int     tcl_pause_cnt = 0;

void
TCLInterface::obtain_tcl_pause()
{
  // otherwise try to obtain tcl_pause_lock, which TCL needs to run
  // if it does not have it, it is frozen
  if (tcl_thread != Thread::self())
  {
    // if we already forced TCL to pause, just increase the counter
    if (tcl_pause_thread == Thread::self())
    {
      tcl_pause_cnt++;
      return;
    }

    // Signal TCL that we want a pause and that
    // it should give up ownership of the pause mutex
    PauseTCLEventMessage pem;
    tclQueue.send(&pem);
    addEventCallback();
    pem.wait_for_message_delivery();
    
    // We know that TCL is stopped when we can obtain 
    // the mutex. TCL will be waiting until we signal
    // it to get back the lock.
    tcl_pause.lock();
    tcl_pause_cnt++;
    tcl_pause_thread = Thread::self();
  }
}

void
TCLInterface::release_tcl_pause()
{

  if (tcl_thread != Thread::self())
  {

    if (tcl_pause_thread != Thread::self())
    {
      std::cerr << "Error try to unlock TCL from thread that did not lock TCL" << std::endl;
      return;
    }

    tcl_pause_cnt--;
    if (tcl_pause_cnt == 0)
    {
      // signal TCL thread pause is over
      tcl_wait_pause.conditionSignal();
      // allow TCL thread to retrieve its tcl lock
      tcl_pause_thread = 0;
      tcl_pause.unlock();
    }
  }
}

void
TCLInterface::print_stats()
{
  StatsTCLEventMessage sem;
  tclQueue.send(&sem);
  addEventCallback();
  sem.wait_for_message_delivery();
}


//! Lock the TCL thread
void
TCLInterface::lock()
{
  //! If we are inside the TCL thread we should not lock
  if (tcl_thread == Thread::self())
  {
    tcl_use_cnt++;
  }
  else
  {
    tcl_lock.lock();
  }
}

//! Unlock the TCL thread
void
TCLInterface::unlock()
{
  if (tcl_thread == Thread::self())
  {
    tcl_use_cnt--;
    if (tcl_use_cnt < 1)
    {
      tcl_use_cnt = 0;
      if (signal_tcl_when_free)
      {
        Tcl_ThreadAlert(tcl_thread_id);
      }
    }
  }
  else
  {
    tcl_lock.unlock();
  }
}

//! Try to lock the TCL thread
bool
TCLInterface::try_lock()
{
  //! If we are inside the TCL Thread we should not lock
  if (tcl_thread == Thread::self())
  {
    //! Inside we can always lock
    tcl_use_cnt++;
    return (true);
  }
  else
  {
    if (tcl_lock.tryLock())
    {
      return (true);
    }
    return (false);
  }
}

void
TCLInterface::mark_as_tcl_thread()
{
  // Set our ID
  tcl_thread = Thread::self();
  // Lock the thread
  tcl_pause.lock();
}

bool
TCLInterface::is_tcl_thread()
{
  return(tcl_thread == Thread::self());
}


bool
TCLInterface::tcl_in_use()
{
  return(tcl_use_cnt > 0);
}

void 
TCLInterfaceHelper::set_tcl_thread_id(Tcl_ThreadId thread_id)
{ 
  tcl_thread_id = thread_id;
  signal_tcl_when_free = true;
}

Tcl_ThreadId 
TCLInterfaceHelper::get_tcl_thread_id()
{ 
  return(tcl_thread_id);
}

bool
TCLInterface::string_is_list_element(const std::string& str)
{
  if (!str.length()) return false;
  return (str[str.length()-1] == ']');
}

bool
TCLInterface::string_is_map(const std::string& str)
{
  if (!str.length()) return false;
  return (str[str.length()-1] == ')');
}

std::string
TCLInterface::get_map_key_from_string(const std::string &varname)
{
  ASSERT(string_is_map(varname));
  ASSERT(varname.size() >= 4); // smallest possible name: ie a(b) = 4 chars
  std::string::size_type open_paren = varname.find_first_of("(");
  std::string::size_type close_paren = varname.find_last_of(")");
  ASSERT(open_paren && close_paren);
  std::string key = varname.substr(open_paren+1, close_paren-open_paren-1);
  
  if (string_is_map(key) || string_is_list_element(key))
  {
    std::string new_key;
    if (!get(key, new_key)) 
    {
      return key;
    }
    key = new_key;
  }
  return key;
}


std::string
TCLInterface::get_map_name_from_string(const std::string &varname)
{
  ASSERT(string_is_map(varname));
  ASSERT(varname.size() >= 4); // smallest possible name: ie a(b) = 4 chars
  
  std::string::size_type open_paren = varname.find_first_of("(");
  return varname.substr(0, open_paren);
}

int
TCLInterface::pop_last_list_index_from_string(std::string &varname)
{
  ASSERT(varname.size() >= 4); // smallest possible name: ie i[0] = 4 chars
  std::string::size_type open_bracket = varname.find_last_of("[");
  std::string::size_type close_bracket = varname.find_last_of("]");
  ASSERT(open_bracket && close_bracket);
  
  int i = -1;
  from_string(varname.substr(open_bracket+1, close_bracket-1),i);
  if (open_bracket > 0)
  {
    varname = varname.substr(0, open_bracket);
  }
  
  return i;
}

bool
TCLEventMessage::execute()
{
  std::cerr << "ERROR IN EVENT CALL\n";
  return (true);
}


bool 
GetTCLEventMessage::execute()
{
  code() = TCLInterface::get(var_,value_,ctx_);
  stat_num_tcl_get_callbacks++;
  return (true);
}

bool 
AsyncGetStringTCLEventMessage::execute()
{
  TCLInterface::async_get(var_,value_,ctx_);
  return (false);
}

bool 
AsyncGetDoubleTCLEventMessage::execute()
{
  TCLInterface::async_get(var_,value_,ctx_);
  return (false);
}

bool 
AsyncGetIntTCLEventMessage::execute()
{
  TCLInterface::async_get(var_,value_,ctx_);
  return (false);
}


bool
SetStringTCLEventMessage::execute()
{  
  code() = TCLInterface::set(var_,value_,ctx_);
  stat_num_tcl_set_callbacks++;
  return (true);
}

bool
SetDoubleTCLEventMessage::execute()
{  
  code() = TCLInterface::set(var_,value_,ctx_);
  stat_num_tcl_set_callbacks++;
  return (true);
}

bool
AsyncSetStringTCLEventMessage::execute()
{  
  TCLInterface::async_set(var_,value_,ctx_);
  return (false);
}

bool
AsyncSetDoubleTCLEventMessage::execute()
{  
  TCLInterface::async_set(var_,value_,ctx_);
  return (false);
}

bool 
AddCommandTCLEventMessage::execute()
{
  TCLCommandData* command_data=new TCLCommandData;

  command_data->object=callback_;
  command_data->userdata=userdata_;

  if (Tcl_CreateCommand(the_interp, ccast_unsafe(command_), 
                      &do_command, command_data, 0))
    code() = TCL_OK;
  else
    code() = !TCL_OK;

  return (true);
}

bool 
DeleteCommandTCLEventMessage::execute()
{
  if (Tcl_DeleteCommand(the_interp, ccast_unsafe(command_)))
    code() = TCL_OK;
  else
    code() = !TCL_OK;

  return (true);
}


bool 
CommandTCLEventMessage::execute()
{
  code() = TCLInterface::eval(command_,result(),ctx_);
  return (true);
}

bool 
AsyncCommandTCLEventMessage::execute()
{
  TCLInterface::async_execute(command_,ctx_);
  return (false);
}

bool
PauseTCLEventMessage::execute()
{
  return (true);
}

bool
SynchronizeTCLEventMessage::execute()
{
  return (true);
}

bool
StatsTCLEventMessage::execute()
{
  std::cerr << "stat_num_tcl_callbacks = "<< stat_num_tcl_callbacks<<"\n";
  std::cerr << "stat_num_tcl_sync_callbacks = "<<stat_num_tcl_sync_callbacks<<"\n";
  std::cerr << "stat_num_tcl_pause_callbacks = "<<stat_num_tcl_pause_callbacks<<"\n";
  std::cerr << "stat_num_tcl_set_callbacks = "<<stat_num_tcl_set_callbacks<<"\n";
  std::cerr << "stat_num_tcl_get_callbacks = "<<stat_num_tcl_get_callbacks<<"\n";
  std::cerr << "stat_num_tcl_command_callbacks = "<<stat_num_tcl_command_callbacks<<"\n";
  stat_num_tcl_callbacks = 0;
  stat_num_tcl_sync_callbacks = 0;
  stat_num_tcl_pause_callbacks = 0;
  stat_num_tcl_set_callbacks = 0;
  stat_num_tcl_get_callbacks = 0;
  stat_num_tcl_command_callbacks = 0;
  return (true);
}

}
