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
 *  OpenGL.cc: Render geometry using opengl
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   December 1994
 *
 */

#ifndef SCIRun_src_Dataflow_Modules_Render_OpenGL_h
#define SCIRun_src_Dataflow_Modules_Render_OpenGL_h

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <vector>


#include <Dataflow/Modules/Render/ViewWindow.h>
#include <Dataflow/Modules/Render/ViewScene.h>
#include <Dataflow/Modules/Render/Ball.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Core/Datatypes/Image.h>
#include <Core/Thread/FutureValue.h>
#include <Core/Thread/Runnable.h>
#include <Core/Thread/Thread.h>
#include <Core/Geom/GeomObj.h>
#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/Light.h>
#include <Core/Geom/Lighting.h>
#include <Core/Geom/GeomRenderMode.h>
#include <Core/Geom/View.h>
#include <Core/Geom/GeomCone.h>
#include <Core/Geom/GeomCylinder.h>
#include <Core/Geom/GeomSphere.h>
#include <Core/Geom/GeomTri.h>
#include <Core/Geom/GeomText.h>
#include <Core/Geom/GeomLine.h>
#include <Dataflow/GuiInterface/TCLInterface.h>
#include <Core/Math/MiscMath.h>
#include <Core/Util/Timer.h>

#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT
#include <boost/shared_ptr.hpp>
#include <boost/timer.hpp>

#include <queue>

#include <ClientMessageHandler.h>
#include <SCIRunMessageComposer.h>
#include <TCPMessageCommunicator.h>
#include <Core/Codec/VP8Encoder.h>
#include <Core/Codec/PNGEncoder.h>

#endif

namespace BioMesh3d
{
	class RTPMediaCommunication;
	class MediaCommunicationBase;
	class PNGMediaEncoder;
}

namespace SCIRun {

	namespace CompressionTypes
	{
		enum CompressionType
		{
			MPEG4_E,
			JPEG_E,
			PNG_E
		};
	}

class OpenGLHelper;
class GuiArgs;

class RenderWindowMsg : public UsedWithLockingHandle<Mutex>
{
  public:

    RenderWindowMsg() :
      UsedWithLockingHandle<Mutex>("RenderWindow Lock"),
      wait_("RenderWindowMsg Semaphore"),
      done_(false),
      message_(0)
    {}

    explicit RenderWindowMsg(int message) :
      UsedWithLockingHandle<Mutex>("RenderWindow Lock"),
      wait_("RenderWindowMsg Semaphore"),
      done_(false),
      message_(message)
    {}

    RenderWindowMsg(int message, int x, int y) :
      UsedWithLockingHandle<Mutex>("RenderWindow Lock"),
      wait_("RenderWindowMsg Semaphore"),
      done_(false),
      message_(message),
      send_pick_x_(x),
      send_pick_y_(y)
    {}

    RenderWindowMsg(int message, const std::string& fname, const std::string& type, int x, int y) :
      UsedWithLockingHandle<Mutex>("RenderWindow Lock"),
      wait_("RenderWindowMsg Semaphore"),
      done_(false),
      message_(message),
      file_name_(fname),
      file_type_(type),
      image_res_x_(x),
      image_res_y_(y)
    {}

    RenderWindowMsg(int message, int datamask, FutureValue<GeometryData*>* result) :
      UsedWithLockingHandle<Mutex>("RenderWindow Lock"),
      wait_("RenderWindowMsg Semaphore"),
      done_(false),
      message_(message),
      data_mask_(datamask),
      geom_data_(result)
    {}
    
    inline int message()
      { return (message_); }
    
    inline int& send_pick_x()
      { return (send_pick_x_); }

    inline int& send_pick_y()
      { return (send_pick_y_); }
    
    inline int& ret_pick_index()
      { return (ret_pick_index_); }

    inline GeomHandle& ret_pick_obj()
      { return (ret_pick_obj_); }

    inline GeomPickHandle& ret_pick_pick()
      { return (ret_pick_pick_); }
    
    inline void signal_done()
      { 
        lock.lock();
        done_ = true;
        wait_.conditionBroadcast();
        lock.unlock();
      }
      
    inline void wait_signal()
      { 
        lock.lock();
        if (!done_) wait_.wait(lock);
        lock.unlock();
      }
      
    inline int data_mask()
      { return (data_mask_); }
      
    inline FutureValue<GeometryData*>* geom_data()
      { return (geom_data_); }
      
    inline std::string& file_name()
      { return (file_name_); }

    inline std::string& file_type()
      { return (file_type_); }

    inline int& image_res_x()
      { return (image_res_x_); }

    inline int& image_res_y()
      { return (image_res_y_); }

  private:
  
    //! A semaphore that needs to be signaled when done
    ConditionVariable wait_;
    bool              done_;
    
    //! The message to the RenderWindow
    int               message_;
  
    //! Additional information to be send for picking
    int               send_pick_x_;
    int               send_pick_y_;
  
    //! Where to return picking information
    int               ret_pick_index_;
    GeomHandle        ret_pick_obj_;
    GeomPickHandle    ret_pick_pick_;    
    
    //! For GetData
    int               data_mask_;
    FutureValue<GeometryData*>*    geom_data_;

    //! For writing images
    std::string       file_name_;
    std::string       file_type_;
    
    int               image_res_x_;
    int               image_res_y_;

};

typedef LockingHandle<RenderWindowMsg> RenderWindowMsgHandle;

struct Frustum {
  double znear;
  double zfar;
  double left;
  double right;
  double bottom;
  double top;
  double width;
  double height;
};

struct HiRes {
  double nrows;
  double ncols;
  int row;
  int col;
  int resx;
  int resy;
};

class TkOpenGLContext;

#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT

// TODO: This is horrible, we should remove this
namespace DataBufferLength
{
enum DataBufferLength
{
	// HACK -- assume no screen is > 1280x1024
	MAX_BUF_LENGTH_E = 1310720, // 1310720 = 1280x1024 
	ENCODED_BUF_LENGTH_E = 1310720
};
} // end namespace BufferLength


class OpenGL : public BioMesh3d::ClientMessageHandler, 
	public BioMesh3d::SCIRunMessageComposer
{
#else // BUILD_BIOMESH3D_REMOTE_SUPPORT
class OpenGL {
#endif
public:
  OpenGL( ViewScene *, ViewWindow *);
  ~OpenGL();
  void                  kill_helper();
  void                  start_helper();
  void                  redraw_loop();
  
  void                  get_pick(int, int, GeomHandle&, GeomPickHandle&, int&);
  void                  schedule_redraw();
  void                  redraw();
  
  void                  redraw(double tbeg, double tend,
                               int ntimesteps, double frametime);
  void                  getData(int datamask,
                                FutureValue<GeometryData*>* result);
  bool                  compute_depth(const View& view,
                                      double& near, double& far);
  bool                  compute_fog_depth(const View& view,
                                          double& near, double& far,
                                          bool visible_only);
  void                  saveImage(const std::string& fname,
                                  const std::string& type,
                                  int x, int y);

  void                  saveImage(const std::string& fname,
                                  const std::string& type = "png");

  // Adds a DO_SYNC_FRAME to the send_mailbox_.
  void                  scheduleSyncFrame();
  
  // Compute world space point under cursor (x,y).  If successful,
  // set 'p' to that value & return true.  Otherwise, return false.
  bool                  pick_scene(int x, int y, Point *p);

  // Public Member Variables
  int                   xres_;
  int                   yres_;
  bool                  doing_image_p_;
  bool                  doing_movie_p_;
  std::string           movie_frame_extension_;  // Currently "png" or "ppm".
  int                   current_movie_frame_;
  bool                  add_timestamp_to_movie_frame_name_;
  std::string           movie_name_;
  // True if we want only to dump a movie frame when a DO_SYNC_FRAME
  // message is received.
  bool                  doing_sync_frame_;
  // If doing_sync_frame_ is set, only write a frame if
  // dump_sync_frame_ is true.  This will prevent dumping frames on
  // regular redraws.
  bool                  dump_sync_frame_;

  TkOpenGLContext *     tk_gl_context_;
  TkOpenGLContext *     old_tk_gl_context_;
  std::string           myname_;

#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT
  void initialize_codecs();
  void set_communicator();

  void handle_keyboard( char key_pressed, bool ctrl, bool alt, bool shift );
  void handle_mouse( short _x, short _y, 
	  bool left_click, bool middle_click, bool right_click,
	  bool ctrl, bool alt, bool shift, unsigned long timestamp );
  void handle_clipping( unsigned short plane, bool enable, 
	  bool show_frame, bool reverse_normal, float x, float y, float z, float d );
  void handle_resize( unsigned short width, unsigned short height );
  void handle_exit_scirun();
  void handle_quality_adjust( unsigned int quality, unsigned short gop );
  void handle_latency( unsigned long timestamp );
  
#endif

private:
  void                  redraw_frame();
  void                  setFrustumToWindowPortion();
  void                  deriveFrustum();
  void                  redraw_obj(ViewScene*, ViewWindow*, GeomHandle obj);
  void                  request_obj_state(ViewScene*, ViewWindow*, GeomHandle obj);
  void                  pick_draw_obj(ViewScene*, ViewWindow*, GeomHandle obj);
  void                  dump_image(const std::string&, const std::string& type = "raw");
  void                  put_scanline(int, int, Color* scanline, int repeat=1);
  void                  real_get_pick(int, int, GeomHandle&, 
                                      GeomPickHandle&, int&);
  void                  render_and_save_image( int x, int y,
                                               const std::string& fname,
                                               const std::string& type = "ppm" );
  void                  real_getData(int datamask, 
                                     FutureValue<GeometryData*>* result);
  void                  render_rotation_axis(const View &view,
                                             bool do_stereo, int i, 
                                             const Vector &eyesep);

  void                  render_scale_bar(const View &view,
                                         bool do_stereo, int i, 
                                         const Vector &eyesep);
                                         
  void                  render_clip_frames();

#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT
  void client_communicator_disconnected();
  void communication_receive_callback();
  void send_message( BioMesh3d::MessageHandle msg );
  void close_unconnected_scirun();

  void send_high_quality( const boost::system::error_code& ec, Core::Image3Handle image );
  void populate_key_map();
#endif

  // Private Member Variables
  Runnable *            helper_;
  Thread *              helper_thread_;


  ViewScene*            viewer_;
  ViewWindow*           view_window_;
  DrawInfoOpenGL*       drawinfo_;
  Frustum               frustum_;
  HiRes                 hi_res_;
  bool                  dead_;
  bool                  do_hi_res_;
  int                   max_gl_lights_;
  int                   animate_num_frames_;
  double                animate_time_begin_;
  double                animate_time_end_;
  double                animate_framerate_;
  double                znear_;
  double                zfar_;
  double                current_time_;
  unsigned int          frame_count_;
  unsigned int          total_frame_count_; // total number of frames since start
  
  //TODO: This code looks bad
  // HACK -- support data for get_pixel_depth, assume no screen is > 1280x1024
  float                 depth_buffer_[1310720]; // 1310720 = 1280x1024 
  GLdouble              modelview_matrix_[16];
  GLdouble              projection_matrix_[16];
  GLint                 viewport_matrix_[4];
  View                  cached_view_;

  // Mouse Picking variables
  int                   send_pick_x_;
  int                   send_pick_y_;
  int                   ret_pick_index_;
  GeomHandle            ret_pick_obj_;
  GeomPickHandle        ret_pick_pick_;

  // Thread Communication Mailboxes
  Mailbox<RenderWindowMsgHandle>  mailbox_;
  std::vector<RenderWindowMsgHandle> reply_;


  bool                   allocated_frame_buffer_;
  GLuint                 render_buffer_id_;
  GLuint                 render_buffer_depth_id_;
  GLuint                 frame_buffer_id_;

  bool                    disable_inertia_;

#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT
  bool m_socket_initialized_;

  size_t biomesh_frame_count_;

  Core::VP8EncoderHandle vp8_encoder_;
  Core::PNGEncoderHandle png_encoder_;

  BioMesh3d::TCPMessageCommunicatorHandle m_client_communicator_;

  std::vector<boost::shared_ptr<boost::thread> > io_workers_;
  boost::asio::io_service::work* io_service_work_;

  bool m_mouse_event_ready_;
  
  boost::shared_ptr<boost::asio::deadline_timer> m_client_connected_timer_;
    boost::asio::io_service m_io_service_;
  
   boost::shared_ptr<boost::asio::deadline_timer> m_redraw_timer_;

  std::map<std::string, std::string> m_key_map;

  boost::mutex m_initialized_mutex_;
  bool m_initialized_;
#endif
};

} // End namespace SCIRun

#endif // of #ifndef SCIRun_src_Dataflow_Modules_Render_OpenGL_h
