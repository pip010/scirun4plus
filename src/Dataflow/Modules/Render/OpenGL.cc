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

#include <sci_gl.h>
#include <sci_glx.h>
#include <sci_values.h>
#include <tcl.h>
#include <tk.h>

#include <sci_defs/bits_defs.h>
#include <sci_defs/opengl_defs.h>

#include <string.h>
#include <teem/nrrd.h>
#include <png.h>

// For the TEEM lock
#include <Core/Datatypes/NrrdData.h>

#include <Core/Util/StringUtil.h>
#include <Core/Util/Environment.h>
#include <Core/Util/FileUtils.h>
#include <Core/Geom/GeomViewerItem.h>
#include <Core/Geom/GeomQuads.h>
#include <Core/Geom/GeomResourceManager.h>


#include <Dataflow/Modules/Render/OpenGL.h>
#include <Dataflow/GuiInterface/TCLInterface.h>
#include <Dataflow/GuiInterface/TkOpenGLContext.h>
#include <Dataflow/Network/Network.h>
#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>

#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT
#include <CommunicatorBase.h>
#include <RTPMediaCommunication.h>
#include <TCPMessageCommunicator.h>

#include <Message.h>
#endif

#ifdef _WIN32
#  include <Core/Thread/Time.h>
#  undef near
#  undef far
#  undef min
#  undef max
#  ifndef BUILD_SCIRUN_STATIC
#    define SCISHARE __declspec(dllimport)
#  else
#    define SCISHARE
#  endif
#else
#  define SCISHARE
#endif

extern "C" SCISHARE Tcl_Interp* the_interp;

namespace SCIRun {

#define DO_REDRAW     0
#define DO_PICK       1
#define DO_GETDATA    2
#define REDRAW_DONE   4
#define PICK_DONE     5
#define DO_IMAGE      6
#define IMAGE_DONE    7
#define DO_SYNC_FRAME 8
#define EXECUTE_NETWORK_CALLBACK 9


	int CAPTURE_Z_DATA_HACK = 0;

	static const int pick_buffer_size = 512;
	static const double pick_window = 10.0;


	OpenGL::OpenGL( ViewScene *viewer, ViewWindow *vw) :
	xres_(0),
		yres_(0),
		doing_image_p_(false),
		doing_movie_p_(false),
		current_movie_frame_(0),
		add_timestamp_to_movie_frame_name_(false),
		movie_name_("./movie.%04d"),
		doing_sync_frame_(false),
		dump_sync_frame_(false),
		tk_gl_context_(0),
		old_tk_gl_context_(0),
		myname_("Not Intialized"),
		// private member variables
		helper_(0),
		helper_thread_(0),
		viewer_(viewer),
		view_window_(vw),
		drawinfo_(new DrawInfoOpenGL),
		dead_(false),
		do_hi_res_(false),
		max_gl_lights_(0),
		animate_num_frames_(0),
		animate_time_begin_(0.0),
		animate_time_end_(0.0),
		animate_framerate_(1.0),
		znear_(0.0),
		zfar_(0.0),
		current_time_(0.0),
		frame_count_(1),
		total_frame_count_(1),
		cached_view_(),
		send_pick_x_(0),
		send_pick_y_(0),
		ret_pick_index_(0),
		ret_pick_obj_(0),
		ret_pick_pick_(0),
		mailbox_("OpenGL renderer send mailbox",50),
		allocated_frame_buffer_(false),
		render_buffer_id_(0),
		frame_buffer_id_(0),
		disable_inertia_(false)
#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT
		,
		m_socket_initialized_( false ),
		io_service_work_( 0 ),
		m_mouse_event_ready_( false ),
		m_initialized_( false )
#endif
	{
		if (sci_getenv("SCI_REGRESSION_TESTING")) disable_inertia_ = true;
#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT

		size_t num_io_workers = 3;


		this->io_service_work_ = new boost::asio::io_service::work( this->m_io_service_ );
	
		for ( size_t j = 0; j < num_io_workers; j++ )
		{
			this->io_workers_.push_back( boost::shared_ptr<boost::thread>( new boost::thread( 
				boost::bind( &boost::asio::io_service::run, &( this->m_io_service_ ) ) ) ) );
		}

		this->m_client_connected_timer_ = boost::shared_ptr<boost::asio::deadline_timer>(
			new boost::asio::deadline_timer( this->m_io_service_ ) );
			
		this->m_client_connected_timer_->expires_from_now( boost::posix_time::minutes( 2 ) );
		this->m_client_connected_timer_->async_wait( boost::bind( &OpenGL::close_unconnected_scirun, this ) );
		
		this->m_redraw_timer_ = boost::shared_ptr<boost::asio::deadline_timer>(
			new boost::asio::deadline_timer( this->m_io_service_ ) );

		populate_key_map();

#endif
	}


	OpenGL::~OpenGL()
	{
		kill_helper();

		RenderWindowMsgHandle rwm;
		while (mailbox_.tryReceive(rwm))
		{
			rwm->signal_done();
		}

		delete drawinfo_;
		drawinfo_ = 0;

		if (tk_gl_context_)
		{
			delete tk_gl_context_;
			tk_gl_context_ = 0;
		}

#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT
		delete this->io_service_work_;
		
		this->m_client_connected_timer_.reset();
		this->m_redraw_timer_.reset();
		
		this->m_io_service_.stop();
		
		for ( size_t j = 0; j < this->io_workers_.size(); j++ )
		{
			this->io_workers_[ j ]->join();
		}

#endif
	}

	//////////////////////////////////////////////
	// Class that runs the opengl rendering

	class OpenGLHelper : public Runnable
	{
	private:
		OpenGL* opengl_;   // main class

	public:
		OpenGLHelper(OpenGL* opengl);
		virtual ~OpenGLHelper();
		virtual void run();
	};

	OpenGLHelper::OpenGLHelper(OpenGL* opengl) :
	opengl_(opengl)
	{
	}


	OpenGLHelper::~OpenGLHelper()
	{
	}

	void
		OpenGLHelper::run()
	{
		opengl_->redraw_loop();
	}


	//////////////////////////////////////////////


	void
		OpenGL::redraw(double tbeg, double tend, int nframes, double framerate)
	{
		// TODO: This is not thread safe
		//  if (dead_) return;
		animate_time_begin_ = tbeg;
		animate_time_end_ = tend;
		animate_num_frames_ = nframes;
		animate_framerate_ = framerate;

		RenderWindowMsgHandle rwm = new RenderWindowMsg(DO_REDRAW);
		mailbox_.send(rwm);
		rwm->wait_signal();
	}

	void
		OpenGL::schedule_redraw()
	{
		RenderWindowMsgHandle rwm = new RenderWindowMsg(DO_REDRAW);
		mailbox_.send(rwm);
	}

	void
		OpenGL::redraw()
	{
		RenderWindowMsgHandle rwm = new RenderWindowMsg(DO_REDRAW);
		mailbox_.send(rwm);
		rwm->wait_signal();
	}


	void
		OpenGL::start_helper()
	{
		helper_ = new OpenGLHelper(this);
		helper_thread_ = new Thread(helper_,std::string("OpenGL: "+myname_).c_str(),
			0, Thread::NotActivated);
		helper_thread_->setStackSize(1024*1024);
		helper_thread_->activate(false);
		helper_thread_->detach();

		Thread::yield();
	}


	void
		OpenGL::kill_helper()
	{
		// kill the helper thread
		dead_ = true;
		if (helper_thread_)
		{
			RenderWindowMsgHandle rwm= new RenderWindowMsg(86);
			mailbox_.send(rwm);
			rwm->wait_signal();
			helper_thread_ = 0;
		}
	}


	void
		OpenGL::redraw_loop()
	{
		// Ensure that we have a proper context
		tk_gl_context_->make_current();

		int resx = -1;
		int resy = -1;
		std::string fname, ftype;
		TimeThrottle throttle;
		throttle.start();
		double newtime = 0.0;

		// This is the frame when inertia is started.  This is needed since
		// we need to compute the number of frames that have been rendered
		// since this initial frame.
		unsigned int inertia_frame_start = 0;

		bool do_sync_frame = false;

		for (;;)
		{
			// Debug information
			//TCLInterface::print_stats();

			// check whether we need to reply to any messages
			if (reply_.size())
			{
				// Tell that all the messages are handled
				for (size_t j=0; j<reply_.size(); j++)
				{
					reply_[j]->signal_done();
				}
				// clear out message reply list
				reply_.clear();
			}

			////////////////////////////////////////
			// Synchronize GUI varibles
			view_window_->gui_inertia_mode_.reset();

			if(view_window_->gui_inertia_mode_.get() && !(disable_inertia_))
			{
				current_time_ = throttle.time();
				if (animate_framerate_ == 0)
				{
					animate_framerate_ = 30;
				}
				double frametime = 1.0 / animate_framerate_;
				const double delta = current_time_ - newtime;
				if (delta > 1.5 * frametime)
				{
					animate_framerate_ = 1.0 / delta;
					frametime = delta;
					newtime = current_time_;
				}
				if (delta > 0.85 * frametime)
				{
					animate_framerate_ *= 0.9;
					frametime = 1.0 / animate_framerate_;
					newtime = current_time_;
				}
				else if (delta < 0.5 * frametime)
				{
					animate_framerate_ *= 1.1;
					if (animate_framerate_ > 30)
					{
						animate_framerate_ = 30;
					}
					frametime = 1.0 / animate_framerate_;
					newtime = current_time_;
				}
				newtime += frametime;
				throttle.wait_for_time(newtime);

				RenderWindowMsgHandle rwm;
				while (mailbox_.tryReceive(rwm))
				{

					if (rwm->message() == 86)
					{
						throttle.stop();
						// Send out replies to all threads waiting for this one
						if (reply_.size())
						{
							for (size_t j=0; j<reply_.size(); j++) reply_[j]->signal_done();
							reply_.clear();
						} 
						// Send a reply to thread that asked us to quit
						rwm->signal_done();
						return;
					}
					else if (rwm->message() == DO_PICK)
					{
						send_pick_x_ = rwm->send_pick_x();
						send_pick_y_ = rwm->send_pick_y();

						real_get_pick(send_pick_x_, send_pick_y_, ret_pick_obj_,
							ret_pick_pick_, ret_pick_index_);
						rwm->ret_pick_index() = ret_pick_index_;
						rwm->ret_pick_obj() = ret_pick_obj_;
						rwm->ret_pick_pick() = ret_pick_pick_;

						rwm->signal_done();
					}
					else if (rwm->message() == DO_GETDATA)
					{
						int data_mask = rwm->data_mask(); 
						real_getData(data_mask, rwm->geom_data());
						rwm->signal_done();
					}
					else if (rwm->message() == DO_IMAGE)
					{
						do_hi_res_ = true;
						fname = rwm->file_name();
						ftype = rwm->file_type();
						resx = rwm->image_res_x();
						resy = rwm->image_res_y();
						rwm->signal_done();
					}
					else if (rwm->message() == DO_SYNC_FRAME)
					{
						do_sync_frame = true;
						rwm->signal_done();
					}
					else
					{
						reply_.push_back(rwm);
					}
				}

				view_window_->gui_inertia_recalculate_.reset();
				view_window_->gui_inertia_recalculate_.request();
				view_window_->gui_inertia_loop_count_.reset();
				view_window_->gui_inertia_loop_count_.request();
				view_window_->gui_inertia_x_.reset();
				view_window_->gui_inertia_x_.request();
				view_window_->gui_inertia_y_.reset();
				view_window_->gui_inertia_y_.request();
				TCLInterface::synchronize();

				if (view_window_->gui_inertia_recalculate_.get()) 
				{
					if (view_window_->gui_inertia_recalculate_.get() == 1) 
					{
						view_window_->gui_inertia_recalculate_.set(0);
						view_window_->ball_->vDown = HVect(0.0, 0.0, 0.0, 1.0);
						view_window_->ball_->vNow = 
							HVect(view_window_->gui_inertia_x_.get()/2.0, 
							view_window_->gui_inertia_y_.get()/2.0, 0.0, 1.0);
						view_window_->ball_->dragging = 1;
						view_window_->ball_->Update();
						view_window_->ball_->qNorm = view_window_->ball_->qNow.Conj();
						const double c = 1.0/view_window_->ball_->qNow.VecMag();
						view_window_->ball_->qNorm.x *= c;
						view_window_->ball_->qNorm.y *= c;
						view_window_->ball_->qNorm.z *= c;
					}

					view_window_->loop_count_ = view_window_->gui_inertia_loop_count_.get();  
					throttle.stop();
					throttle.clear();
					throttle.start();
					current_time_ = throttle.time(); 
					newtime = throttle.time()+frametime;
					inertia_frame_start = total_frame_count_;
					view_window_->gui_view_.reset();
					view_window_->gui_view_.request();
					TCLInterface::synchronize();

					View tmpview(view_window_->gui_view_.get());
					view_window_->rot_view_ = tmpview;
					Vector y_axis = tmpview.up();
					Vector z_axis = tmpview.eyep() - tmpview.lookat();
					Vector x_axis = Cross(y_axis,z_axis);
					x_axis.normalize();
					y_axis.normalize();
					view_window_->eye_dist_ = z_axis.normalize();
					view_window_->prev_trans_.load_frame(x_axis,y_axis,z_axis);
				}


				// you want to just rotate around the current rotation
				// axis - the current quaternion is viewwindow->ball_->qNow       
				// the first 3 components of this

				// This used to use newtime with the angular velocity, but we
				// replaced angular velocity with the number of rendered frames
				// per loop.

				// This needs to compute the number of frames that have been
				// rendered from when we started inertia.  This is the 0 angle
				// in the rotation.
				double angle = ((total_frame_count_-inertia_frame_start)%view_window_->loop_count_/
					(double)(view_window_->loop_count_)*(2*M_PI));
				view_window_->ball_->SetAngle(angle);
				View tmpview(view_window_->rot_view_);
				Transform tmp_trans;
				HMatrix mNow;
				view_window_->ball_->Value(mNow);
				tmp_trans.set(&mNow[0][0]);
				Transform prv = view_window_->prev_trans_;
				prv.post_trans(tmp_trans);
				HMatrix vmat;
				prv.get(&vmat[0][0]);
				Point y_a(vmat[0][1], vmat[1][1], vmat[2][1]);
				Point z_a(vmat[0][2], vmat[1][2], vmat[2][2]);
				tmpview.up(y_a.vector());

				if (view_window_->gui_inertia_mode_.get() == 1)
				{
					tmpview.eyep((z_a*(view_window_->eye_dist_))+tmpview.lookat().vector());
					view_window_->gui_view_.set(tmpview);
				}
				else if (view_window_->gui_inertia_mode_.get() == 2)
				{
					tmpview.lookat(tmpview.eyep()-
						(z_a*(view_window_->eye_dist_)).vector());
					view_window_->gui_view_.set(tmpview);
				}

			}
			else
			{
				for (;;)
				{
					RenderWindowMsgHandle rwm;
					rwm = mailbox_.receive();

					if (rwm->message() == 86)
					{
						throttle.stop();
						// Send out replies to all threads waiting for this one
						if (reply_.size())
						{
							for (size_t j=0; j<reply_.size(); j++) reply_[j]->signal_done();
							reply_.clear();
						} 
						// Send a reply to thread that asked us to quit
						rwm->signal_done();
						return;
					}
					else if (rwm->message() == DO_PICK)
					{
						send_pick_x_ = rwm->send_pick_x();
						send_pick_y_ = rwm->send_pick_y();

						real_get_pick(send_pick_x_, send_pick_y_, ret_pick_obj_,
							ret_pick_pick_, ret_pick_index_);
						rwm->ret_pick_index() = ret_pick_index_;
						rwm->ret_pick_obj() = ret_pick_obj_;
						rwm->ret_pick_pick() = ret_pick_pick_;

						rwm->signal_done();
					}
					else if (rwm->message() == DO_GETDATA)
					{
						int data_mask = rwm->data_mask(); 
						real_getData(data_mask, rwm->geom_data());
						rwm->signal_done();
					}
					else if (rwm->message() == DO_IMAGE)
					{
						do_hi_res_ = true;
						fname = rwm->file_name();
						ftype = rwm->file_type();
						resx =  rwm->image_res_x();
						resy =  rwm->image_res_y();
						reply_.push_back(rwm);
						break;
					}
					else if (rwm->message() == DO_SYNC_FRAME)
					{
						do_sync_frame = true;
						rwm->signal_done();
					}
#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT
					else if (rwm->message() == EXECUTE_NETWORK_CALLBACK)
					{
						rwm->signal_done();
					}
#endif
					else
					{
						reply_.push_back(rwm);
						break;
					}
				}

				newtime = throttle.time();
				throttle.stop();
				throttle.clear();
				throttle.start();
			}

			if (do_hi_res_)
			{
				render_and_save_image(resx, resy, fname, ftype);

				do_hi_res_ = false;
			}

			if (do_sync_frame) 
			{
				dump_sync_frame_ = true;
				// Prevent dumping a frame on the next loop iteration.
				do_sync_frame = false;
			}

            if ( doing_image_p_ || doing_movie_p_ )
            {
                bool doing_movie_p_old = doing_movie_p_;
                bool doing_image_p_old = doing_image_p_;
                int current_movie_frame_old = current_movie_frame_;
                
                doing_image_p_ = false;
                doing_movie_p_ = false;
                redraw_frame();
                doing_image_p_ = doing_image_p_old;  
                doing_movie_p_ = doing_movie_p_old;     
                current_movie_frame_ = current_movie_frame_old;        
            }
            
			redraw_frame();
            

// TODO: This needs to go into redraw_frame()
			dump_sync_frame_ = false;

			if (reply_.size())
			{
				for (size_t j=0; j<reply_.size(); j++) reply_[j]->signal_done();
				reply_.clear();
			} 

		} // end for (;;)
	}

#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT

	// TODO: Rename Send_high_quality_image
	void OpenGL::send_high_quality( const boost::system::error_code& ec, Core::Image3Handle image )
	{
		if ( ec )
		{
			return;
		}
		
		Core::CompressedImage3Handle cimage;
		this->png_encoder_->encode_image( image, cimage );
		if ( cimage )
		{
			this->send_message( compose_image_message( cimage ) );
		}
	}

	void OpenGL::handle_keyboard( char key_pressed, bool ctrl, bool alt, bool shift )
	{
		char modifier_off = '0';
		char modifier_on = '1';
		std::string key_cmd = "";
		key_cmd += ( ctrl ) ? modifier_on : modifier_off;
		key_cmd += ( alt ) ? modifier_on : modifier_off;
		key_cmd += ( shift ) ? modifier_on : modifier_off;
		key_cmd += key_pressed;
		if (m_key_map.find(key_cmd) != m_key_map.end())
		{
			std::string key_str = m_key_map[key_cmd];
			TCLInterface::async_execute( key_str );
			TCLInterface::synchronize();
		}
	}

	void OpenGL::handle_mouse( short _x,
		short _y, bool left_click, bool middle_click, bool right_click,
		bool ctrl, bool alt, bool shift, unsigned long timestamp )
	{
		std::string tcl_command = "SCIRun_Render_ViewScene_0-ViewWindow_0-c ";

		static bool prev_mouse_button_1_down = false;
		static bool prev_mouse_button_2_down = false;
		static bool prev_mouse_button_3_down = false;
		static std::string prev_model_cmd = "";
		if ( left_click )
		{
			prev_model_cmd = "mtranslate "; 
			tcl_command += prev_model_cmd;
			tcl_command += ( prev_mouse_button_1_down ) ? "move " : "start ";
		}
		else if ( middle_click )
		{
			prev_model_cmd = "mrotate ";
			tcl_command += prev_model_cmd;
			tcl_command += ( prev_mouse_button_2_down ) ? "move " : "start ";
		}
		else if ( right_click )
		{
			prev_model_cmd = "mscale ";
			tcl_command += prev_model_cmd;
			tcl_command += ( prev_mouse_button_3_down ) ? "move " : "start ";
		}
		else
		{
			tcl_command += prev_model_cmd + "end ";
		}

		prev_mouse_button_1_down = left_click;
		prev_mouse_button_2_down = middle_click;
		prev_mouse_button_3_down = right_click;

		if ( left_click || middle_click || right_click )
		{

			tcl_command += boost::lexical_cast< std::string >( _x )  + " " + 
				boost::lexical_cast< std::string >( _y );

			TCLInterface::async_execute( tcl_command );
			TCLInterface::synchronize();
		}

	}

	void OpenGL::handle_clipping( unsigned short plane, bool enable, 
		bool show_frame, bool reverse_normal, float x, float y, float z, float d )
	{
		std::string plane_str = boost::lexical_cast< std::string >( plane + 1 );
		std::string clip_cmd = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-clip-visible-" + 
			plane_str + ( enable ? " 1;" : " 0;" );
		if ( enable )
		{
			clip_cmd += "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-clip-frame-" +
				plane_str + ( show_frame ? " 1;" : " 0;" );
			clip_cmd += "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-clip-normal-reverse-" +
				plane_str + ( reverse_normal ? " 1;" : " 0;" );
			std::string x_str = boost::lexical_cast< std::string >( x ).substr( 0, x < 0 ? 5 : 4 );
			std::string y_str = boost::lexical_cast< std::string >( y ).substr( 0, y < 0 ? 5 : 4 );
			std::string z_str = boost::lexical_cast< std::string >( z ).substr( 0, z < 0 ? 5 : 4 );
			std::string d_str = boost::lexical_cast< std::string >( d ).substr( 0, d < 0 ? 5 : 4 );
			clip_cmd += "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-clip-normal-x-" +
				plane_str + " " + x_str + ";";
			clip_cmd += "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-clip-normal-y-" +
				plane_str + " " + y_str + ";";
			clip_cmd += "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-clip-normal-z-" +
				plane_str + " " + z_str + ";";
			clip_cmd += "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-clip-normal-d-" +
				plane_str + " " + d_str + ";";
		}

		clip_cmd += "::SCIRun_Render_ViewScene_0-ViewWindow_0-c clipFrame " + 
			plane_str + ";";
		clip_cmd += "::SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw;";

		TCLInterface::async_execute( clip_cmd );
		TCLInterface::synchronize();
	}

	void OpenGL::handle_resize( unsigned short width, unsigned short height )
	{
		std::string tcl_command = "SCIRun_Render_ViewScene_0-ViewWindow_0 resizeWindow " + 
			boost::lexical_cast< std::string >( width ) + " " + boost::lexical_cast< std::string >( height ); 
		if ( width  > 0 && height > 0 )
		{
			TCLInterface::async_execute( tcl_command );
			TCLInterface::synchronize();
		}
	}

	void OpenGL::handle_exit_scirun()
	{
		exit( 0 );
	}

	void OpenGL::populate_key_map()
	{
		char modifier_off = '0';
		char modifier_on = '1';
		std::string key_cmd;
		key_cmd.resize( 3 );
		key_cmd[0] = modifier_off; //control modifier
		key_cmd[1] = modifier_off; //alt modifier
		key_cmd[2] = modifier_off; //shift modifier

		// First set all the keys without modifiers
		m_key_map[key_cmd + "0"] = "SCIRun_Render_ViewScene_0-ViewWindow_0-c autoview";
		m_key_map[key_cmd + "1"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos x0_z1 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "2"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos x1_z1 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "3"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos y0_z1 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "4"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos y1_z1 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "5"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos z0_x0 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "6"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos z1_x0 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "7"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos x0_y1 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "8"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos x0_y0 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "X"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos closest ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "L"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changelock";
		m_key_map[key_cmd + "O"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changeorientation ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";
		m_key_map[key_cmd + "A"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changeaxes ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";
		m_key_map[key_cmd + "W"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changewire ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";
		m_key_map[key_cmd + "F"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changeflat ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";
		m_key_map[key_cmd + "P"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changeortho ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";
		m_key_map[key_cmd + "K"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changelight ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";
		m_key_map[key_cmd + "D"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changefog ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";
		m_key_map[key_cmd + "B"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changebbox ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";  
		m_key_map[key_cmd + "C"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changeclip ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";
		m_key_map[key_cmd + "U"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changecull ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";
		m_key_map[key_cmd + "S"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changestereo ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";
		m_key_map[key_cmd + "M"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 changescalebar ; SCIRun_Render_ViewScene_0-ViewWindow_0-c redraw";
		m_key_map[key_cmd + "H"] = "SCIRun_Render_ViewScene_0-ViewWindow_0-c gohome";
		m_key_map[key_cmd + "I"] = "SCIRun_Render_ViewScene_0-ViewWindow_0 showhelpwindow";

		// Next set the keys with control modifiers
		key_cmd[0] = modifier_on;

		m_key_map[key_cmd + "0"] = "SCIRun_Render_ViewScene_0-ViewWindow_0-c scaled_autoview";
		m_key_map[key_cmd + "1"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos ViewWindow0 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "2"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos ViewWindow1 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "3"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos ViewWindow2 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "4"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos ViewWindow3 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "5"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos ViewWindow4 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "6"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos ViewWindow5 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "7"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos ViewWindow6 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "8"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos ViewWindow7 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "9"] = "setGlobal SCIRun_Render_ViewScene_0-ViewWindow_0-pos ViewWindow8 ; SCIRun_Render_ViewScene_0-ViewWindow_0-c Views";
		m_key_map[key_cmd + "H"] = "SCIRun_Render_ViewScene_0-ViewWindow_0-c sethome";
	}

	void OpenGL::handle_quality_adjust( unsigned int quality, unsigned short gop )
	{
		this->vp8_encoder_->set_bitrate( quality );
		this->vp8_encoder_->set_keyframe_interval( gop );
	}

	void OpenGL::handle_latency( unsigned long timestamp )
	{
		BioMesh3d::MessageHandle msg = compose_latency_message( timestamp );
		send_message( msg );
	}

	void OpenGL::communication_receive_callback()
	{
		RenderWindowMsgHandle rwm = new RenderWindowMsg(EXECUTE_NETWORK_CALLBACK);
		mailbox_.trySend(rwm);

	}

	void OpenGL::initialize_codecs()
	{
		this->vp8_encoder_ = Core::VP8EncoderHandle( new Core::VP8Encoder );
		this->png_encoder_ = Core::PNGEncoderHandle( new Core::PNGEncoder );
		
		this->vp8_encoder_->set_bitrate( 4000000 );
		this->vp8_encoder_->set_keyframe_interval( 10 );
	}

	void OpenGL::set_communicator()
	{
		unsigned short server_port = 40321; // 40321 is the default
		try
		{
			server_port = boost::lexical_cast< unsigned short >( sci_getenv( "SCIRUN_SERVER_PORT" ) );
		}
		catch ( boost::bad_lexical_cast& err ) 
		{ 
			server_port = 40321;
		}

		BioMesh3d::MediaCommunicationTypes::MediaCommunicationType com_type = 
			BioMesh3d::MediaCommunicationTypes::TCP_E;

		unsigned short begin_port = 43212;
		unsigned short end_port = 43312;
		if ( sci_getenv( "SCIRUN_PORT_RANGE" ) )
		{
			try
			{
				std::vector< std::string > ports_string;
				std::string port_range = std::string( sci_getenv( "SCIRUN_PORT_RANGE" ) );
				boost::algorithm::split( ports_string, port_range,
					boost::is_any_of( "-" ) );
				if ( ports_string.size() > 1 )
				{
					begin_port = boost::lexical_cast< unsigned short >( ports_string[0] );
					end_port = boost::lexical_cast< unsigned short >( ports_string[1] );
				}
			}
			catch ( std::exception& e )
			{
				begin_port = 43212;
				end_port = 43312;
			}
		}

		unsigned short client_port = begin_port;
		// until we find a port we can use
		if ( !BioMesh3d::CommunicatorBase::get_available_port( begin_port, end_port, 
			client_port, com_type ) )
		{
			exit( -1 ); // we can't connect to a port
		}


		this->m_client_communicator_ = BioMesh3d::TCPMessageCommunicatorHandle( 
			new BioMesh3d::TCPMessageCommunicator( "", client_port, 
			BOOST_BIND( &OpenGL::handle_message, this, _1 ) ) );
		
		this->m_client_communicator_->set_close_function( 
			boost::bind( &OpenGL::client_communicator_disconnected, this ) );
		this->m_client_communicator_->run();

		BioMesh3d::TCPMessageCommunicatorHandle server_communicator = 
			BioMesh3d::TCPMessageCommunicatorHandle( 
			new BioMesh3d::TCPMessageCommunicator( "127.0.0.1", server_port, 
			BOOST_BIND( &OpenGL::handle_message, this, _1 ) ) );
		server_communicator->run();
		
		// TODO: This will stall when a connection cannot be established
		while ( !server_communicator->is_connected() )
		{
			boost::this_thread::sleep( boost::posix_time::seconds( 1 ) );
		}
		
		// TODO: This all could fail..., but it is not handled
		BioMesh3d::MessageHandle msg = compose_notify_server_of_port_message( client_port );
		server_communicator->send_message( msg );
	}

	void OpenGL::client_communicator_disconnected()
	{
		exit( 0 );
	}

#endif

	void
		OpenGL::render_and_save_image( int x, int y,
		const std::string & fname, const std::string & ftype )

	{
		// TODO: This code is extremely messy: it needs to be cleaned up

		viewer_->get_network()->add_log("VIEW SCENE: starting to write image file "+fname);
		// First figure out which type we are writing:
		// Natively we can write png, ppm, and raw

		// Determine file extension
		std::string fileext = fname.substr(fname.find_last_of('.')+1);

		// Copy file name so we can alter it (it is const on input)
		std::string filename = fname;
		// std::string filename_convert = fname;

		// Convert file extension to lower case
		fileext = string_tolower(fileext);

		bool write_nrrd = false;
		bool write_ppm = false;
		bool write_png = false;

		if (ftype == "raw") 
		{
			// Indicate that we are writing a raw
			write_nrrd = true;
			// If extension is improper add extension
			if ((fileext != "nhdr")&&(fileext != "nrrd")) filename += std::string(".nhdr");
		}

		if (ftype == "png") 
		{
			// Indicate that we are writing a png
			write_png = true;
			// Make sure extension matches
			if (fileext != "png") filename += std::string(".png");
		}
		if (ftype == "ppm") 
		{
			write_ppm = true;
			// Make sure extension matches
			if (fileext != "ppm") filename += std::string(".ppm");
		}

		if (!write_nrrd && !write_ppm && !write_png)
		{
			if (fileext == "raw") 
			{
				filename += std::string(".nhdr");
				write_nrrd = true;
			}
			else if (fileext == "nrrd") write_nrrd = true;
			else if (fileext == "nhdr") write_nrrd = true;
			else if (fileext == "png") write_png = true;
			else if (fileext == "ppm") write_ppm = true;
			else
			{
				// TODO: Need to wite this to a proper UI window.
				// This messgae disappears on WIndows, hence not really useful
				write_png = true;
				std::cout << "Unsupported image file format.\n";
				std::cout << "Put convert in your file path to support more image formats.\n";
				std::cout << "Using PNG fileformat instead.\n";
				filename += std::string(".png");
			}
		}

		// No other thread should communicate with TCL and TCL should not be issuing
		// UI commands as it can mess up synchronization on Windows and Linux

		TCLInterface::lock();  
		TCLInterface::obtain_tcl_pause();

		// OK filename and format should be known now

		if (tk_gl_context_)
		{
      tk_gl_context_->restackWindow();
		}

		// Make sure our GL context is current
		if (tk_gl_context_ != old_tk_gl_context_)
		{
			old_tk_gl_context_ = tk_gl_context_;
			tk_gl_context_->make_current();
		}

		deriveFrustum();
		// Get Viewport dimensions
		GLint vp[4];
		glGetIntegerv(GL_VIEWPORT, vp);

		// End of critical section: allow communication with TCL again and restart
		// TCL/TK event processing

		TCLInterface::release_tcl_pause();  
		TCLInterface::unlock();

		hi_res_.resx = x;
		hi_res_.resy = y;
		// The EXACT # of screen rows and columns needed to render the image
		hi_res_.ncols = (double)hi_res_.resx/(double)vp[2];
		hi_res_.nrows = (double)hi_res_.resy/(double)vp[3];

		// The # of screen rows and columns that will be rendered to make the image
		const int nrows = (int)ceil(hi_res_.nrows);
		const int ncols = (int)ceil(hi_res_.ncols);

		// Determine the type of image used

		int channel_bytes = 1;
		int num_channels = 3;  
		if (write_nrrd) num_channels = 4;
		int pix_size = channel_bytes * num_channels;


		// Obtain image (Render full image to memory)

		unsigned char* pixels = 0;
		unsigned char* tmp_row = 0;
		try
		{
			// Allocate enough memory to save the full image
			pixels = new unsigned char[hi_res_.resx*hi_res_.resy*pix_size];
			tmp_row = new unsigned char[hi_res_.resx*pix_size];   
		}
		catch (...)
		{
			if (pixels) delete pixels;
			if (tmp_row) delete tmp_row;
			std::cerr << "Not enough memory to generate image\n";
			throw;
		}

		// WORK AROUND FOR MAC, FIX ME
		// THE FIRST DRAW IS BAD. 
		// SOMEHOW WE DO NOT SET ALL STATE VARIABLES. AFTER THIS
		// CALL EVERYTHING SHOULD BE SET PROPERLY

		doing_image_p_ = true; // forces pbuffer if available
		redraw_frame();
		doing_image_p_ = false;

		for (hi_res_.row = nrows - 1; hi_res_.row >= 0; --hi_res_.row)
		{
			int read_height = hi_res_.resy - hi_res_.row * vp[3];
			read_height = (vp[3] < read_height) ? vp[3] : read_height;

			for (hi_res_.col = 0; hi_res_.col < ncols; hi_res_.col++)
			{
				// Render the col and row in the hi_res struct.
				doing_image_p_ = true; // forces pbuffer if available
				redraw_frame();
				doing_image_p_ = false;


				// We need to block TCL/Tk from doing anything while we are reading the 
				// frame buffer.
				TCLInterface::lock();
				TCLInterface::obtain_tcl_pause();
				// Tell OpenGL where to put the data in our pixel buffer
				// Read the data from OpenGL into our memory
				glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);

				glPixelStorei(GL_PACK_ALIGNMENT,1);
				glPixelStorei(GL_PACK_SKIP_PIXELS, hi_res_.col * vp[2]);
				glPixelStorei(GL_PACK_SKIP_ROWS,0);
				glPixelStorei(GL_PACK_ROW_LENGTH, hi_res_.resx);

				int read_width = hi_res_.resx - hi_res_.col * vp[2];
				read_width = (vp[2] < read_width) ? vp[2] : read_width;

				glReadPixels(0,0,read_width, read_height,
					(num_channels == 3) ? GL_RGB : GL_RGBA,
					(channel_bytes == 1) ? GL_UNSIGNED_BYTE : GL_UNSIGNED_SHORT,
					pixels+(hi_res_.resx*(hi_res_.resy - read_height - 
					vp[3]*hi_res_.row)*pix_size));

				TCLInterface::release_tcl_pause();
				TCLInterface::unlock();
			}
			// OpenGL renders upside-down to image_file writing

			unsigned char *top_row, *bot_row;   
			int top, bot;

			for (top = read_height-1, bot = 0; bot < read_height/2; top--, bot++)
			{
				top_row = pixels+(hi_res_.resx*(hi_res_.resy - read_height - vp[3]*hi_res_.row)*pix_size) + hi_res_.resx*top*pix_size;
				bot_row = pixels+(hi_res_.resx*(hi_res_.resy - read_height - vp[3]*hi_res_.row)*pix_size) + hi_res_.resx*bot*pix_size;
				memcpy(tmp_row, top_row, hi_res_.resx*pix_size);
				memcpy(top_row, bot_row, hi_res_.resx*pix_size);
				memcpy(bot_row, tmp_row, hi_res_.resx*pix_size);
			}
		}

		// Need to figure out how we need to do locking properly here. We should freeze
		// SCIRun fully when writing an image
		TCLInterface::lock();
		TCLInterface::obtain_tcl_pause();

		// Set OpenGL back to nice PixelStore values for somebody else
		glPixelStorei(GL_PACK_SKIP_PIXELS,0);
		glPixelStorei(GL_PACK_ROW_LENGTH,0);
		if(GLEW_EXT_framebuffer_object) {
			glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
		}
		CHECK_OPENGL_ERROR();

		TCLInterface::release_tcl_pause();
		TCLInterface::unlock();

		std::string netversion = TCLInterface::eval("getNetVersion");

		// Now write image
		if (write_png)
		{
			NrrdData::lock_teem();
			Nrrd* nrrd = nrrdNew();
			size_t nrrddim[3]; nrrddim[0] = num_channels; nrrddim[1] = x; nrrddim[2] = y;
			nrrdWrap_nva(nrrd,pixels,nrrdTypeUChar,3,nrrddim);
			NrrdIoState *nio = nrrdIoStateNew();
			nio->encoding = nrrdEncodingRaw;
			nio->format = nrrdFormatPNG;
			nio->endian = airMyEndian;

			nrrdKeyValueAdd(nrrd,"netversion",netversion.c_str());

			if (nrrdSave(filename.c_str(),nrrd,nio))
			{
				// Cannot use nrrd error system as it is not thread safe
				// Hence error number may come from somewhere else
				std::cerr << "Error could not save PNG image file\n";
			}
			nrrdIoStateNix(nio);
			nrrdNix(nrrd);
			NrrdData::unlock_teem();
		}

		if (write_ppm)
		{
			NrrdData::lock_teem();
			Nrrd* nrrd = nrrdNew();
			size_t nrrddim[3]; nrrddim[0] = num_channels; nrrddim[1] = x; nrrddim[2] = y;
			nrrdWrap_nva(nrrd,pixels,nrrdTypeUChar,3,nrrddim);
			NrrdIoState *nio = nrrdIoStateNew();
			nio->encoding = nrrdEncodingRaw;
			nio->format = nrrdFormatPNM;
			nio->endian = airMyEndian;

			nrrdKeyValueAdd(nrrd,"netversion",netversion.c_str());

			if (nrrdSave(filename.c_str(),nrrd,nio))
			{
				// Cannot use nrrd error system as it is not thread safe
				// Hence error number may come from somewhere else
				std::cerr << "Error could not save PNG image file\n";
			}
			nrrdIoStateNix(nio);
			nrrdNix(nrrd);
			NrrdData::unlock_teem();
		}


		if (write_nrrd)
		{
			NrrdData::lock_teem();

			Nrrd* nrrd = nrrdNew();
			size_t nrrddim[3]; nrrddim[0] = num_channels; nrrddim[1] = x; nrrddim[2] = y;
			nrrdWrap_nva(nrrd,pixels,nrrdTypeUChar,3,nrrddim);
			NrrdIoState *nio = nrrdIoStateNew();
			nio->encoding = nrrdEncodingRaw;
			nio->format = nrrdFormatNRRD;
			nio->endian = airMyEndian;

			nrrdKeyValueAdd(nrrd,"netversion",netversion.c_str());

			if (nrrdSave(filename.c_str(),nrrd,nio))
			{
				// Cannot use nrrd error system as it is not thread safe
				// Hence error number may come from somewhere else
				std::cerr << "Error could not save PNG image file\n";
			}
			nrrdIoStateNix(nio);
			nrrdNix(nrrd);

			NrrdData::unlock_teem();
		}

		// Free space used to store image
		if (pixels) delete [] pixels;
		if (tmp_row) delete [] tmp_row;

		viewer_->get_network()->add_log("VIEW SCENE: finished writing image file "+fname);

	} // end render_and_save_image()


	void
		OpenGL::redraw_frame()
	{    
		// TODO: FIX THIS - THIS CONSTRUCTION IS GENERALLY BAD AS IT IS NOT ATOMIC
		if (!tk_gl_context_ || dead_ ) return;

		// This lock ensures that no data is modified at the ports of the data
		// As long as we have this lock the number of objects will not change
		// We need to keep this lock through out the full redraw_frame call
		viewer_->geomlock_.readLock();

		// We lock TCL to ensure that no other threads are changing the options
		TCLInterface::lock();

		// Send request for variables to be synchronized
		// These ones will be precached
		// As communication with TCL is asynchronuously we send in all the requests
		// of variables we need and the subsequent synchronize call will ensure that
		// the data was actually obtained.
		// The smaller the number of synchronizes the better performance will be.

		view_window_->gui_ambient_scale_.request();
		view_window_->gui_diffuse_scale_.request();
		view_window_->gui_specular_scale_.request();
		view_window_->gui_shininess_scale_.request();
		view_window_->gui_emission_scale_.request();
		view_window_->gui_line_width_.request();
		view_window_->gui_point_size_.request();
		view_window_->gui_polygon_offset_factor_.request();
		view_window_->gui_polygon_offset_units_.request();
		view_window_->gui_text_offset_.request();
		view_window_->gui_sbase_.request();
		view_window_->gui_sr_.request();
		view_window_->gui_ortho_view_.request();
		view_window_->gui_view_.request();
		view_window_->gui_fog_visibleonly_.request();
		view_window_->gui_fog_start_.request();
		view_window_->gui_fog_end_.request();
		view_window_->gui_fogusebg_.request();
		view_window_->gui_fogcolor_.request();
		view_window_->gui_raxes_.request();
		view_window_->gui_caxes_.request();
		view_window_->gui_bgcolor_.request();
		view_window_->gui_view_.request();
		view_window_->gui_do_stereo_.request();
		view_window_->gui_raxes_.request();
		view_window_->gui_scalebar_.request();

		// Request object state for objects that do not have a fixed set of
		// gui vars  
		view_window_->requestClip();
		view_window_->requestScaleBar();
		view_window_->requestState("global");
		view_window_->requestVisible();

		// We need to synchronize here as do_for_visible relies on the state being
		// synchronized
		TCLInterface::synchronize();

		view_window_->do_for_visible(this, &OpenGL::request_obj_state);

		TCLInterface::synchronize();

		// The next session needs to be done while having full control over the
		// GUI System

		// Obtain_tcl_pause will force TCL into a waiting state, it will be forced
		// to wait until critical adjustments have been made from this thread
		// Besides freezing the TCL thread, this function will ensure that only
		// one thread has access to the GUI system by forcing other threads to wait
		// until this thread releases the tcl_pause resource.

		// This serves two functions:
		// (1) TCL/TK will not make any UI calls in these critical sections
		// (2) No other thread will be modifying the shared OpenGL Context, as certain
		//     operations such as alocating framebuffers etc, can only be done single
		//     threaded as the OpenGL standard does not require an atomic implementation
		//     of many openGL calls.

		// The only disadvantage is that during this critical state TCL is completely
		// frozen (we could unfreeze TCL during for instance volume rendering calls)
		// Hence we cannot make gui variable requests during this critical state,
		// we have to ensure that all the variables that we request are precached with
		// the TCL values before we enter the tcl_pause critical section

		// As unfreezing and freezing of TCL is not the fastest operation code has to
		// be written with this in mind.

		// Every get of a ui variable needs to be a get_cached(), as TCL may have altered
		// the value again. We cannot query TCL for its value as we need to lock TCL
		// for proper OpenGL access.

		TCLInterface::obtain_tcl_pause();

		// Increment the total_frame_count_ now that we know we will draw
		total_frame_count_++;

		// Make sure our GL context is current
		if ((tk_gl_context_ != old_tk_gl_context_))
		{
			tk_gl_context_->make_current();
			old_tk_gl_context_ = tk_gl_context_;

			sci_glew_init();

			GLint data[1];
			glGetIntegerv(GL_MAX_LIGHTS, data);
			max_gl_lights_=data[0];
		}

		// Get the window size
		int xres = tk_gl_context_->width();
		int yres = tk_gl_context_->height();

		if (xres != xres_ || yres_ != yres)
		{
			tk_gl_context_->make_current();
		}

		xres_ = xres;
		yres_ = yres;

		// Delete any pending texture objects just after make current but
		// before we draw anything.
		GeomResourceManager::delete_pending_objects();


		if (xres_ == 1) xres_ = 100;
		if (yres_ == 1) yres_ = 100;

		const bool dump_frame =
			// Saving an image
			doing_image_p_
			// Recording a movie, but we don't care about synchronized frames
			|| (doing_movie_p_ && !doing_sync_frame_)
			// Recording a movie, but we care about synchronized frames
			|| (doing_movie_p_ && doing_sync_frame_ && dump_sync_frame_);

		// Set up a framebuffer object for image or movie making.

#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT

		if ( GLEW_EXT_framebuffer_object )
		{
			glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
			CHECK_OPENGL_ERROR();
		}
		
		glDrawBuffer(GL_BACK);
		glReadBuffer(GL_BACK);
		glClearColor(0,0,0,0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glFinish();
		tk_gl_context_->swap();
				
#else
		if (dump_frame && GLEW_EXT_framebuffer_object)
#endif
		{
			// Only for writing frames are we putting data into a frame buffer
			if (!allocated_frame_buffer_) 
			{
				glGenFramebuffersEXT(1, &frame_buffer_id_);
				CHECK_OPENGL_ERROR();
				
				glGenRenderbuffersEXT(1, &render_buffer_id_);
				CHECK_OPENGL_ERROR();
				
				glGenRenderbuffersEXT(1, &render_buffer_depth_id_);
				CHECK_OPENGL_ERROR();

				//color buffer
				glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, render_buffer_id_);
				CHECK_OPENGL_ERROR();

				glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGBA, 
					xres_, yres_);
				CHECK_OPENGL_ERROR();

				// depth buffer
				glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, render_buffer_depth_id_);
				CHECK_OPENGL_ERROR();
				
				glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT32, 
					xres_, yres_);
				CHECK_OPENGL_ERROR();

				// frame buffer
				glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frame_buffer_id_);
				CHECK_OPENGL_ERROR();

				glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, 
					GL_RENDERBUFFER_EXT, render_buffer_id_);
				CHECK_OPENGL_ERROR();

				glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
					GL_RENDERBUFFER_EXT, render_buffer_depth_id_);
				CHECK_OPENGL_ERROR();

				allocated_frame_buffer_ = true;
			}

			glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frame_buffer_id_);
			CHECK_OPENGL_ERROR();

			GLint w, h;
			glGetRenderbufferParameterivEXT(GL_RENDERBUFFER_EXT, 
				GL_RENDERBUFFER_WIDTH_EXT, &w);
			CHECK_OPENGL_ERROR();
				
			glGetRenderbufferParameterivEXT(GL_RENDERBUFFER_EXT, 
				GL_RENDERBUFFER_HEIGHT_EXT, &h);
			CHECK_OPENGL_ERROR();

			if (w != xres_ || h != yres_) 
			{
				glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, render_buffer_id_);
				CHECK_OPENGL_ERROR();

				glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGBA, xres_, yres_);
				CHECK_OPENGL_ERROR();

				glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, render_buffer_depth_id_);
				CHECK_OPENGL_ERROR();

				glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24, xres_, yres_);
				CHECK_OPENGL_ERROR();
			}
		}

		// Set the area in the framebuffer or the screen where we draw opengl content
		glViewport(0, 0, xres_, yres_);

		// As we precached this value at the top of this function this will not spawn
		// a request to TCL hence we can safely call this function.
		Color bg(view_window_->gui_bgcolor_.get_cached());
		
#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT
		if ( allocated_frame_buffer_ )
		{
#else
		if (allocated_frame_buffer_ && dump_frame)
		{
#endif
			glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
			glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);

		}
		else if(GLEW_EXT_framebuffer_object)
		{
			glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
			glDrawBuffer(GL_BACK);
			glReadBuffer(GL_BACK);
		}
		else
		{
			glDrawBuffer(GL_BACK);
			glReadBuffer(GL_BACK);
		}


		glPixelStorei(GL_PACK_ALIGNMENT,   1);
		glPixelStorei(GL_PACK_ROW_LENGTH,  0);
		glPixelStorei(GL_PACK_IMAGE_HEIGHT,0);
		glPixelStorei(GL_PACK_SKIP_ROWS,   0);
		glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
		glPixelStorei(GL_PACK_SKIP_IMAGES,  0);
		glPixelStorei(GL_UNPACK_ALIGNMENT,   1);
		glPixelStorei(GL_UNPACK_ROW_LENGTH,  0);
		glPixelStorei(GL_UNPACK_IMAGE_HEIGHT,0);
		glPixelStorei(GL_UNPACK_SKIP_ROWS,   0);
		glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
		glPixelStorei(GL_UNPACK_SKIP_IMAGES,  0); 

		// Clear the screen.
		glClearColor(bg.r(), bg.g(), bg.b(), 0);


		// Setup the view...
		// Again this variable (or rather structure of gui variables) has been
		// precached, hence we can request it, without a TCL call being triggered.

		View view(view_window_->gui_view_.get_cached());
		cached_view_ = view;

		const double aspect = double(xres_)/double(yres_);
		const double fovy = RtoD(2*atan(1.0/aspect*tan(DtoR(view.fov()/2.))));

		// Get the stereo parameter again this is pre cached.
		bool do_stereo = view_window_->gui_do_stereo_.get_cached();

		// TODO: There should be an OpenGL class that records capabilities, testing
		// this over and over is pointless as the driver will not change
		if (do_stereo)
		{
			GLboolean supported;
			glGetBooleanv(GL_STEREO, &supported);
			if (!supported)
			{
				do_stereo = false;

				// TODO : this type of error messaging is bad......
				static bool warnonce = true;
				if (warnonce)
				{
					std::cout << "Stereo display selected but not supported.\n";
					warnonce = false;
				}
				//////////////////////////////////////////////////////////
			}
		}

		// TODO: Setup the drawinfo class and fix issues with this class.
		// This class is outdated and needs more OpenGL state. It was originally
		// intended to record OpenGL state so object knew what state to alter.
		// Since its creation a lot of code works around it, causing all types of
		// problems by bad design decisions over the years.
		drawinfo_->reset();

		// Get the precached TCL state  
		drawinfo_->view_ = view;
		drawinfo_->ambient_scale_ = view_window_->gui_ambient_scale_.get_cached();
		drawinfo_->diffuse_scale_ = view_window_->gui_diffuse_scale_.get_cached();
		drawinfo_->specular_scale_ = view_window_->gui_specular_scale_.get_cached();
		drawinfo_->shininess_scale_ = view_window_->gui_shininess_scale_.get_cached();
		drawinfo_->emission_scale_ = view_window_->gui_emission_scale_.get_cached();
		drawinfo_->line_width_ = view_window_->gui_line_width_.get_cached();
		drawinfo_->point_size_ = view_window_->gui_point_size_.get_cached();
		drawinfo_->polygon_offset_factor_ = view_window_->gui_polygon_offset_factor_.get_cached();
		drawinfo_->polygon_offset_units_ = view_window_->gui_polygon_offset_units_.get_cached();
		drawinfo_->text_offset_ = view_window_->gui_text_offset_.get_cached();

		if (compute_depth(view, znear_, zfar_))  
		{    // Set up graphics state.
			glDepthFunc(GL_LEQUAL);
			glEnable(GL_DEPTH_TEST);
			glEnable(GL_NORMALIZE);
			glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
			glEnable(GL_COLOR_MATERIAL);

			view_window_->setState(drawinfo_, "global");
			drawinfo_->pickmode_=0;

			CHECK_OPENGL_ERROR();

			// Do the redraw loop for each time value.
			const double frametime = 
				animate_framerate_ == 0 ? 0 : 1.0 / animate_framerate_;

			TimeThrottle throttle;
			throttle.start();
			Vector eyesep(0.0, 0.0, 0.0);

			// Setup stereo mode by altering the location of the eye point
			// Setting up separation variables
			if (do_stereo)
			{
				// Note: both gui variables are precached at the start of the function
				// and do not make a call to TCL/TK
				const double eye_sep_dist = view_window_->gui_sbase_.get_cached() *
					(view_window_->gui_sr_.get_cached() ? 0.048 : 0.0125);

				Vector u, v;
				view.get_viewplane(aspect, 1.0, u, v);
				u.safe_normalize();

				const double zmid = (znear_+zfar_) / 2.0;
				eyesep = u * eye_sep_dist * zmid;
			}

			for (int t=0; t<animate_num_frames_; t++)
			{
				int n = 1;
				if ( do_stereo ) n = 2;

				// Loop over stereo drawings
				for (int i=0; i<n; i++)
				{
					// TODO: THis code needs to fixed for saving out stereo images
					// The Stereo mode should be able to render to the frame_buffer as well
					// Note: currently stereo requires stereo support from the opengl driver
					// this is not strictly needed for saving stereo images as one could
					// render these directly.

					if ( do_stereo )
					{
						glDrawBuffer( (i == 0) ? GL_BACK_LEFT : GL_BACK_RIGHT);
					}

					// Clearing buffers again
					glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
					glViewport(0, 0, xres_, yres_);

					// const double modeltime = t * dt + animate_time_begin_;
					// view_window_->set_current_time(modeltime);

					// render normal
					// Setup view.
					glMatrixMode(GL_PROJECTION);
					glLoadIdentity();

					// This gui variable has been cached        
					if (view_window_->gui_ortho_view_.get())
					{
						drawinfo_->view_.set_is_ortho(true);
						const double len = (view.lookat() - view.eyep()).length();
						const double yval = tan(fovy * M_PI / 360.0) * len;
						const double xval = yval * aspect;
						glOrtho(-xval, xval, -yval, yval, znear_, zfar_);
					}
					else
					{
						drawinfo_->view_.set_is_ortho(false);
						gluPerspective(fovy, aspect, znear_, zfar_);
					}

					glMatrixMode(GL_MODELVIEW);

					glLoadIdentity();
					Point eyep(view.eyep());
					Point lookat(view.lookat());

					// Setup the correct eye point:
					// for stereo the eye is moved to the left for the left eye and
					// to the right for the right eye.
					// again the gui_sr_ variable is precached
					if (do_stereo)
					{
						if (i==0)
						{
							eyep -= eyesep;
							if (!view_window_->gui_sr_.get_cached())
							{
								lookat-=eyesep;
							}
						}
						else
						{
							eyep += eyesep;
							if (!view_window_->gui_sr_.get_cached())
							{
								lookat += eyesep;
							}
						}
					}

					Vector up(view.up());
					gluLookAt(eyep.x(), eyep.y(), eyep.z(),
						lookat.x(), lookat.y(), lookat.z(),
						up.x(), up.y(), up.z());

					// TODO: The tiling in SCIRun is not really working
					// Stickies need to know about tiling
					if (do_hi_res_)
					{
						setFrustumToWindowPortion();
					}


					// Set up Lighting
					glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
					const Lighting& lighting = viewer_->lighting_;

					int idx=0;
					int ii;

					for (ii=0; ii<static_cast<int>(lighting.lights.size()); ii++)
					{
						LightHandle light = lighting.lights[ii];
						light->opengl_setup(view, drawinfo_, idx);
					}

					for (ii=0; ii<idx && ii<max_gl_lights_; ii++)
					{
						glEnable((GLenum)(GL_LIGHT0 + ii));
					}
					for (;ii<max_gl_lights_;ii++)
					{
						glDisable((GLenum)(GL_LIGHT0 + ii));
					}

					// Now set up the fog stuff.
					// The fog parameters are again precached
					double fognear, fogfar;
					compute_fog_depth(view, fognear, fogfar,
						view_window_->gui_fog_visibleonly_.get_cached());

					glFogi(GL_FOG_MODE, GL_LINEAR);
					const float fnear = fognear + (fogfar - fognear) * view_window_->gui_fog_start_.get_cached();
					glFogf(GL_FOG_START, fnear);

					const double ffar = fognear + (fogfar - fognear) / Max(view_window_->gui_fog_end_.get_cached(), 0.001);
					glFogf(GL_FOG_END, ffar);

					GLfloat bgArray[4];
					if (view_window_->gui_fogusebg_.get_cached())
					{
						bgArray[0] = bg.r();
						bgArray[1] = bg.g();
						bgArray[2] = bg.b();
					}
					else
					{
						Color fogcolor(view_window_->gui_fogcolor_.get_cached());
						bgArray[0] = fogcolor.r();
						bgArray[1] = fogcolor.g();
						bgArray[2] = fogcolor.b();
					}       
					bgArray[3] = 1.0;
					glFogfv(GL_FOG_COLOR, bgArray);

					// now make the ViewWindow setup its clipping planes...
					// Note: that setClip needs requestClip to be called first to synchronize
					// current settings with TCL
					view_window_->setClip(drawinfo_);
					view_window_->setMouse(drawinfo_);

					// UNICAM addition
					glGetDoublev (GL_MODELVIEW_MATRIX, modelview_matrix_);
					glGetDoublev (GL_PROJECTION_MATRIX, projection_matrix_);
					glGetIntegerv(GL_VIEWPORT, viewport_matrix_);

					// set up point size, line size, and polygon offset
					glPointSize(drawinfo_->point_size_);
					glLineWidth(drawinfo_->line_width_);
					glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
					if (drawinfo_->polygon_offset_factor_ ||
						drawinfo_->polygon_offset_units_)
					{
						glPolygonOffset(drawinfo_->polygon_offset_factor_,
							drawinfo_->polygon_offset_units_);
						glEnable(GL_POLYGON_OFFSET_FILL);
					}
					else
					{
						glDisable(GL_POLYGON_OFFSET_FILL);
					}

					// Draw it all.
					// current_time_ = modeltime;
					view_window_->do_for_visible(this, &OpenGL::redraw_obj);
					CHECK_OPENGL_ERROR();

					// Check clip planes again, but this time draw their outlining 
					// geometry if clipping plane frames are turned on.  Do this after 
					// drawing all of the other objects.
					if(drawinfo_->clip_planes_)
					{
						render_clip_frames();
					}

					if (view_window_->gui_raxes_.get_cached())
					{
						render_rotation_axis(view, do_stereo, i, eyesep);
					}

					if (view_window_->gui_scalebar_.get_cached())
					{
						render_scale_bar(view,do_stereo,i,eyesep);
					}
				}

				// TODO: This code cannot work properly as the buffer is of a pre selected
				// size, this code will cause problems in the near future  
				// Save z-buffer data.
				if (CAPTURE_Z_DATA_HACK)
				{
					CAPTURE_Z_DATA_HACK = 0;
					glReadPixels(0, 0, xres_, yres_, GL_DEPTH_COMPONENT, GL_FLOAT,
						depth_buffer_ );
				}
				////////////////////////////////////////////////////////////////////  

				// Wait for the right time before swapping buffers
				const double realtime = t * frametime;
				if (animate_num_frames_>1)
				{
					throttle.wait_for_time(realtime);
				}


#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT
				{
				
					if ( !this->m_socket_initialized_ && sci_getenv("SCIRUN_SERVER_PORT") )
					{
						this->initialize_codecs();
						this->set_communicator();
						this->m_socket_initialized_ = true;
					}
					
					bool skip_frame = false;
					if ( this->vp8_encoder_ && this->m_client_communicator_ && this->m_client_communicator_->is_connected() ) // We should encode using the encoder
					{
						glPixelStorei(GL_PACK_ALIGNMENT, 1);

						// YES, This code is a hack. It seems on MAC there are problems with the driver.....
						glFinish();
						unsigned char temp_buffer[400];
						glReadPixels(0,0,4,4,GL_BGRA,GL_UNSIGNED_BYTE, temp_buffer );
						glFinish();
						
						Core::BufferHandle buffer = Core::Buffer::New( xres_ * yres_ * 3 );
						glReadPixels(0, 0, xres_, yres_, GL_RGB, GL_UNSIGNED_BYTE, buffer->get_data_ptr() );
						
						Core::Image3Handle image = Core::Image3::New( xres_, yres_, buffer );
						Core::CompressedImage3Handle cimage;
						
						image->flip_vertical();
						
						this->vp8_encoder_->encode_image( image, cimage );
						send_message( this->compose_image_message( cimage ) );
						
						// TODO: Not enough threads for timer, add its own io_service object

						this->m_redraw_timer_->cancel();
						this->m_redraw_timer_->expires_from_now( boost::posix_time::seconds( 2 ) );
						this->m_redraw_timer_->async_wait( boost::bind( &OpenGL::send_high_quality, this, _1, image ) );

					}
				}
#endif

				// Show the pretty picture.
				if (!( dump_frame || doing_image_p_ )) 
				{
					tk_gl_context_->swap();
				}
			}
			throttle.stop();

		}
		else // !compute_depth()
		{
			// Just show the cleared screen
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			if (!(frame_buffer_id_ && dump_frame))
			{
				tk_gl_context_->swap();
			}
		}

		TCLInterface::release_tcl_pause();
		TCLInterface::unlock();

		// TODO: This code looks supspicous
		// I doubt whether it works correctly
		if (dump_frame && doing_movie_p_)
		{
			CHECK_OPENGL_ERROR();

			std::string::size_type pos = movie_name_.find_last_of("%0");

			if ( pos == std::string::npos ||
				movie_name_[pos+2] != 'd' ||
				movie_name_.find("%") != movie_name_.find_last_of("%") )
			{
				std::string message = "Bad File Name Format - Illegal Frame Format, no C style formating found for the frame number.";

				// THIS FUNCTION NEEDS TCL
				view_window_->setMovieMessage( message, true );
				TCLInterface::synchronize();
			}
			else
			{
				  char fname[2048];
#ifdef __WIN32__
				  _snprintf(fname, 2048, movie_name_.c_str(), current_movie_frame_);
#else
				  snprintf(fname, 2048, movie_name_.c_str(), current_movie_frame_);
#endif
				  std::string fullpath = std::string(fname);
				  std::string fileext = fullpath.substr(fullpath.find_last_of('.')+1);
				  if (fileext != movie_frame_extension_)
				  {
					fullpath = std::string(fname) + "." + movie_frame_extension_;
				  }
				  
				  std::string message = "Creating image file " + fullpath;
				  view_window_->setMovieMessage( message );

				  TCLInterface::lock();
				  TCLInterface::obtain_tcl_pause();

				  dump_image(fullpath, movie_frame_extension_);

				  current_movie_frame_++;

				  TCLInterface::release_tcl_pause();
				  TCLInterface::unlock();
				  TCLInterface::synchronize();

				  view_window_->setMovieFrame(current_movie_frame_,"global");
			}
		}

		viewer_->geomlock_.readUnlock();
	}


	void
		OpenGL::get_pick(int x, int y,
		GeomHandle& pick_obj, GeomPickHandle& pick_pick,
		int& pick_index)
	{
		RenderWindowMsgHandle rwm = new RenderWindowMsg(DO_PICK,x,y);

		mailbox_.send(rwm);
		rwm->wait_signal();

		pick_obj = rwm->ret_pick_obj();
		pick_pick = rwm->ret_pick_pick();
		pick_index = rwm->ret_pick_index();
	}


	void
		OpenGL::real_get_pick(int x, int y,
		GeomHandle& pick_obj, GeomPickHandle& pick_pick,
		int& pick_index)
	{
		// TODO: FIX THIS - THIS CONSTRUCTION IS GENERALLY BAD AS IT IS NOT ATOMIC
		if (!tk_gl_context_ || dead_ ) return;

		// This lock ensures that no data is modified at the ports of the data
		// As long as we have this lock the number of objects will not change
		// We need to keep this lock through out the full redraw_frame call
		viewer_->geomlock_.readLock();

		// We lock TCL to ensure that no other threads are changing the options
		TCLInterface::lock();

		// Send request for variables to be synchronized
		// These ones will be precached
		// As communication with TCL is asynchronuously we send in all the requests
		// of variables we need and the subsequent synchronize call will ensure that
		// the data was actually obtained.
		// The smaller the number of synchronizes the better performance will be.

		view_window_->gui_ambient_scale_.request();
		view_window_->gui_diffuse_scale_.request();
		view_window_->gui_specular_scale_.request();
		view_window_->gui_shininess_scale_.request();
		view_window_->gui_emission_scale_.request();
		view_window_->gui_line_width_.request();
		view_window_->gui_point_size_.request();
		view_window_->gui_polygon_offset_factor_.request();
		view_window_->gui_polygon_offset_units_.request();
		view_window_->gui_text_offset_.request();
		view_window_->gui_sbase_.request();
		view_window_->gui_sr_.request();
		view_window_->gui_ortho_view_.request();
		view_window_->gui_view_.request();
		view_window_->gui_fog_visibleonly_.request();
		view_window_->gui_fog_start_.request();
		view_window_->gui_fog_end_.request();
		view_window_->gui_fogusebg_.request();
		view_window_->gui_fogcolor_.request();
		view_window_->gui_raxes_.request();
		view_window_->gui_caxes_.request();
		view_window_->gui_bgcolor_.request();
		view_window_->gui_view_.request();
		view_window_->gui_do_stereo_.request();
		view_window_->gui_raxes_.request();
		view_window_->gui_scalebar_.request();

		// Request object state for objects that do not have a fixed set of
		// gui vars  
		view_window_->requestClip();
		view_window_->requestScaleBar();
		view_window_->requestState("global");
		view_window_->do_for_visible(this, &OpenGL::request_obj_state);
		TCLInterface::synchronize();

		// The next session needs to be done while having full control over the
		// GUI System

		// Obtain_tcl_pause will force TCL into a waiting state, it will be forced
		// to wait until critical adjustments have been made from this thread
		// Besides freezing the TCL thread, this function will ensure that only
		// one thread has access to the GUI system by forcing other threads to wait
		// until this thread releases the tcl_pause resource.

		// This serves two functions:
		// (1) TCL/TK will not make any UI calls in these critical sections
		// (2) No other thread will be modifying the shared OpenGL Context, as certain
		//     operations such as alocating framebuffers etc, can only be done single
		//     threaded as the OpenGL standard does not require an atomic implementation
		//     of many openGL calls.

		// The only disadvantage is that during this critical state TCL is completely
		// frozen (we could unfreeze TCL during for instance volume rendering calls)
		// Hence we cannot make gui variable requests during this critical state,
		// we have to ensure that all the variables that we request are precached with
		// the TCL values before we enter the tcl_pause critical section

		// As unfreezing and freezing of TCL is not the fastest operation code has to
		// be written with this in mind.

		TCLInterface::obtain_tcl_pause();

		pick_obj = 0;
		pick_pick = 0;
		pick_index = 0x12345678;
		// Make ourselves current

		// Make sure our GL context is current
		if ((tk_gl_context_ != old_tk_gl_context_))
		{
			old_tk_gl_context_ = tk_gl_context_;
			tk_gl_context_->make_current();
		}

		// Setup the view...
		View view(view_window_->gui_view_.get());

		// Compute znear and zfar.
		double znear;
		double zfar;
		if (compute_depth(view, znear, zfar))
		{

			GLuint pick_buffer[pick_buffer_size];
			glSelectBuffer(pick_buffer_size, pick_buffer);
			glRenderMode(GL_SELECT);
			glInitNames();
#ifdef SCI_64BITS
			glPushName(0);
			glPushName(0);
			glPushName(0);
			glPushName(0);
			glPushName(0);
#else
			glPushName(0); //for the pick
			glPushName(0); //for the object
			glPushName(0); //for the object's face index
#endif

			// Picking
			const double aspect = double(xres_)/double(yres_);
			// XXX - UNICam change-- should be '1.0/aspect' not 'aspect' below
			const double fovy = RtoD(2*atan(1.0/aspect*tan(DtoR(view.fov()/2.))));
			glViewport(0, 0, xres_, yres_);
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			GLint viewport[4];
			glGetIntegerv(GL_VIEWPORT, viewport);
			gluPickMatrix(x, viewport[3]-y, pick_window, pick_window, viewport);
			if (view_window_->gui_ortho_view_.get())
			{
				drawinfo_->view_.set_is_ortho(true);
				const double len = (view.lookat() - view.eyep()).length();
				const double yval = tan(fovy * M_PI / 360.0) * len;
				const double xval = yval * aspect;
				glOrtho(-xval, xval, -yval, yval, znear, zfar);
			}
			else
			{
				drawinfo_->view_.set_is_ortho(false);
				gluPerspective(fovy, aspect, znear, zfar);
			}

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			Point eyep(view.eyep());
			Point lookat(view.lookat());
			Vector up(view.up());
			gluLookAt(eyep.x(), eyep.y(), eyep.z(),
				lookat.x(), lookat.y(), lookat.z(),
				up.x(), up.y(), up.z());

			drawinfo_->lighting_=0;
			drawinfo_->set_drawtype(DrawInfoOpenGL::Flat);
			drawinfo_->pickmode_=1;

			// Draw it all.
			view_window_->do_for_visible(this,
				(ViewWindowVisPMF)&OpenGL::pick_draw_obj);

#ifdef SCI_64BITS
			glPopName();
			glPopName();
			glPopName();
#else
			glPopName();
			glPopName();
#endif

			glFlush();
			int hits = glRenderMode(GL_RENDER);

			GLuint min_z;
#ifdef SCI_64BITS
			unsigned long long hit_obj=0;
			//    GLuint hit_obj_index = 0x12345678;
			unsigned long long hit_pick=0;
			//    GLuint hit_pick_index = 0x12345678;  // need for object indexing
#else
			GLuint hit_obj = 0;
			//GLuint hit_obj_index = 0x12345678;  // need for object indexing
			GLuint hit_pick = 0;
			//GLuint hit_pick_index = 0x12345678;  // need for object indexing
#endif

			if (drawinfo_->show_bbox_) hits = 0;

			if (hits >= 1)
			{
				int idx = 0;
				min_z = 0;
				bool have_one = false;
				for (int h=0; h<hits; h++)
				{
					int nnames = pick_buffer[idx++];
					GLuint z=pick_buffer[idx++];
					if (nnames > 1 && (!have_one || z < min_z))
					{
						min_z = z;
						have_one = true;
						idx++; // Skip Max Z
#ifdef SCI_64BITS
						idx += nnames - 5; // Skip to the last one.
						const unsigned int ho1 = pick_buffer[idx++];
						const unsigned int ho2 = pick_buffer[idx++];
						hit_pick = ((static_cast<unsigned long long>(ho1))<<32) | ho2;
						//hit_obj_index = pick_buffer[idx++];
						const unsigned int hp1 = pick_buffer[idx++];
						const unsigned int hp2 = pick_buffer[idx++];
						hit_obj = ((static_cast<unsigned long long>(hp1))<<32)|hp2;
						//hit_pick_index = pick_buffer[idx++];
						idx++;
#else
						// hit_obj=pick_buffer[idx++];
						// hit_obj_index=pick_buffer[idx++];
						//for (int i=idx; i<idx+nnames; ++i) cerr << pick_buffer[i] << "\n";
						idx += nnames - 3; // Skip to the last one.
						hit_pick = pick_buffer[idx++];
						hit_obj = pick_buffer[idx++];
						idx++;
						//hit_pick_index=pick_buffer[idx++];
#endif
					}
					else
					{
						idx += nnames + 1;
					}
				}

				if ((hit_obj != 0)&&(hit_pick != 0))
				{
					pick_obj = (GeomObj*)hit_obj;
					pick_pick = (GeomPick*)hit_pick;
					pick_obj->getId(pick_index); //(int)hit_pick_index;
				}
			}
		}

		TCLInterface::release_tcl_pause();
		TCLInterface::unlock();

		viewer_->geomlock_.readUnlock();
	}


	// Dump an image in various formats.
	// TODO: replace with calls to Teem.
	void
		OpenGL::dump_image(const std::string& fname, const std::string& ftype)
	{

	  CHECK_OPENGL_ERROR();

	  GLint vp[4];
	  glGetIntegerv(GL_VIEWPORT, vp);

	  CHECK_OPENGL_ERROR();

	  const size_t pix_size = 3;  // for RGB
	  const size_t pix_width = vp[2];
	  const size_t pix_height = vp[3];

	  const size_t num_pixels = pix_size * pix_width * pix_height;
	  std::vector<unsigned char> pixels(num_pixels);
	  
	  glFinish();
	  
	  glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
	  glPixelStorei(GL_PACK_ALIGNMENT,   1);
	  glPixelStorei(GL_PACK_ROW_LENGTH,  0);
	  glPixelStorei(GL_PACK_IMAGE_HEIGHT,0);
	  glPixelStorei(GL_PACK_SKIP_ROWS,   0);
	  glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
	  glPixelStorei(GL_PACK_SKIP_IMAGES,  0);
	  glPixelStorei(GL_UNPACK_ALIGNMENT,   1);
	  glPixelStorei(GL_UNPACK_ROW_LENGTH,  0);
	  glPixelStorei(GL_UNPACK_IMAGE_HEIGHT,0);
	  glPixelStorei(GL_UNPACK_SKIP_ROWS,   0);
	  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
	  glPixelStorei(GL_UNPACK_SKIP_IMAGES,  0);  
	  
	  // YES, This code is a hack. It seems on MAC there are problems with the driver.....
	  glFinish();
	  unsigned char buffer[400];
	  glReadPixels(0,0,4,4,GL_BGRA,GL_UNSIGNED_BYTE, buffer );
	  glFinish();
	  
	  glReadPixels(0, 0, pix_width, pix_height, GL_RGB, GL_UNSIGNED_BYTE, &pixels[0]);
	  CHECK_OPENGL_ERROR();

	  // Flip image
	  for (size_t j=0; j< (pix_height/2); j++)
	  {
		unsigned char* a = &(pixels[0]) + j*(pix_size*pix_width);
		unsigned char* b = &(pixels[0]) + (pix_height-1-j)*(pix_size*pix_width);
		unsigned char c;
		for (size_t i=0; i<(pix_size*pix_width); i++)
		{
		  c = a[i]; a[i] = b[i]; b[i] = c;
		}
	  }

	  // First figure out which type we are writing:
	  // Natively we can write png, ppm, and raw

	  // Determine file extension
	  std::string fileext = fname.substr(fname.find_last_of('.')+1);

	  // Copy file name so we can alter it (it is const on input)
	  std::string filename = fname;
	  // std::string filename_convert = fname;
	  
	  // Convert file extension to lower case
	  fileext = string_tolower(fileext);

	  bool write_nrrd = false;
	  bool write_ppm = false;
	  bool write_png = false;
	  
	  if (ftype == "raw") 
	  {
		// Indicate that we are writing a raw
		write_nrrd = true;
		// If extension is improper add extension
		if ((fileext != "nhdr")&&(fileext != "nrrd")) filename += std::string(".nhdr");
	  }
	  
	  if (ftype == "png") 
	  {
		// Indicate that we are writing a png
		write_png = true;
		// Make sure extension matches
		if (fileext != "png") filename += std::string(".png");
	  }
	  
	  if (ftype == "ppm") 
	  {
		write_ppm = true;
		// Make sure extension matches
		if (fileext != "ppm") filename += std::string(".ppm");
	  }
	  
	  if (!write_nrrd && !write_ppm && !write_png)
	  {
		if (fileext == "raw") 
		{
		  filename += std::string(".nhdr");
		  write_nrrd = true;
		}
		else if (fileext == "nrrd") write_nrrd = true;
		else if (fileext == "nhdr") write_nrrd = true;
		else if (fileext == "png") write_png = true;
		else if (fileext == "ppm") write_ppm = true;
		else
		{
		  // TODO: Need to wite this to a proper UI window.
		  // This messgae disappears on WIndows, hence not really useful
		  write_png = true;
		  std::cout << "Unsupported image file format.\n";
		  std::cout << "Put convert in your file path to support more image formats.\n";
		  std::cout << "Using PNG fileformat instead.\n";
		  filename += std::string(".png");
		}
	  }
      
      // Now write image
      if (write_png)
      {
        Nrrd* nrrd = nrrdNew();
        size_t nrrddim[3]; nrrddim[0] = 3; nrrddim[1] = pix_width; nrrddim[2] = pix_height;
        nrrdWrap_nva(nrrd,&pixels[0],nrrdTypeUChar,3,nrrddim);
        NrrdIoState *nio = nrrdIoStateNew();
        nio->encoding = nrrdEncodingRaw;
        nio->format = nrrdFormatPNG;
        nio->endian = airMyEndian;
        
        if (nrrdSave(filename.c_str(),nrrd,nio))
        {
          // Cannot use nrrd error system as it is not thread safe
          // Hence error number may come from somewhere else
          std::cerr << "Error could not save PNG image file\n";
        }
        nrrdIoStateNix(nio);
        nrrdNix(nrrd);
      }

      if (write_ppm)
      {
        Nrrd* nrrd = nrrdNew();
        size_t nrrddim[3]; nrrddim[0] = 3; nrrddim[1] = pix_width; nrrddim[2] = pix_height;
        nrrdWrap_nva(nrrd,&pixels[0],nrrdTypeUChar,3,nrrddim);
        NrrdIoState *nio = nrrdIoStateNew();
        nio->encoding = nrrdEncodingRaw;
        nio->format = nrrdFormatPNM;
        nio->endian = airMyEndian;
        
        if (nrrdSave(filename.c_str(),nrrd,nio))
        {
          // Cannot use nrrd error system as it is not thread safe
          // Hence error number may come from somewhere else
          std::cerr << "Error could not save PNG image file\n";
        }
        nrrdIoStateNix(nio);
        nrrdNix(nrrd);
      }


      if (write_nrrd)
      {
        Nrrd* nrrd = nrrdNew();
        size_t nrrddim[3]; nrrddim[0] = 3; nrrddim[1] = pix_width; nrrddim[2] = pix_height;
        nrrdWrap_nva(nrrd,&pixels[0],nrrdTypeUChar,3,nrrddim);
        NrrdIoState *nio = nrrdIoStateNew();
        nio->encoding = nrrdEncodingRaw;
        nio->format = nrrdFormatNRRD;
        nio->endian = airMyEndian;

        if (nrrdSave(filename.c_str(),nrrd,nio))
        {
          // Cannot use nrrd error system as it is not thread safe
          // Hence error number may come from somewhere else
          std::cerr << "Error could not save PNG image file\n";
        }
        nrrdIoStateNix(nio);
        nrrdNix(nrrd);
      }
	}

	void
		OpenGL::put_scanline(int y, int width, Color* scanline, int repeat)
	{
		float* pixels = new float[width*3];
		float* p = pixels;
		int i;
		for (i=0; i<width; i++)
		{
			*p++ = scanline[i].r();
			*p++ = scanline[i].g();
			*p++ = scanline[i].b();
		}
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
		glTranslated(-1.0, -1.0, 0.0);
		glScaled(2.0 / xres_, 2.0 / yres_, 1.0);
		glDepthFunc(GL_ALWAYS);
		glDrawBuffer(GL_FRONT);
		for (i=0; i<repeat; i++)
		{
			glRasterPos2i(0, y + i);
			glDrawPixels(width, 1, GL_RGB, GL_FLOAT, pixels);
		}
		glDepthFunc(GL_LEQUAL);
		glDrawBuffer(GL_BACK);
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
		delete[] pixels;
	}



	void
		OpenGL::pick_draw_obj(ViewScene* viewer, ViewWindow*, GeomHandle obj)
	{
#ifdef SCI_64BITS
		unsigned long long o = reinterpret_cast<unsigned long long>(obj.get_rep());
		unsigned int o1 = (o>>32)&0xffffffff;
		unsigned int o2 = o&0xffffffff;
		glPopName();
		glPopName();
		glPopName();
		glPushName(o1);
		glPushName(o2);
		glPushName(0x12345678);
#else
		glPopName();
		glPushName((GLuint)(obj.get_rep()));
		glPushName(0x12345678);
#endif
		obj->draw(drawinfo_, viewer->default_material_.get_rep(), current_time_);
	}

	void
		OpenGL::redraw_obj(ViewScene* viewer, ViewWindow* viewwindow, GeomHandle obj)
	{
		GeomViewerItem *gvi  = dynamic_cast<GeomViewerItem *>(obj.get_rep());
		viewwindow->setDI(drawinfo_, gvi->getString());
		obj->draw(drawinfo_, viewer->default_material_.get_rep(), current_time_);
	}

	void
		OpenGL::request_obj_state(ViewScene* viewer, ViewWindow* viewwindow, GeomHandle obj)
	{
		GeomViewerItem *gvi  = dynamic_cast<GeomViewerItem *>(obj.get_rep());
		viewwindow->requestDI(gvi->getString());
	}

	void
		OpenGL::deriveFrustum()
	{
		double pmat[16];
		glGetDoublev(GL_PROJECTION_MATRIX, pmat);
		const double G = (pmat[10]-1)/(pmat[10]+1);
		frustum_.znear = -(pmat[14]*(G-1))/(2*G);
		frustum_.zfar = frustum_.znear*G;

		if (view_window_->gui_ortho_view_.get())
		{
			drawinfo_->view_.set_is_ortho(true);
			frustum_.left = (pmat[8]-1)/pmat[0];
			frustum_.right = (pmat[8]+1)/pmat[0];
			frustum_.bottom = (pmat[9]-1)/pmat[5];
			frustum_.top = (pmat[9]+1)/pmat[5];
			frustum_.width = frustum_.right - frustum_.left;
			frustum_.height = frustum_.top - frustum_.bottom;
		}
		else
		{
			drawinfo_->view_.set_is_ortho(false);
			frustum_.left = frustum_.znear*(pmat[8]-1)/pmat[0];
			frustum_.right = frustum_.znear*(pmat[8]+1)/pmat[0];
			frustum_.bottom = frustum_.znear*(pmat[9]-1)/pmat[5];
			frustum_.top = frustum_.znear*(pmat[9]+1)/pmat[5];
			frustum_.width = frustum_.right - frustum_.left;
			frustum_.height = frustum_.top - frustum_.bottom;
		}
	}



	void
		OpenGL::setFrustumToWindowPortion()
	{
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		if (view_window_->gui_ortho_view_.get())
		{
			drawinfo_->view_.set_is_ortho(true);
			glOrtho(frustum_.left + frustum_.width / hi_res_.ncols * hi_res_.col,
				frustum_.left + frustum_.width / hi_res_.ncols * (hi_res_.col+1),
				frustum_.bottom + frustum_.height / hi_res_.nrows * hi_res_.row,
				frustum_.bottom + frustum_.height / hi_res_.nrows *(hi_res_.row+1),
				znear_, zfar_);
		}
		else
		{
			drawinfo_->view_.set_is_ortho(false);
			glFrustum(frustum_.left + frustum_.width / hi_res_.ncols * hi_res_.col,
				frustum_.left + frustum_.width / hi_res_.ncols * (hi_res_.col+1),
				frustum_.bottom + frustum_.height / hi_res_.nrows* hi_res_.row,
				frustum_.bottom + frustum_.height /hi_res_.nrows*(hi_res_.row+1),
				frustum_.znear, frustum_.zfar);
		}
	}


	void
		OpenGL::saveImage(const std::string& fname,
		const std::string& type,
		int x, int y)
	{
		RenderWindowMsgHandle rwm = new RenderWindowMsg(DO_IMAGE,fname,type,x,y);
		mailbox_.send(rwm);
		rwm->wait_signal();
	}

	void
		OpenGL::saveImage(const std::string& fname,
		const std::string& type)
	{
		RenderWindowMsgHandle rwm = new RenderWindowMsg(DO_IMAGE,fname,type,xres_,yres_);
		mailbox_.send(rwm);
		rwm->wait_signal();
	}

	void
		OpenGL::scheduleSyncFrame()
	{
		RenderWindowMsgHandle rwm = new RenderWindowMsg(DO_SYNC_FRAME);
		mailbox_.send(rwm);
	}

	void
		OpenGL::getData(int datamask, FutureValue<GeometryData*>* result)
	{
		RenderWindowMsgHandle rwm = new RenderWindowMsg(DO_GETDATA,datamask,result);
		mailbox_.send(rwm);
	}


	void
		OpenGL::real_getData(int datamask, FutureValue<GeometryData*>* result)
	{
		GeometryData* res = new GeometryData;
		if (datamask&GEOM_VIEW)
		{
			res->view=new View(cached_view_);
			res->xres=xres_;
			res->yres=yres_;
			res->znear=znear_;
			res->zfar=zfar_;
		}
		if (datamask&(GEOM_COLORBUFFER|GEOM_DEPTHBUFFER/*CollabVis*/|GEOM_MATRICES))
		{
			TCLInterface::lock();
		}
		if (datamask&GEOM_COLORBUFFER)
		{
			ColorImage* img = res->colorbuffer = new ColorImage(xres_, yres_);
			float* data=new float[xres_*yres_*3];

			glReadPixels(0, 0, xres_, yres_, GL_RGB, GL_FLOAT, data);

			float* p = data;
			for (int y=0; y<yres_; y++)
			{
				for (int x=0; x<xres_; x++)
				{
					img->put_pixel(x, y, Color(p[0], p[1], p[2]));
					p += 3;
				}
			}
			delete[] data;
		}
		if (datamask&GEOM_DEPTHBUFFER)
		{
			unsigned int* data=new unsigned int[xres_*yres_*3];

			glReadPixels(0, 0,xres_,yres_, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, data);
		}

		if (datamask&(GEOM_COLORBUFFER|GEOM_DEPTHBUFFER/*CollabVis*/|GEOM_MATRICES))
		{
			//    CHECK_OPENGL_ERROR();
			TCLInterface::unlock();
		}

		if (datamask&(GEOM_VIEW_BOUNDS))
		{
			view_window_->get_bounds_all(res->view_bounds_);
		}
		result->send(res);
	}

	// Return world-space depth to point under pixel (x, y).
	bool
		OpenGL::pick_scene( int x, int y, Point *p )
	{
		// y = 0 is bottom of screen (not top of screen, which is what X
		// events reports)
		y = (yres_ - 1) - y;
		int index = x + (y * xres_);
		double z = depth_buffer_[index];
		if (p)
		{
			// Unproject the window point (x, y, z).
			GLdouble world_x, world_y, world_z;
			gluUnProject(x, y, z,
				modelview_matrix_, projection_matrix_, viewport_matrix_,
				&world_x, &world_y, &world_z);

			*p = Point(world_x, world_y, world_z);
		}

		// if z is close to 1, then assume no object was picked
		return (z < .999999);
	}



	bool
		OpenGL::compute_depth(const View& view, double& znear, double& zfar)
	{
		znear = DBL_MAX;
		zfar =- DBL_MAX;
		BBox bb;
		view_window_->get_bounds(bb);
		if (bb.valid())
		{
			// We have something to draw.
			Point min(bb.min());
			Point max(bb.max());



			Point eyep(view.eyep());
			Vector dir(view.lookat()-eyep);
			const double dirlen2 = dir.length2();
			if (dirlen2 < 1.0e-6 || dirlen2 != dirlen2)
			{
				return false; 
			}
			dir.safe_normalize();
			const double d = -Dot(eyep, dir);
			for (int i=0;i<8;i++)
			{
				const Point p((i&1)?max.x():min.x(),
					(i&2)?max.y():min.y(),
					(i&4)?max.z():min.z());
				const double dist = Dot(p, dir) + d;

				znear = Min(znear, dist);
				zfar = Max(zfar, dist);
			}

			znear -= 0.01*znear;
			zfar  += 0.01*znear;

			if (znear <= 0.0)
			{
				if (zfar <= 0.0)
				{
					// Everything is behind us - it doesn't matter what we do.
					znear = 1.0;
					zfar = 2.0;
				}
				else
				{
					znear = zfar * 0.001;
				}
			}
			return true;
		}
		else
		{
			return false;
		}
	}


	bool
		OpenGL::compute_fog_depth(const View &view, double &znear, double &zfar,
		bool visible_only)
	{
		znear = DBL_MAX;
		zfar = -DBL_MAX;
		BBox bb;
		if (visible_only)
		{
			view_window_->get_bounds(bb);
		}
		else
		{
			view_window_->get_bounds_all(bb);
		}
		if (bb.valid())
		{
			// We have something to draw.
			Point eyep(view.eyep());
			Vector dir(view.lookat()-eyep);
			const double dirlen2 = dir.length2();
			if (dirlen2 < 1.0e-6 || dirlen2 != dirlen2)
			{
				return false;
			}
			dir.safe_normalize();
			const double d = -Dot(eyep, dir);

			// Compute distance to center of bbox.
			const double dist = Dot(bb.center(), dir);
			// Compute bbox view radius.
			const double radius = bb.diagonal().length() * dir.length2() * 0.5;

			znear = d + dist - radius;
			zfar = d + dist + radius;

			return true;
		}
		else
		{
			return false;
		}
	}


	// i is the frame number, usually refers to left or right when do_stereo
	// is set.

	void
		OpenGL::render_scale_bar(const View &view,
		bool do_stereo, int i, const Vector &eyesep)
	{

		// Get Width/Height orginal viewport
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);

		const double xscale = 1.0;

		// Compute where we need to draw
		const int xsize = static_cast<int>(viewport[2] * xscale);
		const int ysize = static_cast<int>(viewport[3] * 0.10);
		glViewport(viewport[2] - xsize, 0, xsize, ysize);

		// Now compute projection, so that projection in original image
		// and scalebar viewport is the same

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();

		const double aspect = double(xsize)/double(ysize);
		const double fovy = RtoD(2*atan(1.0/aspect*tan(DtoR(view.fov()/2.))));

		double zdist = (view.eyep() - view.lookat()).length();

		if (view_window_->gui_ortho_view_.get())
		{
			const double len = (view.lookat() - view.eyep()).length();
			const double yval = tan(fovy * M_PI / 360.0) * len;
			const double xval = yval * aspect;
			glOrtho(-1.8*xval, 0.2*xval, -yval, yval, -0.05, 0.05);
			zdist = 0.01;
		}
		else
		{
			const double len = (view.lookat() - view.eyep()).length();
			const double yval = tan(fovy * M_PI / 360.0) * len;
			const double xval = yval * aspect;
			glFrustum(-1.8*xval, 0.2*xval, -yval, yval, znear_, zfar_);

			//    gluPerspective(fovy, aspect, znear_, zfar_);
		}

		GeomHandle bar_obj = view_window_->createScaleBar();

		// Create model view based on distance between object and eye
		// Only the scaling we need to copy from the view

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		Vector up(0.0,1.0,0.0);
		Vector lookat(0.0,0.0,0.0);
		Vector eyep(0.0,0.0,zdist);

		Vector sep(1.0,0.0,0.0);
		sep = sep*eyesep.length();

		if (do_stereo)
		{
			if (i == 0)
			{
				eyep -= sep;
				if (!view_window_->gui_sr_.get())
				{
					lookat -= sep;
				}
			}
			else
			{
				eyep += sep;
				if (!view_window_->gui_sr_.get())
				{
					lookat += sep;
				}
			}
		}

		gluLookAt(eyep.x(), eyep.y(), eyep.z(),
			lookat.x(), lookat.y(), lookat.z(),
			up.x(), up.y(), up.z());
		glScalef(1.0/xscale,1.0/xscale,1.0);

		if (do_hi_res_)
		{
			// Draw in upper right hand corner of total image, not viewport image.
			const int xsize = static_cast<int>(viewport[2] * xscale);
			const int ysize = static_cast<int>(viewport[3] * 0.15);
			glViewport(viewport[2] - xsize, 0, xsize, ysize);
		}

		view_window_->setState(drawinfo_, "global");

		// Disable fog for the orientation axis.
		const bool fog = drawinfo_->fog_;
		if (fog) { glDisable(GL_FOG); }
		drawinfo_->fog_ = false;

		// Set up Lighting
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		const Lighting& l = viewer_->lighting_;
		int idx=0;
		int ii;
		for (ii=0;ii<l.lights.size();ii++)
		{
			LightHandle light=l.lights[ii];
			light->opengl_setup(view, drawinfo_, idx);
		}
		//  for (ii=0;ii<idx && ii<max_gl_lights_;ii++)
		//    glEnable((GLenum)(GL_LIGHT0+ii));
		for (ii = 0;ii<max_gl_lights_;ii++)
			glDisable((GLenum)(GL_LIGHT0+ii));

		// Disable clipping planes for the orientation icon.
		std::vector<char> cliplist(6, 0);
		for (ii = 0; ii < 6; ii++)
		{
			if (glIsEnabled((GLenum)(GL_CLIP_PLANE0+ii)))
			{
				glDisable((GLenum)(GL_CLIP_PLANE0+ii));
				cliplist[ii] = 1;
			}
		}

		// Use depthrange to force the icon to move forward.
		// Ideally the rest of the scene should be drawn at 0.05 1.0,
		// so there was no overlap at all, but that would require
		// mucking about in the picking code.
		glDepthRange(0.0 , 0.05);

		bar_obj->draw(drawinfo_, 0, current_time_);
		glDepthRange(0.0, 1.0);

		drawinfo_->fog_ = fog;  // Restore fog state.
		if (fog) { glEnable(GL_FOG); }

		// restore matrix projection
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();

		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();

		// Restore original viewport
		glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);

		// Reenable clipping planes.
		for (ii = 0; ii < 6; ii++)
		{
			if (cliplist[ii])
			{
				glEnable((GLenum)(GL_CLIP_PLANE0+ii));
			}
		}

		for (ii=0;ii<l.lights.size();ii++)
		{
			LightHandle light=l.lights[ii];
			light->opengl_setup(view, drawinfo_, idx);
		}
		for (ii=0;ii<idx && ii<max_gl_lights_;ii++)
			glEnable((GLenum)(GL_LIGHT0+ii));

	}

	// i is the frame number, usually refers to left or right when do_stereo
	// is set.

	void
		OpenGL::render_rotation_axis(const View &view,
		bool do_stereo, int i, const Vector &eyesep)
	{
		if (!(view_window_->orientation_axis_.get_rep()))
			view_window_->orientation_axis_ = view_window_->createGenAxes();
		GeomHandle axis_obj = view_window_->orientation_axis_;

		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);

		const int xysize = Min(viewport[2], viewport[3]) / 4;
		glViewport(viewport[2] - xysize, viewport[3] - xysize, xysize, xysize);
		const double aspect = 1.0;

		// fovy 16 eyedist 10 is approximately the default axis view.
		// fovy 32 eyedist 5 gives an exagerated perspective.
		const double fovy = 32.0;
		const double eyedist = 5.0;
		const double znear = eyedist - 2.0;
		const double zfar = eyedist + 2.0;

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluPerspective(fovy, aspect, znear, zfar);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		Vector oldeye(view.eyep().asVector() - view.lookat().asVector());
		oldeye.safe_normalize();
		Point eyep((oldeye * eyedist).asPoint());
		Point lookat(0.0, 0.0, 0.0);
		if (do_stereo)
		{
			if (i == 0)
			{
				eyep -= eyesep;
				if (!view_window_->gui_sr_.get())
				{
					lookat -= eyesep;
				}
			}
			else
			{
				eyep += eyesep;
				if (!view_window_->gui_sr_.get())
				{
					lookat += eyesep;
				}
			}
		}

		Vector up(view.up());
		gluLookAt(eyep.x(), eyep.y(), eyep.z(),
			lookat.x(), lookat.y(), lookat.z(),
			up.x(), up.y(), up.z());

		if (do_hi_res_)
		{
			// Draw in upper right hand corner of total image, not viewport image.
			const int xysize = Min(hi_res_.resx, hi_res_.resy) / 4;
			const int xoff = hi_res_.resx - hi_res_.col * viewport[2];
			const int yoff = hi_res_.resy - hi_res_.row * viewport[3];
			glViewport(xoff - xysize, yoff - xysize, xysize, xysize);
		}

		view_window_->setState(drawinfo_, "global");

		// Disable fog for the orientation axis.
		const bool fog = drawinfo_->fog_;
		if (fog) 
		{ 
			glDisable(GL_FOG); 
		}
		drawinfo_->fog_ = false;

		// Set up Lighting
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		const Lighting& l = viewer_->lighting_;
		int idx=0;
		int ii;
		for (ii=0;ii<l.lights.size();ii++)
		{
			LightHandle light=l.lights[ii];
			light->opengl_setup(view, drawinfo_, idx);
		}
		for (ii=0;ii<idx && ii<max_gl_lights_;ii++)
			glEnable((GLenum)(GL_LIGHT0+ii));
		for (;ii<max_gl_lights_;ii++)
			glDisable((GLenum)(GL_LIGHT0+ii));

		// Disable clipping planes for the orientation icon.
		std::vector<bool> cliplist(6, false);
		for (ii = 0; ii < 6; ii++)
		{
			if (glIsEnabled((GLenum)(GL_CLIP_PLANE0+ii)))
			{
				glDisable((GLenum)(GL_CLIP_PLANE0+ii));
				cliplist[ii] = true;
			}
		}

		// Use depthrange to force the icon to move forward.
		// Ideally the rest of the scene should be drawn at 0.05 1.0,
		// so there was no overlap at all, but that would require
		// mucking about in the picking code.
		glDepthRange(0.0, 0.05);
		axis_obj->draw(drawinfo_, 0, current_time_);
		glDepthRange(0.0, 1.0);

		drawinfo_->fog_ = fog;  // Restore fog state.
		if (fog) { glEnable(GL_FOG); }

		glMatrixMode(GL_PROJECTION);
		glPopMatrix();

		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();

		glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);

		// Reenable clipping planes.
		for (ii = 0; ii < 6; ii++)
		{
			if (cliplist[ii])
			{
				glEnable((GLenum)(GL_CLIP_PLANE0+ii));
			}
		}
	}

	void
		OpenGL::render_clip_frames()
	{
		// TODO: This code needs some updates: it currently generates the clip frame

		unsigned int i,ii;

		// first check to see if their are any frames on..
		bool frame_on = false;
		for(i = 0; i < 6; i++)
		{
			if(view_window_->viewwindow_clip_frames_draw_[i])
			{
				frame_on = true;
				break;
			}
		}

		if(!frame_on) return;

		// a frame is on, let's draw it.
		view_window_->setState(drawinfo_, "global");

		GeomGroup *frames = new GeomGroup;
		Material  *frame_color = new Material(Color(0.0, 0.0, 0.0),
			Color(.44, .50, .86),
			Color(.5,.5,.5), 20);
		GeomTranspQuads *quads = new GeomTranspQuads;
		const Color black(0,0,0), gray(0.3,0.3,0.3);
		Material *screen = new Material(black, gray, gray, 5);
		screen->transparency = 0.85;

		for(i = 0; i < view_window_->viewwindow_clip_frames_.size(); i++)
		{
			if( view_window_->viewwindow_clip_frames_draw_[i] )
			{
				ViewWindowClipFrame *cf = view_window_->viewwindow_clip_frames_[i];

				for(ii = 0; ii < cf->edges_.size(); ii++)
				{
					frames->add(new GeomSphere(*(cf->corners_[ii])));
					frames->add(new GeomCylinder(*(cf->edges_[ii])));
				}

				// keep all quads in one list so that they are sorted properly.
				quads->add(cf->verts_[0],  screen,
					cf->verts_[1],  screen,
					cf->verts_[2],  screen,
					cf->verts_[3],  screen);

			}
		}
		// add the quads to the frame
		frames->add(quads);

		// turn clip planes off
		std::vector<bool> cliplist(6, false);
		for (ii = 0; ii < 6; ii++)
		{
			if (glIsEnabled((GLenum)(GL_CLIP_PLANE0+ii)))
			{
				glDisable((GLenum)(GL_CLIP_PLANE0+ii));
				cliplist[ii] = true;
			}
		}

		// Put everything in a GeometryHandle.  This should
		// clean up all of the GeomObjects we just made.
		GeomHandle objs = new GeomMaterial(frames, frame_color);
		// now draw it.
		objs->draw(drawinfo_,viewer_->default_material_.get_rep(), current_time_);


		// turn clip planes back on
		for (ii = 0; ii < 6; ii++)
		{
			if (cliplist[ii])
			{
				glEnable((GLenum)(GL_CLIP_PLANE0+ii));
			}
		}
	}

#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT

	void OpenGL::send_message( BioMesh3d::MessageHandle msg )
	{
		this->m_client_communicator_->send_message( msg );
	}

	void OpenGL::close_unconnected_scirun()
	{	
		// if we have never connected then let's exit after some time
		if ( !this->m_client_communicator_ || !this->m_client_communicator_->is_connected() )
		{
			exit( 0 );
		}
	}

#endif

} // End namespace SCIRun
