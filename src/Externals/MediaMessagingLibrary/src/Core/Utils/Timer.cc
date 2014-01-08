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

// Boost includes
#include <boost/thread.hpp>
#include <boost/signals2.hpp>

// Core includes
#include <Core/Utils/Timer.h>
#include <Core/Utils/Log.h>

namespace Core
{


class TimerPrivate
{
	// Member variables
public:
	// Pointer back to main thread
	Timer* timer_;

	// Mutex for timer variables
	boost::mutex timer_mutex_;
	
	// CV for signaling thread that the state has changed
	boost::condition_variable signal_timer_wait_;

	// CV signalling that thread was properly destroyed
	boost::condition_variable timer_done_wait_;

	// Whether timer is on
	bool timer_on_;

	// Whether the timer is done and needs to be desctructed
	bool timer_done_;
	
	// Whether the timer thread is done
	bool timer_thread_done_;

	// The time out
	double timer_wait_;
	
	// The thread running the times
	boost::shared_ptr<boost::thread> thread_;
	
public:
	TimerPrivate() :
		timer_( 0 ),
		timer_on_( false ),
		timer_done_( false ),
		timer_thread_done_( false ),
		timer_wait_( 10.0 )
	{}
	
	
public:
	// RUN_THREAD:
	// This function is the main logic that the separate timer thread use, it only wakes up
	// when it has to trigger the next event.
	void run_thread();

};

void TimerPrivate::run_thread()
{
	CORE_LOG_DEBUG( "Started timer thread." );
	boost::mutex::scoped_lock lock( this->timer_mutex_ );

	for ( ;; )
	{
		if ( this->timer_done_ ) 
		{
			// Signal that we are about to leave the thread
			// A signal the waiting thread that is destroying the timer
			// The main thread will only be release when the lock is released
			// when this function goes out of scope.
			this->timer_thread_done_ = true;
			this->timer_done_wait_.notify_all();
			break;
		}
	
		if ( this->timer_on_ )
		{
			boost::posix_time::time_duration duration = 
				boost::posix_time::seconds( this->timer_wait_ ); 
			this->signal_timer_wait_.timed_wait(lock,  duration );
			
			if ( this->timer_on_ )
			{
				CORE_LOG_DEBUG( "Trigger timer." );
				this->timer_->timer_signal_();
			}
		}
		else
		{
			this->signal_timer_wait_.wait( lock );
		}
	
	}
	CORE_LOG_DEBUG( "Finished timer thread." );
}

Timer::Timer() :
	private_( new TimerPrivate )
{
	this->private_->timer_ = this;
}

Timer::~Timer()
{
	// Kill the thread safely
	if ( this->private_->thread_ )
	{
		boost::mutex::scoped_lock lock( this->private_->timer_mutex_ );
		
		// Tell the thread to stop
		this->private_->timer_done_ = true;
		this->private_->signal_timer_wait_.notify_all();

		while ( this->private_->timer_thread_done_ )
		{
			this->private_->timer_done_wait_.wait( lock );
		}
	}
}

bool Timer::start( double interval )
{
	boost::mutex::scoped_lock lock( this->private_->timer_mutex_ );
	
	if ( !this->private_->thread_ )
	{
		this->private_->timer_on_ = true;
		this->private_->timer_wait_ = interval;
		
		this->private_->thread_ = boost::shared_ptr< boost::thread>(
			new boost::thread( boost::bind( &TimerPrivate::run_thread, this->private_ ) ) );
	}
	else
	{
		this->private_->timer_on_ = true;
		this->private_->timer_wait_ = interval;
		this->private_->signal_timer_wait_.notify_all();
	}
	
	return true;
}
	
	
bool Timer::stop()
{	
	boost::mutex::scoped_lock lock( this->private_->timer_mutex_ );
	this->private_->timer_on_ = false;
	this->private_->signal_timer_wait_.notify_all();	

	return true;
}

} // end namespace Core
