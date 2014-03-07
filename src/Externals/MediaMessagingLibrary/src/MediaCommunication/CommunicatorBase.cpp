/*
 For more information, please see: http://software.sci.utah.edu

 The MIT License

 Copyright (c) 2011 Scientific Computing and Imaging Institute,
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
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>


#include "CommunicatorBase.h"

namespace BioMesh3d
{

class CommunicatorBasePrivate
{
	// -- member variables --
public:
	// IO Service Object: that manages all io interaction
	boost::asio::io_service io_service_;
	
	// An object that records whether there still is work to be done, so the threads
	// will not exit. The destructor of this object will signal the io_service object
	// that all threads should end when done.
	boost::asio::io_service::work* io_service_work_;
	
	// Mutex protecting this class.
	boost::mutex thread_mutex_;
	
	// Condition variable for checking whether the receiving thread has done
	// We need to make sure this is done before we can destroy the private class
	// As the thread needs to have ended before io_servive can be released.
	boost::condition_variable receive_thread_done_;

	// Condition variable that signals that thread has been started
	boost::condition_variable thread_started_;
	
	// Amount of threads for running io service
	int thread_count_;
	
	// Function that is called when the connection is closed
	CommunicatorBase::close_function_type close_function_;
	
	// Function that is called when the connection is established
	CommunicatorBase::connected_function_type connected_function_;
	
	// -- private functions --
public:

	// Function that is executed on the receive thread
	void receive_thread();
};

void CommunicatorBasePrivate::receive_thread()
{	
	{
		// Need to lock here to ensure that the other thread is waiting
		boost::mutex::scoped_lock lock( this->thread_mutex_ );
	
		this->thread_started_.notify_all();
	}

	// Start the asynchronous IO thread
	this->io_service_.run();

	{
		// Need to lock here to ensure that the other thread is waiting
		boost::mutex::scoped_lock lock( this->thread_mutex_ );
		
		this->thread_count_--;
		
		// Wake up the thread that is waiting for this one to end
		receive_thread_done_.notify_all();
	}
}

CommunicatorBase::CommunicatorBase( std::string ip_address, unsigned short port ) : 
	m_is_connected( false ),
	m_server_ip_address( ip_address ), 
	m_port( port ), 
	m_is_client_( ip_address.size() > 0 ),
	m_packet_size_( NetworkBufferSizes::MAX_PACKET_SIZE_E ), 
	m_stop_( false ),
	// TODO: Need to remove this variable
	m_message_sent_(false),
	private_( new CommunicatorBasePrivate )
{	
	// Tell the io_service object that there is still work to be done
	this->private_->io_service_work_ = new boost::asio::io_service::work( 
		this->private_->io_service_ );
	
	this->private_->thread_count_ = 0;

	// This lock will ensure that the counter is properly protected
	boost::mutex::scoped_lock lock( this->private_->thread_mutex_ );
	
	// Create threads to handle the messages
	for ( int j = 0; j < 2; j++ )
	{
		boost::thread( boost::bind(
			&CommunicatorBasePrivate::receive_thread, this->private_ ) );
		this->private_->thread_count_++;
		
		this->private_->thread_started_.wait( lock );
	}
}	

CommunicatorBase::~CommunicatorBase()
{
	// Need to lock before killing the thread
	boost::mutex::scoped_lock thread_lock( this->private_->thread_mutex_ );

	// Delete the work object to signal that we no longer need the thread
	delete this->private_->io_service_work_;

	// Wait for the thread to finish so we can clean up the private class
	while ( this->private_->thread_count_ > 0 )
	{
		this->private_->receive_thread_done_.wait( thread_lock );
	}
}


boost::asio::io_service& CommunicatorBase::get_io_service()
{
	return this->private_->io_service_;
}



bool CommunicatorBase::get_available_port( unsigned short port_start, unsigned short port_end,
	unsigned short& available_port, MediaCommunicationTypes::MediaCommunicationType communication_type )
{
	boost::asio::io_service test_io_service;

	available_port = port_start;
	available_port = ( communication_type == MediaCommunicationTypes::RTP_E && 
		available_port % 2 != 0) ? 
		available_port + 1 : available_port;
	while ( available_port <= port_end )
	{
		try
		{
			if ( communication_type == MediaCommunicationTypes::TCP_E )
			{
				boost::asio::ip::tcp::socket test_socket(test_io_service, 
					boost::asio::ip::tcp::endpoint( boost::asio::ip::tcp::v4(), 
					available_port ) );
				test_socket.close();
			}
			else if ( communication_type == MediaCommunicationTypes::RTP_E )
			{
				boost::asio::ip::udp::socket test_socket(test_io_service, 
					boost::asio::ip::udp::endpoint( boost::asio::ip::udp::v4(), 
					available_port ) );
				test_socket.close();
			}
			return true;
		}
		catch ( std::exception& e )
		{
			std::string exception_string = e.what();
			available_port = ( communication_type == MediaCommunicationTypes::RTP_E ) ? 
				available_port + 2 : available_port + 1;
		}
	}
	return false;
}

bool CommunicatorBase::is_connected()
{
	lock_type lock( this->get_mutex() );
	return this->m_is_connected;
}

bool CommunicatorBase::set_is_connected( bool is_connected )
{
	boost::mutex::scoped_lock lock( this->m_is_connected_lock );
	this->m_is_connected = is_connected;
	return is_connected;
}

void CommunicatorBase::set_close_function( close_function_type close_function )
{
	lock_type lock( this->get_mutex() );
	this->private_->close_function_ = close_function;
}

void CommunicatorBase::set_connected_function( connected_function_type connected_function )
{
	lock_type lock( this->get_mutex() );
	this->private_->connected_function_ = connected_function;
}

void CommunicatorBase::reset_close_function()
{
	lock_type lock( this->get_mutex() );
	this->private_->close_function_.clear();
}

void CommunicatorBase::reset_connected_function()
{
	lock_type lock( this->get_mutex() );
	this->private_->connected_function_.clear();
}



void CommunicatorBase::call_communication_close_functor()
{
	close_function_type close_function;
	{
		// While the system is locked grab the call back function
		// and remove the function
		lock_type lock( this->get_mutex() );
		close_function = this->private_->close_function_;
		this->private_->close_function_.clear();
	}

	// If the function is defined call it.
	if ( close_function )
	{
		close_function();
	}
}

void CommunicatorBase::call_communication_connected_functor()
{
	connected_function_type connected_function;
	{
		// While the system is locked grab the call back function
		// and remove the function
		lock_type lock( this->get_mutex() );
		connected_function = this->private_->connected_function_;
		this->private_->connected_function_.clear();
	}

	// If the function is defined call it.
	if ( connected_function )
	{
		connected_function();
	}
}

}
