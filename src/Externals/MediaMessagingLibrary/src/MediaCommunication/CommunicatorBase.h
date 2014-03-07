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

#ifndef MEDIA_SERVER_BASE_H
#define MEDIA_SERVER_BASE_H

// Boost includes
#include <boost/function.hpp>
#include <boost/thread/condition.hpp>
#include <boost/noncopyable.hpp>
#include <boost/enable_shared_from_this.hpp>

#include <boost/cstdint.hpp>
#include <boost/unordered_map.hpp>
#include <boost/array.hpp>

#include <Message.h>

// Core includes
#include <Core/Utils/EnumClass.h>
#include <Core/Utils/Lockable.h>

namespace BioMesh3d
{

	namespace NetworkBufferSizes
	{
		enum NetworkBufferSize
		{
			MAX_PACKET_SIZE_E = 10000,
			BUFFER_SIZE_E = 1000000
		};
	}

	namespace MediaCommunicationTypes
	{
		enum MediaCommunicationType
		{
			TCP_E,
			RTP_E
		};
	}

// Forward declarations
class CommunicatorBase;
class CommunicatorBasePrivate;

typedef boost::shared_ptr< CommunicatorBase > CommunicatorBaseHandle;
typedef boost::shared_ptr< CommunicatorBasePrivate > CommunicatorBasePrivateHandle;

class CommunicatorBase : public boost::enable_shared_from_this< CommunicatorBase >,
	public Core::RecursiveLockable
{
public:
	// CLOSE_FUNCTION_TYPE
	// Function that is called when socket is closed
	typedef boost::function< void () > close_function_type;

	// CONNECTED_FUNCTION_TYPE
	// Function that is called when connection is established
	typedef boost::function< void () > connected_function_type;


public:
	CommunicatorBase( std::string ip_address, unsigned short port );
	virtual ~CommunicatorBase();

	bool is_connected();
	bool set_is_connected( bool is_connected );

	
	void set_close_function( close_function_type close_functor );
	void set_connected_function( connected_function_type connected_functor );

	void reset_close_function();
	void reset_connected_function();
	
	static bool get_available_port( unsigned short port_start, 
		unsigned short port_end, unsigned short& available_port, 
		MediaCommunicationTypes::MediaCommunicationType communication_type );

	/*
	 * Run returns immediately. 
	 */
	virtual void run()=0;
	virtual void stop()=0;


protected:

	// GET_IO_SERVICE
	// Get the IO service object, needed for creating sockets and setting up receiver functions
	boost::asio::io_service& get_io_service();

	void call_communication_close_functor();
	void call_communication_connected_functor();
	virtual void register_receive()=0;

	bool m_is_connected;
	boost::mutex m_is_connected_lock;

	std::string m_server_ip_address;
	unsigned short m_port;
	bool m_is_client_;
	boost::uint16_t m_packet_size_;

	boost::array< unsigned char, NetworkBufferSizes::BUFFER_SIZE_E > m_buf_;

	boost::recursive_mutex m_stop_mutex_;
	bool m_stop_;

	boost::mutex m_send_message_mutex_;
	boost::condition m_send_message_condition_;
	bool m_message_sent_;

private:
	CommunicatorBasePrivateHandle private_;
};

} // end namespace MediaCommunication


#endif