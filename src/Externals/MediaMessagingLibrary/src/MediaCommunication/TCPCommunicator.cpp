#include <boost/asio.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <string>
#include <iostream>
#include <fstream>

#include "TCPCommunicator.h"
#include <Message.h>
#include <MessageParser.h>


namespace BioMesh3d
{
	std::map< unsigned short, TCPAcceptorHandle > TCPCommunicator::s_tcp_acceptor_map_;
	boost::mutex TCPCommunicator::s_tcp_acceptor_map_mutex_;

	TCPCommunicator::TCPCommunicator( std::string ip_address, unsigned short port ) :
		CommunicatorBase( ip_address, port )
	{
		boost::mutex::scoped_lock acceptor_lock( s_tcp_acceptor_map_mutex_ );

		// we have to have a separate acceptor for each port we are listening on, if there isn't an acceptor on that port,
		// then create on and add it to the map using the port number as the key
		if ( !m_is_client_ && s_tcp_acceptor_map_.find( port ) == s_tcp_acceptor_map_.end() )
		{
			// Construct an acceptor opened on the given endpoint. Creates an 
			// acceptor and automatically opens it to listen for new connections 
			// on the specified endpoint.
			s_tcp_acceptor_map_[port] = TCPAcceptorHandle( new boost::asio::ip::tcp::acceptor( 
				this->get_io_service(), 
				boost::asio::ip::tcp::endpoint( boost::asio::ip::tcp::v4(), port ) ) );
		}
		this->m_socket_handle_ = TCPSocketHandle( 
			new boost::asio::ip::tcp::socket( this->get_io_service() ) );
		set_is_connected( false );
	}

	TCPCommunicator::~TCPCommunicator()
	{
	}

	std::string TCPCommunicator::get_endpoint_ip_address()
	{
		std::string address = "";
		if ( this->m_socket_handle_ )
		{
			boost::system::error_code ec;
			address = this->m_socket_handle_->remote_endpoint( ec ).address().to_string();
		}
		return address;
	}

	void TCPCommunicator::stop() // stop all socket operations, all the async functions will return and when all the shared pointers of this are gone, the destructor will be called
	{
		boost::recursive_mutex::scoped_lock lock( this->m_stop_mutex_ );
		if ( !this->m_stop_ )
		{
			this->m_stop_ = true;
			set_is_connected( false );
			call_communication_close_functor();
			this->reset_connected_function();
			this->reset_close_function();
			this->m_socket_handle_->close();
		}
	}

	void TCPCommunicator::run()
	{
		if ( this->m_is_client_ ) // if we are the client, let's use DNS and try and connect
		{
			boost::asio::ip::tcp::resolver resolver( this->m_socket_handle_->get_io_service() );
			boost::asio::ip::tcp::resolver::query query( this->m_server_ip_address, 
				boost::lexical_cast< std::string >( this->m_port ) );
			boost::system::error_code e_code;
			boost::asio::ip::tcp::resolver::iterator iterator = resolver.resolve( query, e_code );
			if ( e_code ) { return; }
			boost::asio::ip::tcp::endpoint endpoint = *iterator;
			this->m_socket_handle_->async_connect( endpoint, boost::bind( 
				&TCPCommunicator::handle_connect, 
				boost::dynamic_pointer_cast< TCPCommunicator >( shared_from_this() ), 
				boost::asio::placeholders::error, ++iterator) );
		}
		else // if we are a server, we register an accept callback and wait
		{
			boost::mutex::scoped_lock acceptor_lock( s_tcp_acceptor_map_mutex_ );
			TCPAcceptorHandle acceptor_handle = s_tcp_acceptor_map_[this->m_port];
			acceptor_handle->async_accept( *this->m_socket_handle_, 
				boost::bind( &TCPCommunicator::handle_accept, 
				boost::dynamic_pointer_cast< TCPCommunicator >( shared_from_this() ), _1 ) );
		}
	}

	void TCPCommunicator::handle_connect( const boost::system::error_code& error, 
		boost::asio::ip::tcp::resolver::iterator endpoint_iterator )
	{
		if ( !error ) // if everything is fine, then register the receive function
		{
			register_receive();
			set_is_connected( true );
			call_communication_connected_functor();
		}
		else if ( endpoint_iterator != boost::asio::ip::tcp::resolver::iterator() ) // that iterator failed, try the next one
		{
			this->m_socket_handle_->close();
			boost::asio::ip::tcp::endpoint endpoint = *endpoint_iterator;
			this->m_socket_handle_->async_connect( endpoint, boost::bind( 
				&TCPCommunicator::handle_connect, 
				boost::dynamic_pointer_cast< TCPCommunicator >( shared_from_this() ), 
				boost::asio::placeholders::error, ++endpoint_iterator ) );
		}
		else // if there is an error (and all enpoints fail) then do the callback (if any) and return without registering another function
		{
			call_communication_close_functor();
		}
	}


	void TCPCommunicator::handle_accept( const boost::system::error_code& error )
	{	
	
		if ( error )// if there is an error then do the callback (if any) and return without registering another function
		{
			call_communication_close_functor();
			return;
		}
		register_receive();
		set_is_connected( true );
		call_communication_connected_functor();
	}

}// end namespace MediaServer
