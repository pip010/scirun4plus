#include <MessageParser.h>

#include "TCPMessageCommunicator.h"

namespace BioMesh3d
{


	TCPMessageCommunicator::TCPMessageCommunicator( std::string ip_address, unsigned short port, 
		boost::function< void ( MessageHandle ) > handle_message_functor ) :
		TCPCommunicator( ip_address, port ),
		m_handle_message_functor_( handle_message_functor )
	{
		this->buffer_.reserve( 1 << 20 );
		this->buffer_.resize( sizeof( boost::uint32_t ) );
	}

	TCPMessageCommunicator::~TCPMessageCommunicator()
	{
	}

	void TCPMessageCommunicator::register_receive()
	{
		boost::asio::async_read( *this->m_socket_handle_, 
			boost::asio::buffer( &this->buffer_[ 0 ], sizeof( boost::uint32_t ) ), 
			boost::asio::transfer_all(),
			boost::bind( &TCPMessageCommunicator::handle_receive, 
			boost::dynamic_pointer_cast< TCPMessageCommunicator >( shared_from_this() ), _1, _2 ) );
	}

	void TCPMessageCommunicator::handle_receive( const boost::system::error_code& error, 
		size_t bytes_transferred )
	{
		if ( error ) // if there is an error then do the callback (if any) and return without registering another receive
		{
			stop();
			return;
		}
		if ( bytes_transferred == sizeof( boost::uint32_t ) ) // We need to read the first 4 bytes so we know how many more bytes there are in this message
		{
			unsigned int message_size = MessageParser::get_message_size( &this->buffer_[ 0 ] );
			this->buffer_.resize( message_size );
			boost::asio::async_read( *this->m_socket_handle_, 
				boost::asio::buffer( ( &this->buffer_[ 0 ] + sizeof( boost::uint32_t ) ), 
				( message_size - sizeof( boost::uint32_t ) ) ), boost::asio::transfer_all(),
				boost::bind( &TCPMessageCommunicator::handle_receive, 
				boost::dynamic_pointer_cast< TCPMessageCommunicator >( shared_from_this() ), _1, _2 ) ); // we register the async_rec with the message size - sizeof 
		}
		else // now we have read the entire message object and it's in the m_receive_buffer_, we parse it and call the callback on it.
		{
			MessageHandle msg = MessageParser::parse_message( &this->buffer_[ 0 ] );
			this->m_handle_message_functor_( msg ); // we are calling the callback directly so we might later fix this so we have another thread that handles this.
			boost::asio::async_read( *this->m_socket_handle_, 
				boost::asio::buffer( &this->buffer_[ 0 ] , sizeof( boost::uint32_t ) ), 
				boost::asio::transfer_all(),
				boost::bind( &TCPMessageCommunicator::handle_receive, 
				boost::dynamic_pointer_cast< TCPMessageCommunicator >( shared_from_this() ), _1, _2 ) );
		}
	}

	void TCPMessageCommunicator::send_message( MessageHandle msg )
	{
		std::string serialized_params = msg->serialize_parameters();
		unsigned int send_size =  MESSAGE_HEADER_SIZE + serialized_params.size() + 1 + 
			msg->get_binary_size();

		std::vector< boost::uint8_t > raw_buffer;
		raw_buffer.resize( send_size );

		// the message size is the header size + the param string + binary size + 1
		MessageParser::set_message_size( &raw_buffer[0], send_size );
		MessageParser::set_message_type( &raw_buffer[0], msg->get_message_type() );
		memcpy( &raw_buffer[0] + MESSAGE_HEADER_SIZE, serialized_params.c_str(), 
			serialized_params.size() );
		raw_buffer[MESSAGE_HEADER_SIZE + serialized_params.size()] = 0;

		// TODO: Clean this up a bit
		memcpy( reinterpret_cast< char* >( &raw_buffer[0] + MESSAGE_HEADER_SIZE + 
			serialized_params.size() + 1 ), 
			msg->get_binary_buf(), msg->get_binary_size() );

		boost::mutex::scoped_lock send_lock( this->m_send_message_mutex_ ); // we want to lock this and wait so the buffer doesn't go out of scope
		boost::asio::async_write( *this->m_socket_handle_, 
			boost::asio::buffer( &raw_buffer[0], send_size ), 
			boost::bind( &TCPMessageCommunicator::handle_write, 
			boost::dynamic_pointer_cast< TCPMessageCommunicator >( shared_from_this() ), _1 ) );

		this->m_send_message_condition_.wait( send_lock );
	}

	void TCPMessageCommunicator::handle_write( const boost::system::error_code& error )
	{
		boost::mutex::scoped_lock lock( this->m_send_message_mutex_ );
		this->m_send_message_condition_.notify_all();
	}

}