#include "MessageParser.h"

namespace BioMesh3d
{


	MessageTypes::MessageType MessageParser::get_message_type( boost::uint8_t* buf )
	{
		unsigned short message_type = 0;
		retrieve_header_parameter( buf, HeaderOffsets::TYPE_E, message_type );
		message_type = ntohs( message_type );
		return static_cast< MessageTypes::MessageType >( message_type );
	}

	void MessageParser::set_message_type( boost::uint8_t* buf, MessageTypes::MessageType message_type )
	{
		MessageTypes::MessageType msg_type = ( 0 < message_type && message_type < MessageTypes::MESSAGE_TYPE_SIZE_E ) 
			? message_type : MessageTypes::INVALID_E;
		unsigned short temp_type = static_cast< MessageTypes::MessageType >( htons( msg_type ) );
		add_header_parameter( buf, HeaderOffsets::TYPE_E, temp_type );
	}

	boost::uint32_t MessageParser::get_message_size( boost::uint8_t* buf )
	{
		boost::uint32_t message_size = 0;
		retrieve_header_parameter( buf, HeaderOffsets::SIZE_E, message_size );
		message_size = ntohl( message_size );
		return static_cast< boost::uint32_t >( message_size );
	}

	void MessageParser::set_message_size( boost::uint8_t* buf, boost::uint32_t message_size )
	{
		boost::uint32_t temp_msg_size = htonl( message_size );
		add_header_parameter( buf, HeaderOffsets::SIZE_E, 
			static_cast< boost::uint32_t >( temp_msg_size ) );
	}

} // end namespace Message

