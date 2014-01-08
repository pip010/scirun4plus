#ifndef MESSAGE_PARSER_H
#define MESSAGE_PARSER_H

#include <boost/noncopyable.hpp>
#include "Message.h"

namespace BioMesh3d
{
	class MessageParser : public boost::noncopyable
	{
	public:

		static MessageTypes::MessageType get_message_type( boost::uint8_t* buf );
		static void set_message_type( boost::uint8_t* buf, MessageTypes::MessageType message_type );
		static boost::uint32_t get_message_size( boost::uint8_t* buf );
		static void set_message_size( boost::uint8_t* buf, boost::uint32_t message_size );

		// TODO: Move to cc file
		static MessageHandle parse_message( unsigned char* buffer )
		{
			MessageHandle message = MessageHandle( new Message() );
			//message->reset_message();
			message->set_message_type( MessageParser::get_message_type( buffer ) );
			std::string params = std::string( reinterpret_cast< char* >( buffer + MESSAGE_HEADER_SIZE ) );
			if ( !params.empty() )
			{
				std::vector< std::string > param_list;
				boost::algorithm::split( param_list, params, boost::algorithm::is_any_of( "&=" ) );
				for ( unsigned int i = 0; i < param_list.size(); i += 2 )
				{
					if ( i + 1 < param_list.size() )
					{
						try
						{
							message->add_message_parameter( boost::lexical_cast< int >( param_list[i] ), param_list[i + 1] );
						}
						catch ( boost::bad_lexical_cast cast_error )
						{
							message->reset_message();
							message->set_message_type( MessageTypes::INVALID_E );
							std::cout << "There was a problem parsing the message: " << 
								cast_error.what() << std::endl;
						}
					}
				}
			}
			
			// TODO: Create buffer and insert into message
			unsigned int binary_data_size =
				MessageParser::get_message_size( buffer ) - 
				( MESSAGE_HEADER_SIZE + params.size() + 1 );
			if ( 0 < binary_data_size )
			{
				message->copy_binary_data( buffer + MESSAGE_HEADER_SIZE + params.size() + 1, binary_data_size );
			}
			return message;
		}
	private:
		template < class T >
		static void retrieve_header_parameter( boost::uint8_t* buf, HeaderOffsets::Offset offset, T& value )
		{
			value = *reinterpret_cast< T* >( buf + offset );
		}

		template < class T >
		static void add_header_parameter( boost::uint8_t* buf, HeaderOffsets::Offset offset, T value )
		{
			*reinterpret_cast< T* >( buf + offset ) = value;
		}
	};
}

#endif
