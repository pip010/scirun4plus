#include <sstream>

#include "MessageHandler.h"
#include "Message.h"

namespace BioMesh3d
{
	MessageHandler::MessageHandler()
	{
	}

	MessageHandler::~MessageHandler()
	{
	}

	void MessageHandler::register_callback( MessageTypes::MessageType message_type, 
		boost::function< void ( MessageHandle ) > message_callback )
	{
		boost::mutex::scoped_lock lock( m_map_mutex_ );
		m_message_handler_callbacks_.insert(
			std::pair<MessageTypes::MessageType, boost::function< void ( MessageHandle ) > >(message_type, message_callback) );
	}

	void MessageHandler::unregister_callback( MessageTypes::MessageType message_type )
	{
		boost::mutex::scoped_lock lock( m_map_mutex_ );
		std::map< MessageTypes::MessageType, boost::function< void ( MessageHandle ) > >::iterator iter;
		iter = m_message_handler_callbacks_.find( message_type );
		if ( iter != m_message_handler_callbacks_.end() )
		{
			m_message_handler_callbacks_.erase( iter );
		}
	}
}// end namespace Message