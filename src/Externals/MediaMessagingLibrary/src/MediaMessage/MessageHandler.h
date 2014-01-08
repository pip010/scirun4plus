#ifndef MESSAGE_HANDLER_H
#define MESSAGE_HANDLER_H

#include <boost/function.hpp>
#include <boost/noncopyable.hpp>
#include <boost/thread.hpp>

#include <map>

#include "Message.h"

namespace BioMesh3d
{
	class MessageHandler : public boost::noncopyable
	{
	public:
		MessageHandler();
		~MessageHandler();

		void handle_message( MessageHandle message )
		{
			std::cerr << "Received message of type " << message->get_message_type() << std::endl;
			boost::function< void ( MessageHandle ) > funct;
			{ // lock m_message_handler_callbacks_ and get the function callback
				boost::mutex::scoped_lock lock( m_map_mutex_ );
				if ( m_message_handler_callbacks_.find( message->get_message_type() ) == m_message_handler_callbacks_.end() )
				{
					return;
				}
				funct = m_message_handler_callbacks_[message->get_message_type()];
			}
			funct( message );
		}

	protected:
		virtual void register_handlers()=0;
		virtual void unregister_handlers()=0;

		void register_callback( MessageTypes::MessageType message_type, boost::function< void ( MessageHandle ) > message_callback );
		void unregister_callback( MessageTypes::MessageType message_type );

	private:
		boost::mutex m_map_mutex_;
		std::map< MessageTypes::MessageType, boost::function< void ( MessageHandle ) > > m_message_handler_callbacks_;
	};
}
#endif