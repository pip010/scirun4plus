#ifndef TCP_MESSAGE_COMMUNICATOR
#define TCP_MESSAGE_COMMUNICATOR

#include <boost/shared_ptr.hpp>

#include "TCPCommunicator.h"

namespace BioMesh3d
{
	class TCPMessageCommunicator;
	typedef boost::shared_ptr< TCPMessageCommunicator > TCPMessageCommunicatorHandle;

	class TCPMessageCommunicator : public TCPCommunicator
	{
	public:
		TCPMessageCommunicator( std::string ip_address, unsigned short port, 
			boost::function< void ( MessageHandle ) > handle_message_functor );
		~TCPMessageCommunicator();

		void send_message( MessageHandle msg );
	private:
		virtual void register_receive();
		void handle_write( const boost::system::error_code& error );
		void handle_receive( const boost::system::error_code& error, size_t bytes_transferred );

		boost::function< void ( MessageHandle ) > m_handle_message_functor_;
		std::vector< unsigned char > buffer_;
	};
}

#endif