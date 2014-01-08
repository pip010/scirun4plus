#ifndef TCP_MEDIA_SERVER_H
#define TCP_MEDIA_SERVER_H

#include <boost/asio.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/signals2.hpp>
#include <boost/thread/condition.hpp>

#include <string>
#include <map>

#include "CommunicatorBase.h"
#include <Message.h>


namespace BioMesh3d
{
	class TCPCommunicator;
	typedef boost::shared_ptr< TCPCommunicator > TCPHandle;

	typedef boost::shared_ptr< boost::asio::ip::tcp::socket > TCPSocketHandle;
	typedef boost::shared_ptr< boost::asio::ip::tcp::acceptor > TCPAcceptorHandle;

	class TCPCommunicator : public CommunicatorBase
	{
	public:
		TCPCommunicator( std::string server_ip_address, unsigned short server_port );
		~TCPCommunicator();

		virtual void run();
		void stop();

		std::string get_endpoint_ip_address();

	protected:
		void send_bytes( unsigned char* bytes, size_t bytes_size );

		TCPSocketHandle m_socket_handle_;

	private:
		void handle_connect( const boost::system::error_code& error,
			boost::asio::ip::tcp::resolver::iterator endpoint_iterator );
		void handle_accept( const boost::system::error_code& error );

		static boost::mutex s_tcp_acceptor_map_mutex_;
		static std::map< unsigned short, TCPAcceptorHandle > s_tcp_acceptor_map_;
	};
}// end namespace MediaCommunication

#endif