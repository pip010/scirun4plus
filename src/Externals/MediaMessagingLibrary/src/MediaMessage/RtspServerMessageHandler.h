#ifndef RTSP_SERVER_MESSAGE_HANDLER_H
#define RTSP_SERVER_MESSAGE_HANDLER_H

#include <boost/asio.hpp>

namespace BioMesh3d
{
	class RtspServerMessageHandler 
	{
	public:
		RtspServerMessageHandler(){}
		~RtspServerMessageHandler(){}

	protected:
		virtual void handle_setup_response( boost::asio::streambuf& buf );
		virtual void handle_teardown_response( boost::asio::streambuf& buf );
			
	private:

	};

}

#endif