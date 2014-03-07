#ifndef RTSP_CLIENT_MESSAGE_HANDLER_H
#define RTSP_CLIENT_MESSAGE_HANDLER_H

namespace BioMesh3d
{

	class RtspClientMessageHandler 
	{
	public:
		RtspClientMessageHandler(){}
		~RtspClientMessageHandler(){}

	protected:
		//virtual void handle_setup_request( const std::string& setup_request );
		//virtual void handle_teardown_request( const std::string& teardown_response );
		virtual void handle_setup_request( boost::asio::streambuf& buf );
		virtual void handle_teardown_request( boost::asio::streambuf& buf );
		
	private:
		//void register_handlers();
		//void unregister_handlers();

	};
}

#endif