#ifndef RTSP_MEDIA_COMMUNICATION_H
#define RTSP_MEDIA_COMMUNICATION_H

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "TCPCommunicator.h"

namespace BioMesh3d
{
	class RTSPMediaCommunicator;
	typedef boost::shared_ptr< RTSPMediaCommunicator > RTSPHandle;

	class RTSPMediaCommunicator : public TCPCommunicator
	{
	public:
		RTSPMediaCommunicator();
		RTSPMediaCommunicator( const std::string& ip_address );
		~RTSPMediaCommunicator();

		virtual void run();
		void send_biomesh_visualization_setup( const std::string& label_map, const std::string& mesh, 
			const std::string& username, unsigned short stage_number, 
			MediaCommunicationTypes::MediaCommunicationType com_type,
			boost::function< void ( unsigned short, unsigned short ) > biomesh_ok_callback ); // client uses this

		void send_ok( unsigned short scirun_port ); //server uses this
		void set_setup_callback( boost::function< void (  std::string, std::string, unsigned short ) > setup_callback_ );
	private:
		void send_teardown();
		void send_error();

		void handle_receive( const boost::system::error_code& error, size_t bytes_transferred );
		void handle_write( const boost::system::error_code& error ) {}
		bool parse_setup( const std::string& rtsp_message );
		bool parse_ok( const std::string& rtsp_message );
		void register_receive();
		
		unsigned short m_client_rtp_port_;
		unsigned short m_client_rtcp_port_;
		unsigned short m_scirun_rtp_port_;
		unsigned short m_scirun_rtcp_port_;

		std::string m_label_map_name_;
		std::string m_mesh_name_;
		std::string m_username_;
		std::string m_ip_address_;
		unsigned short m_stage_number_;
		MediaCommunicationTypes::MediaCommunicationType m_communication_type_;

		static unsigned int s_c_seq_;
		unsigned int m_session_;

		TCPHandle m_tcp_communicator_;

		boost::function< void ( unsigned short, unsigned short ) > m_ok_callback_;
		boost::function< void (  std::string, std::string, unsigned short ) > m_setup_callback_;
		boost::asio::streambuf m_stream_data_;
	};
} // end namespace 
#endif

//
//	class RtspMediaCommunication : public TCPMediaCommunication
//	{
//	public:
//		RtspMediaCommunication( const std::string& label_map_name, const std::string& mesh_name,
//			 unsigned short stage_number, const std::string& ip_address );
//		RtspMediaCommunication( boost::function< void ( std::string, unsigned short, 
//			std::string, unsigned short& ) > func );
//		~RtspMediaCommunication();
//
//		bool pinhole_firewall_get_port( unsigned short& client_rtp_port, 
//			unsigned short& client_rtcp_port, unsigned short& server_rtp_port,
//			unsigned short& server_rtcp_port );
//
//	private:
//		void register_callbacks();
//		static void handle_receive_rtsp( const boost::system::error_code& e_code, 
//			size_t bytes_transferred, RtspMediaCommunication* communicator );
//
//		void setup_rtsp( const std::string& rtsp_message );
//
//		static int s_cseq_;
//
//		boost::asio::streambuf m_streambuf_;
//		unsigned short m_local_port_;
//
//		std::string m_label_map_name_;
//		std::string m_mesh_name_;
//		unsigned short m_stage_number_;
//
//		boost::mutex m_response_mutex_;
//		boost::condition m_response_condition_;
//		std::string m_response_string_;
//
//		boost::mutex m_rtsp_session_id_mutex_;
//		int m_rtsp_session_id_;
//
//		boost::function< void ( std::string, unsigned short, 
//			std::string, unsigned short& ) > m_setup_callback_;
//	}; 
//

