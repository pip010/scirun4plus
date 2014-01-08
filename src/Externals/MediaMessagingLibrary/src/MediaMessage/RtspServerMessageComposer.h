#ifndef RTSP_SERVER_MESSAGE_COMPOSER_H
#define RTSP_SERVER_MESSAGE_COMPOSER_H

namespace BioMesh3d
{
	class RtspServerMessageComposer
	{
	public:
		RtspServerMessageComposer(){}
		~RtspServerMessageComposer(){}

	protected:
		std::string compose_setup_response_message( 
			const unsigned int session_id,
			const unsigned int cseq,
			const unsigned int client_rtp_port,
			const unsigned int client_rtcp_port,
			const std::string& client_ip,
			const unsigned int server_rtp_port,
			const unsigned int server_rtcp_port,
			const std::string& server_ip
			//const std::string& full_project_name,
			//const unsigned short stage_number
			)
		{
			/*std::string info = full_project_name + "." 
				+ boost::lexical_cast<std::string>(stage_number);*/

			std::string cseq_str = "CSeq: ";
			cseq_str += boost::lexical_cast<std::string> (cseq) + "\r\n";

			std::string session_str = "Session: "; 
			session_str += boost::lexical_cast<std::string> (session_id) +  "\r\n";

			std::string transport_info = "Transport: RTP/UDP;destination=" + 
				client_ip + ";" + "source=" + server_ip + ";" + "client_port=" +
				boost::lexical_cast<std::string> (client_rtp_port) + "-" +
				boost::lexical_cast<std::string> (client_rtcp_port) + ";" +
				"server_port=" + boost::lexical_cast<std::string> (server_rtp_port) +
				"-" + boost::lexical_cast<std::string> (server_rtcp_port) + " " + "\r\n";

			std::string setup_response_msg = "RTSP/1.0 200 OK\r\n";

			setup_response_msg += cseq_str;
			setup_response_msg += session_str;
			setup_response_msg += transport_info;
			setup_response_msg += "\r\n";

			return setup_response_msg;
		}


		std::string compose_teardown_response_message( 
			const unsigned int cseq,
			const unsigned int session_id )
		{
			std::string cseq_str = "CSeq: ";
			cseq_str += boost::lexical_cast<std::string> (cseq) + "\r\n";

			std::string session_str = "Session: ";
			session_str += boost::lexical_cast<std::string> (session_id) +  "\r\n";

			std::string teardown_response_msg = "RTSP/1.0 200 OK\r\n";

			teardown_response_msg += cseq_str;
			teardown_response_msg += session_str;
			teardown_response_msg += "\r\n";

			return teardown_response_msg;
		}


	};

}// end namespace Message

#endif