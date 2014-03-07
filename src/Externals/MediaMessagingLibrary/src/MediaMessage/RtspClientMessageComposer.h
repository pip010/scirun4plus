#ifndef RTSP_CLIENT_MESSAGE_COMPOSER_H
#define RTSP_CLIENT_MESSAGE_COMPOSER_H

namespace BioMesh3d
{
	class RtspClientMessageComposer
	{
	public:
		RtspClientMessageComposer(){}
		~RtspClientMessageComposer(){}


		std::string compose_setup_request_message( 
			const unsigned int cseq,
			const unsigned int client_rtp_port,
			const unsigned int client_rtcp_port,
			const unsigned int server_rtp_port,
			const unsigned int server_rtcp_port,
			const std::string& full_project_name,
			const unsigned short stage_number )
			{
				std::string info = full_project_name + "." 
					+ boost::lexical_cast<std::string>(stage_number);

				std::string cseq_str = "CSeq: ";
				cseq_str += (boost::lexical_cast<std::string> (cseq)) + "\r\n";

				std::string transport_info = "Transport: RTP/UDP;client_port=" +
					boost::lexical_cast<std::string> (client_rtp_port) + "-" +
					boost::lexical_cast<std::string> (client_rtcp_port) + "\r\n";

				std::string user_agent_str = "User_Agent: ";
				user_agent_str +=
					boost::lexical_cast<std::string> (server_rtp_port) +
					"-" + boost::lexical_cast<std::string> (server_rtcp_port) +
					"\r\n";

				std::string setup_request_msg = "SETUP " + info + " "
					+ "RTSP/1.0" + "\r\n";

				setup_request_msg += cseq_str;
				setup_request_msg += transport_info;
				setup_request_msg += user_agent_str;
				setup_request_msg += "\r\n";

				return setup_request_msg;
			}


		std::string compose_teardown_request_message( 
			const unsigned int cseq,
			const unsigned int session_id )
		{
			std::string cseq_str = "CSeq: ";
			cseq_str += boost::lexical_cast<std::string> (cseq) + "\r\n";

			std::string session_str = "Session: ";
			session_str += boost::lexical_cast<std::string> (session_id) +  "\r\n";

			std::string teardown_request_msg = "TEARDOWN RTSP/1.0\r\n";

			teardown_request_msg += cseq_str;
			teardown_request_msg += session_str;
			teardown_request_msg += "\r\n";

			return teardown_request_msg;
		}


	};

}// end namespace Message

#endif