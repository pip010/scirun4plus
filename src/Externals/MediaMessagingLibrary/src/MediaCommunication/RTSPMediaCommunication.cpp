#include <boost/bind.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "RTSPMediaCommunication.h"

namespace BioMesh3d
{

	unsigned int RTSPMediaCommunicator::s_c_seq_ = 0;


	RTSPMediaCommunicator::RTSPMediaCommunicator() :
		TCPCommunicator( "", 8554 )
	{

	}


	RTSPMediaCommunicator::RTSPMediaCommunicator( const std::string& ip_address ) :
		TCPCommunicator( ip_address, 8554 )
	{
		++s_c_seq_;
	}

	RTSPMediaCommunicator::~RTSPMediaCommunicator()
	{
	}

	void RTSPMediaCommunicator::run()
	{
		TCPCommunicator::run();
	}

	void RTSPMediaCommunicator::send_biomesh_visualization_setup( const std::string& label_map, const std::string& mesh, const std::string& username, 
		unsigned short stage_number, MediaCommunicationTypes::MediaCommunicationType com_type,
		boost::function< void ( unsigned short, unsigned short ) > biomesh_ok_callback )
	{
		this->m_ok_callback_ = biomesh_ok_callback;
		this->m_client_rtp_port_ = 35123;
		CommunicatorBase::get_available_port( 35123, 65000, this->m_client_rtp_port_, com_type );
		this->m_client_rtcp_port_ = this->m_client_rtp_port_ + 1;

		std::string url_message = label_map + "?" + mesh + "?" + username + "?" + boost::lexical_cast< std::string >( stage_number ) + "?" + 
			boost::lexical_cast< std::string >( (int)com_type ) + "?";

		std::string setup_message = "SETUP " + url_message + " RTSP/1.0\r\n" 
			"CSeq: " + boost::lexical_cast< std::string >( this->s_c_seq_ ) + "\r\n"
			"Transport: RTP/AVP;unicast;client_port=" +	boost::lexical_cast< std::string >( this->m_client_rtp_port_ ) + "-" +
			boost::lexical_cast< std::string >( this->m_client_rtcp_port_ ) + "\r\n\r\n";
		boost::asio::async_write( *this->m_socket_handle_, 
			boost::asio::buffer( setup_message.c_str(), setup_message.size() ), 
			boost::bind( &RTSPMediaCommunicator::handle_write, boost::dynamic_pointer_cast< RTSPMediaCommunicator >( shared_from_this() ), _1 ) );
	}

	void RTSPMediaCommunicator::send_ok( unsigned short scirun_port )
	{
		this->m_scirun_rtp_port_ = scirun_port;
		this->m_scirun_rtcp_port_ = scirun_port + 1;
		time_t t;
		time( &t );
		char* time_char = asctime( gmtime( &t ) );
		std::string time_string( time_char );
		std::string rtsp_ok = "RTSP/1.0 200 OK\r\nCSeq: " + boost::lexical_cast< std::string >( this->s_c_seq_ ) + "\r\nDate: " + time_string + "\r\n" +
			"Session: " + boost::lexical_cast< std::string >( this->m_session_ ) + "\r\nTransport: RTP/AVP;unicast;\r\n" +
			"client_port=" + boost::lexical_cast< std::string >( this->m_client_rtp_port_ ) + "-" + boost::lexical_cast< std::string >( this->m_client_rtcp_port_ ) +
			";server_port=" + boost::lexical_cast< std::string >( this->m_scirun_rtp_port_ ) + "-" + boost::lexical_cast< std::string >( this->m_scirun_rtcp_port_ ) + "\r\n\r\n";
		boost::asio::async_write( *this->m_socket_handle_, 
			boost::asio::buffer( rtsp_ok.c_str(), rtsp_ok.size() ), 
			boost::bind( &RTSPMediaCommunicator::handle_write, boost::dynamic_pointer_cast< RTSPMediaCommunicator >( shared_from_this() ), _1 ) );
	}

	void RTSPMediaCommunicator::send_error()
	{
		std::string rtsp_error = "RTSP/1.0 400 Bad Request\r\nCSeq: " + boost::lexical_cast< std::string >( this->s_c_seq_ ) + "\r\n";
		rtsp_error += "Session: " + boost::lexical_cast< std::string >( this->m_session_ ) + "\r\n\r\n";
		boost::asio::async_write( *this->m_socket_handle_, 
			boost::asio::buffer( rtsp_error.c_str(), rtsp_error.size() ), 
			boost::bind( &RTSPMediaCommunicator::handle_write, boost::dynamic_pointer_cast< RTSPMediaCommunicator >( shared_from_this() ), _1 ) );
	}

	void RTSPMediaCommunicator::send_teardown()
	{
		std::string teardown_string = "TEARDOWN rtsp://sci.edu RTSP/1.0\r\n";
		teardown_string += "CSeq: " + boost::lexical_cast< std::string >( s_c_seq_ ) + "\r\n";
		teardown_string += "Session: " + boost::lexical_cast< std::string >( this->m_session_ ) + "\r\n\r\n";
		boost::asio::async_write( *this->m_socket_handle_, 
			boost::asio::buffer( teardown_string.c_str(), teardown_string.size() ), 
			boost::bind( &RTSPMediaCommunicator::handle_write, boost::dynamic_pointer_cast< RTSPMediaCommunicator >( shared_from_this() ), _1 ) );
	}

	void RTSPMediaCommunicator::register_receive()
	{
		boost::asio::async_read_until(*this->m_socket_handle_, this->m_stream_data_,
			boost::regex( "\r\n\r\n" ),
			boost::bind( &RTSPMediaCommunicator::handle_receive, 
			boost::dynamic_pointer_cast< RTSPMediaCommunicator >( shared_from_this() ), _1, _2 ) );
	}

	void RTSPMediaCommunicator::handle_receive( const boost::system::error_code& error, size_t bytes_transferred )
	{
		if ( error )// if there is an error then do the callback (if any) and return without registering another function
		{
			this->m_ok_callback_.clear();
			call_communication_close_functor();
			stop();
			return;
		}
		const char* rtsp_buf = boost::asio::buffer_cast<const char*>( this->m_stream_data_.data() );
		this->m_stream_data_.consume( bytes_transferred );
		std::string RTSPmessage( rtsp_buf );
		std::vector< std::string > message_split;
		boost::algorithm::split( message_split, RTSPmessage, boost::is_any_of( " " ) );
		assert( message_split.size() > 4 );
		if ( message_split[0] == "SETUP" ) // the client sent us a setup message
		{
			parse_setup( RTSPmessage );
			//if ( !this->m_setup_callback_.empty() )
			//{
			//	this->m_setup_callback_( boost::dynamic_pointer_cast< RTSPMediaCommunicator >( shared_from_this() ), parse_setup( RTSPmessage ) );
			//}
		}
		else if ( message_split[0] == "TEARDOWN" ) // the server sent us a TEARDOWN
		{
			//if ( !this->m_teardown_callback_.empty() )
			//{
			//	this->m_teardown_callback_();
			//}
		}
		else if ( message_split[1] == "200" ) // the server sent us an OK
		{
			parse_ok( RTSPmessage );
			//if ( !this->m_ok_callback_.empty() )
			//{
			//	this->m_ok_callback_( parse_setup( RTSPmessage ) );
			//}
		}
		else if ( message_split[1] == "400" ) // the server sent us an error
		{
			//if ( !this->m_error_callback_.empty() )
			//{
			//	this->m_error_callback_();
			//}
		}
		register_receive();
	}

	bool RTSPMediaCommunicator::parse_setup( const std::string& rtsp_message )
	{
		std::string label_map_name, mesh_name, username = "";
//		unsigned short stage_number = 0;
		try
		{
			boost::regex reg( "SETUP |\\?|CSeq: |client_port=|-|\r\n" );
			boost::sregex_token_iterator iter( rtsp_message.begin(), rtsp_message.end(), reg, -1 );
			boost::sregex_token_iterator end_iter;
			std::advance( iter, 1 );
			this->m_label_map_name_ = *iter;
			std::advance( iter, 1 );
			this->m_mesh_name_ += *iter;
			std::advance( iter, 1 );
			this->m_username_ += *iter;
			std::advance( iter, 1 );
			this->m_stage_number_ = boost::lexical_cast< unsigned short >( *iter );
			std::advance( iter, 1 );
			this->m_communication_type_ = ( MediaCommunicationTypes::MediaCommunicationType )boost::lexical_cast< unsigned int >( *iter );
			std::advance( iter, 3 );
			this->s_c_seq_ = boost::lexical_cast< int >( *iter );
			std::advance( iter, 2 );
			this->m_client_rtp_port_ = boost::lexical_cast< unsigned short >( *iter );
			std::advance( iter, 1 );
			this->m_client_rtcp_port_ = boost::lexical_cast< unsigned short >( *iter );
			this->m_setup_callback_( this->m_label_map_name_, this->m_mesh_name_, this->m_stage_number_ );
		}
		catch ( ... )
		{
			return false;
		}
		return true;
	}
	bool RTSPMediaCommunicator::parse_ok( const std::string& rtsp_message )
	{
		boost::regex reg( "client_port=|-|;server_port=|-|\r\n" );
		boost::sregex_token_iterator iter( rtsp_message.begin(), rtsp_message.end(), reg, -1 );
		boost::sregex_token_iterator end_iter;
		std::vector< std::string > rtsp_message_split;
		for ( ; iter != end_iter; ++iter)
		{
			rtsp_message_split.push_back( *iter );
		}
		try
		{
			this->m_client_rtp_port_ = boost::lexical_cast< unsigned short >( rtsp_message_split[6] );
			this->m_client_rtcp_port_ = boost::lexical_cast< unsigned short >( rtsp_message_split[7] );
			this->m_scirun_rtp_port_ = boost::lexical_cast< unsigned short >( rtsp_message_split[8] );
			this->m_scirun_rtcp_port_ = boost::lexical_cast< unsigned short >( rtsp_message_split[9] );
			this->m_ok_callback_( this->m_client_rtp_port_, this->m_scirun_rtp_port_ );
		}
		catch ( boost::bad_lexical_cast err )
		{
			// maybe you could use a error callback so you can notify the client
		}
		return true;
	}

	void RTSPMediaCommunicator::set_setup_callback( boost::function< void ( std::string, std::string, unsigned short ) > setup_callback_ )
	{
		this->m_setup_callback_ = setup_callback_;
	}

}
//
//	int RtspMediaCommunication::s_cseq_( 1 );
//
//	RtspMediaCommunication::RtspMediaCommunication( const std::string& label_map_name, 
//		const std::string& mesh_name, unsigned short stage_number, const std::string& ip_address ) :
//			TCPMediaCommunication( ip_address, 554 ), m_label_map_name_( label_map_name ),
//			m_mesh_name_( mesh_name ), m_stage_number_( stage_number )
//	{
//		this->m_streambuf_.prepare( this->m_recieve_buffer_size_ );
//		MediaCommunicationBase::get_available_port( 35123, 65000, this->m_local_port_, 
//			MediaCommunicationTypes::RTP_E );
//		if ( this->m_is_client_ ) 
//		{
//			this->m_rtsp_session_id_ = 1;
//		}
//	}
//
//	RtspMediaCommunication::RtspMediaCommunication( boost::function< void ( std::string, unsigned short, 
//		std::string, unsigned short& ) > func ) :	TCPMediaCommunication( "", 554 ), 
//		m_setup_callback_( func )
//	{
//		this->m_streambuf_.prepare( this->m_recieve_buffer_size_ );
//		MediaCommunicationBase::get_available_port( 35123, 65000, this->m_local_port_, 
//			MediaCommunicationTypes::RTP_E );
//		this->m_rtsp_session_id_ = 1;
//	}
//
//	RtspMediaCommunication::~RtspMediaCommunication()
//	{
//		boost::mutex::scoped_lock destructor_lock( this->m_destructor_mutex_ );
//		this->m_destructing_ = true;
//	}
//
//	bool RtspMediaCommunication::pinhole_firewall_get_port( unsigned short& client_rtp_port, 
//		unsigned short& client_rtcp_port, unsigned short& server_rtp_port, unsigned short& server_rtcp_port )
//	{
//		client_rtp_port = this->m_local_port_;
//		client_rtcp_port = this->m_local_port_ + 1;
//		while ( !is_connected() )
//		{
//			boost::this_thread::sleep( boost::posix_time::milliseconds( 500 ) );
//		}
//		std::string setup_message = "SETUP " + m_label_map_name_ + "?" + m_mesh_name_ + "?" +
//			boost::lexical_cast< std::string >( this->m_stage_number_ ) + "? RTSP/1.0\r\nCSeq: " + 
//			boost::lexical_cast< std::string >( s_cseq_ ) + "\r\nTransport: RTP/AVP;unicast;client_port=" +
//			boost::lexical_cast< std::string >( this->m_local_port_ ) + "-" +
//			boost::lexical_cast< std::string >( this->m_local_port_ + 1 ) + "\r\n\r\n";
//		this->send_raw_bytes( ( unsigned char* )setup_message.c_str(), setup_message.size() );
//		boost::mutex::scoped_lock response_mutex( this->m_response_mutex_ );
//		this->m_response_condition_.timed_wait( response_mutex, boost::posix_time::seconds( 500 ) );
//		if ( this->m_response_string_.empty() ) { return false; } // if there was a timeout
//		try
//		{
//			boost::regex reg( "Session: |server_port=|-|\r\n" );
//			boost::sregex_token_iterator iter( this->m_response_string_.begin(), 
//				this->m_response_string_.end(), reg, -1 );
//			{
//				boost::mutex::scoped_lock session_lock( this->m_rtsp_session_id_mutex_ );
//				std::advance( iter, 5 );
//				this->m_rtsp_session_id_ = boost::lexical_cast< int >( *iter );
//			}
//			{
//				std::advance( iter, 4 );
//				client_rtp_port = boost::lexical_cast< unsigned short >( *iter );
//				std::advance( iter, 1 );
//				client_rtcp_port = boost::lexical_cast< unsigned short >( *iter );
//			}
//		}
//		catch ( ... )
//		{
//			return false;
//		}
//		return true;
//	}
//
//	void RtspMediaCommunication::register_callbacks()
//	{
//		this->m_socket->async_receive( boost::asio::buffer( this->m_receive_buffer_, this->m_recieve_buffer_size_ ),
//			boost::bind( &RtspMediaCommunication::handle_receive_rtsp, _1, _2, this ) );
//		//boost::asio::async_read_until( *this->m_socket, 
//		//	this->m_streambuf_,	"\r\n\r\n",
//		//	boost::bind( &RtspMediaCommunication::handle_receive_rtsp, _1, _2, this ) );
//	}
//
//	void RtspMediaCommunication::handle_receive_rtsp( const boost::system::error_code& e_code, 
//		size_t bytes_transferred, RtspMediaCommunication* communicator )
//	{
//		boost::mutex::scoped_lock destructor_lock( communicator->m_destructor_mutex_ );
//		if ( e_code || !is_communicator_valid( communicator ) ) // if there is an error then stop the socket
//		{
//			std::cout << e_code.message() << std::endl;
//			if ( is_communicator_valid( communicator ) )
//			{
//				communicator->stop();
//			}
//			return;
//		}
//		boost::asio::streambuf::const_buffers_type bufs = communicator->m_streambuf_.data();
//		std::string rtsp_message( ( char* )communicator->m_receive_buffer_, bytes_transferred );
//		if ( communicator->m_is_client_ ) // each time we get a message and we are a client then trigger the condition variable
//		{
//			boost::mutex::scoped_lock response_mutex( communicator->m_response_mutex_ );
//			communicator->m_response_string_ = rtsp_message;
//			communicator->m_response_condition_.notify_all();
//		}
//		else
//		{
//			std::vector< std::string > message_split;
//			boost::algorithm::split( message_split, rtsp_message, boost::is_any_of( " " ) );
//			if ( message_split.size() > 0 && strcmp( message_split[0].c_str(), "SETUP" ) == 0 )
//			{
//				communicator->setup_rtsp( rtsp_message );
//			}
//		}
//		communicator->register_callbacks();
//	}
//
//	void RtspMediaCommunication::setup_rtsp( const std::string& rtsp_message )
//	{
//		bool rtsp_error = false;
//		unsigned short client_rtp_port = 0;
//		unsigned short client_rtcp_port = 0;
//		unsigned short server_rtp_port = 0;
//		unsigned short server_rtcp_port = 0;
//		int c_seq = 0;
//		std::string full_project_name = "";
//		unsigned short stage_number = 0;
//		std::string user_name = "";
//		try
//		{
//			boost::regex reg( "SETUP |\\?|CSeq: |client_port=|-|\r\n" );
//			boost::sregex_token_iterator iter( rtsp_message.begin(), rtsp_message.end(), reg, -1 );
//			boost::sregex_token_iterator end_iter;
//			std::advance( iter, 1 );
//			full_project_name = *iter;
//			std::advance( iter, 1 );
//			full_project_name += "/" + *iter;
//			std::advance( iter, 1 );
//			stage_number = boost::lexical_cast< unsigned short >( *iter );
//			std::advance( iter, 3 );
//			c_seq = boost::lexical_cast< int >( *iter );
//			std::advance( iter, 2 );
//			client_rtp_port = boost::lexical_cast< unsigned short >( *iter );
//			std::advance( iter, 1 );
//			client_rtcp_port = boost::lexical_cast< unsigned short >( *iter );
//			this->m_setup_callback_( full_project_name, stage_number, get_ip_address(), server_rtp_port );
//			server_rtcp_port = server_rtp_port + 1;
//			rtsp_error = server_rtp_port > 0;
//		}
//		catch ( ... )
//		{
//			rtsp_error = true;
//		}
//		if ( !rtsp_error )
//		{
//			time_t t;
//			time( &t );
//			char* time_char = asctime( gmtime( &t ) );
//			std::string time_string( time_char );
//			server_rtcp_port = server_rtp_port + 1;
//			std::string rtp_ok = "RTSP/1.0 200 OK\r\nCSeq: " + boost::lexical_cast< std::string >( c_seq ) + "\r\nDate: " + time_string + "\r\n" +
//				"Session: " + boost::lexical_cast< std::string >( this->m_rtsp_session_id_ ) + "\r\nTransport: RTP/AVP;unicast;\r\n" +
//				"client_port=" + boost::lexical_cast< std::string >( client_rtp_port ) + "-" + boost::lexical_cast< std::string >( client_rtcp_port ) +
//				";server_port=" + boost::lexical_cast< std::string >( server_rtp_port ) + "-" + boost::lexical_cast< std::string >( server_rtcp_port ) + "\r\n\r\n";
//			send_raw_bytes( ( unsigned char* )rtp_ok.c_str(), rtp_ok.size() );
//		}
//	}
//
//} // end namespace 
//
