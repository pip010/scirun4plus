//#include <iostream>
//#include <fstream>
//
//#include <boost/asio.hpp>
//
//#ifdef WIN32
//#include <crtdbg.h> //for _ASSERT
//#else
//#include <cassert>
//#define _ASSERT assert
//#endif
////JRTP includes
//
//
//#include "rtppacket.h"
//#include "rtppacketbuilder.h"
//#include "rtpudpv4transmitter.h"
//#include "rtpipv4address.h"
//#include "rtpsessionparams.h"
//#include "rtperrors.h"
//#include "rtpsourcedata.h"
//#include "rtppacketbuilder.h"
//
//#include "MediaCommunicationBase.h"
//#include "RTPMediaCommunication.h"
//
//#include <Message.h>
//#include <MessageParser.h>
//
////#define RTP_SUPPORT_PROBATION 1
//
//namespace BioMesh3d
//{
//
//	unsigned short RTPMediaCommunication::s_client_port_ = 9000;
//
//	RTPMediaCommunication::RTPMediaCommunication( const std::string& ip_address, unsigned short port,
//		boost::function< void ( MessageHandle ) > handle_message_functor ) : 
//		MediaCommunicationBase( ip_address, port, handle_message_functor ),
//		m_peer_ip_address_( ip_address ), m_port_( port ),
//		m_last_timestamp_displayed_( 0 )
//	{
//
//#ifdef WIN32
//		WSADATA wsaData;
//		WORD wVersionRequested = MAKEWORD( 2, 2 );
//		WSAStartup(wVersionRequested, &wsaData); //WSACleanup called in OnDestroy
//#endif
//
//		for ( int i = 0; i < 3; ++i )
//		{
//			BinaryBuffer* buf = new BinaryBuffer();
//			this->m_binary_buffers_.insert( this->m_binary_buffers_.begin(), buf );
//		}
//		m_mem_mgr_ = new MediaMemoryManager();
//		SetMemoryManager( m_mem_mgr_ );
//
//		this->m_local_ip_address_ = get_local_ip_address();
//
//	}
//
//	RTPMediaCommunication::~RTPMediaCommunication()
//	{
//		{
//			boost::mutex::scoped_lock lock( this->m_send_ping_mutex_ );
//			this->m_send_ping_ = false;
//		}
//		Destroy();
//		delete this->trans_params_;
//		delete this->sess_params_;
//
//#ifdef WIN32
//		WSACleanup();
//#endif
//	}
//
//	std::string RTPMediaCommunication::get_local_ip_address()
//	{
//		boost::asio::io_service io_service;
//		boost::asio::ip::tcp::resolver resolver( io_service );
//		boost::asio::ip::tcp::resolver::query query( boost::asio::ip::host_name(), "" );
//		boost::asio::ip::tcp::resolver::iterator iter = resolver.resolve( query );
//		boost::asio::ip::tcp::resolver::iterator end;
//		std::string local_ip_address = "";
//		while ( iter != end )
//		{
//			boost::asio::ip::tcp::endpoint ep = *iter;
//			local_ip_address = ep.address().to_string();
//			++iter;
//		}
//		return local_ip_address;
//	}
//
//	void RTPMediaCommunication::run()
//	{
//		this->m_thread_ = new boost::thread( boost::bind( &RTPMediaCommunication::init_service, this ) );
//	}
//
//	void RTPMediaCommunication::init_service()
//	{
//		bool is_client = this->m_peer_ip_address_ != "";
//
//		this->sess_params_ = new RTPSessionParams();
//		this->trans_params_= new RTPUDPv4TransmissionParams();
//
//		this->sess_params_->SetOwnTimestampUnit( 1.0 / 100.0 ); //30 video frames per second
//		this->sess_params_->SetUsePollThread( 1 ); //background thread to call virtual callbacks - set by default, but just to be sure
//		this->sess_params_->SetMaximumPacketSize( MAX_PACKET_SIZE );
//
//		bool socket_bind_error = true;
//		int status;
//		while ( socket_bind_error )
//		{
//			unsigned short usage_port = ( is_client ) ? this->s_client_port_ : this->m_port_;
//			this->trans_params_->SetPortbase( usage_port ); // this is the port this computer uses
//		
//			status = this->Create( *this->sess_params_, this->trans_params_ );
//			if ( report_error( status ) )
//			{
//				std::cerr << RTPGetErrorString(status) << std::endl;
//				std::cerr << "error in creating session...\n";
//				this->s_client_port_ += 2;
//			}
//			else
//			{
//				socket_bind_error = false;
//			}
//		}
//
//		port_writer.open( "/home/collab/dillonl/temp/network_info.txt", std::ios::app );
//		port_writer << "--- Start Log ---" << std::endl;
//
//		SetDefaultPayloadType(96);
//		SetDefaultMark(false);
//		SetMaximumPacketSize( MAX_PACKET_SIZE );
//		SetDefaultTimestampIncrement( 0 );
//
//		this->m_send_ping_ = is_client;
//		if ( is_client ) //  client
//		{
//			port_writer << "I'm a ThinClient" << std::endl;
//			//JOIN THE MULTICAST
//			unsigned long intIP = inet_addr( this->m_peer_ip_address_.c_str() );
//			_ASSERT(intIP != INADDR_NONE);
//			intIP = ntohl(intIP); //put in host byte order
//
//			RTPIPv4Address server_addr( intIP, this->m_port_ ); //the port opened on the server
//
//			status = AddDestination( server_addr );
//			if (status < 0)
//			{
//				std::cerr << RTPGetErrorString(status) << std::endl;
//				exit(-1);
//			}
//			//port_writer << "SCIRun's ip address is: " << this->m_peer_ip_address_ << "SCIRun is listening on port: " << boost::lexical_cast< std::string >( this->m_port_ ) << std::endl;
//		}
//		else
//		{
//			port_writer << "I'm a SCIRun" << std::endl;
//		}
//		unsigned short tmp_port = ( is_client ) ? this->s_client_port_ : this->m_port_;
//		port_writer << "My ip address: " << this->m_local_ip_address_ << " I'm listening on port: " << tmp_port << std::endl;
//		port_writer.close();
//	}
//
//	int RTPMediaCommunication::report_error( int errCode )
//	{
//		int isErr = ( errCode < 0 );
//		if ( isErr )
//		{
//			std::string stdErrStr = RTPGetErrorString(errCode);
//			std::cerr <<  errCode << "-->" << stdErrStr.c_str() << std::endl;
//		}
//		return isErr;
//	}
//
//	void RTPMediaCommunication::OnNewSource( RTPSourceData *dat )
//	{
//		if (dat->IsOwnSSRC())
//		{
//			return;
//		}
//
//		uint32_t ip;
//		uint16_t port;
//
//		if (dat->GetRTPDataAddress() != 0)
//		{
//			const RTPIPv4Address *addr = (const RTPIPv4Address *)(dat->GetRTPDataAddress());
//			ip = addr->GetIP();
//			port = addr->GetPort();
//		}
//		else if (dat->GetRTCPDataAddress() != 0)
//		{
//			const RTPIPv4Address *addr = (const RTPIPv4Address *)(dat->GetRTCPDataAddress());
//			ip = addr->GetIP();
//			port = addr->GetPort() - 1;
//		}
//		else
//		{
//			return;
//		}
//
//		RTPIPv4Address dest(ip,port);
//		AddDestination(dest);
//
//		struct in_addr inaddr;
//		inaddr.s_addr = htonl(ip);
//
//		port_writer.open( "/home/collab/dillonl/temp/network_info.txt", std::ios::app );
//		port_writer << "Destination ip address is: " << std::string(inet_ntoa(inaddr)) << " Destination is listening on port: " << port << std::endl;
//		port_writer << "--- End Log ---" << std::endl;
//		port_writer.close();
//
//		std::cout << "Adding destination " << std::string(inet_ntoa(inaddr)) << ":" << port << std::endl;
//		this->set_is_connected( true );
//		
//		//add new receiver buffer
//		m_receivers_buffer_[ dat->GetSSRC() ] =
//			boost::shared_ptr< ReceiverBuffer > ( new ReceiverBuffer() );		
//	}
//
//	void RTPMediaCommunication::OnBYEPacket( RTPSourceData *srcdat )
//	{
//		this->set_is_connected( false );
//
//		if (srcdat->IsOwnSSRC())
//		{
//			return;
//		}
//
//		uint32_t ip;
//		uint16_t port;
//
//		if (srcdat->GetRTPDataAddress() != 0)
//		{
//			const RTPIPv4Address *addr = (const RTPIPv4Address *)(srcdat->GetRTPDataAddress());
//			ip = addr->GetIP();
//			port = addr->GetPort();
//		}
//		else if (srcdat->GetRTCPDataAddress() != 0)
//		{
//			const RTPIPv4Address *addr = (const RTPIPv4Address *)(srcdat->GetRTCPDataAddress());
//			ip = addr->GetIP();
//			port = addr->GetPort()-1;
//		}
//		else
//		{
//			return;
//		}
//
//		RTPIPv4Address dest(ip,port);
//		DeleteDestination(dest);
//
//		struct in_addr inaddr;
//		inaddr.s_addr = htonl(ip);
//		std::cout << "Deleting destination " << std::string(inet_ntoa(inaddr)) << ":" << port << std::endl;
//	}
//
//
//	void RTPMediaCommunication::OnRemoveSource(RTPSourceData *dat)
//	{
//		if (dat->IsOwnSSRC())
//		{
//			return;
//		}
//		if (dat->ReceivedBYE())
//		{
//			return;
//		}
//
//		uint32_t ip;
//		uint16_t port;
//
//		if (dat->GetRTPDataAddress() != 0)
//		{
//			const RTPIPv4Address *addr = (const RTPIPv4Address *)(dat->GetRTPDataAddress());
//			ip = addr->GetIP();
//			port = addr->GetPort();
//		}
//		else if (dat->GetRTCPDataAddress() != 0)
//		{
//			const RTPIPv4Address *addr = (const RTPIPv4Address *)(dat->GetRTCPDataAddress());
//			ip = addr->GetIP();
//			port = addr->GetPort()-1;
//		}
//		else
//		{
//			return;
//		}
//
//		RTPIPv4Address dest(ip,port);
//		DeleteDestination(dest);
//
//		struct in_addr inaddr;
//		inaddr.s_addr = htonl(ip);
//		std::cout << "Deleting destination " << std::string(inet_ntoa(inaddr)) << ":" << port << std::endl;
//
//		//delete new receiver buffer
//		m_receivers_buffer_.erase( dat->GetSSRC() ); 
//	}
//
//	void RTPMediaCommunication::send_message( MessageHandle msg )
//	{
//		if ( this->is_connected() )
//		{
//			int max_buf_size = MaxNetworkBufferSize::MAX_BUFFER_SIZE_E;
//			boost::mutex::scoped_lock lock( this->m_mutex_ );
//
//			std::string serialized_params = msg->serialize_parameters();
//			size_t send_size =  MESSAGE_HEADER_SIZE + serialized_params.size() + 1 + msg->get_binary_size();
//			unsigned short fragments = ( ( send_size + MaxNetworkBufferSize::MAX_BUFFER_SIZE_E - 1 ) / MaxNetworkBufferSize::MAX_BUFFER_SIZE_E );
//			MessageParser::set_message_size( this->m_send_buffer_, send_size ); // the message size is the header size + the param string + binary size + 1
//			MessageParser::set_message_type( this->m_send_buffer_, msg->get_message_type() );
//			memcpy( this->m_send_buffer_ + MESSAGE_HEADER_SIZE, serialized_params.c_str(), serialized_params.size() );
//			this->m_send_buffer_[MESSAGE_HEADER_SIZE + serialized_params.size()] = 0;
//			memcpy( reinterpret_cast< char* >( this->m_send_buffer_ + MESSAGE_HEADER_SIZE + serialized_params.size() + 1 ), 
//				msg->get_binary_buf(), msg->get_binary_size() );
//
//			size_t bytes_transferred = 0;
//			size_t tmp_send_size;
//			bool marker;
//			uint8_t sequence_number = 1;
//			unsigned char send_size_hdr[sizeof( size_t )];
//			*reinterpret_cast< size_t* >( send_size_hdr ) = send_size;
//			while ( bytes_transferred < send_size )			
//			{
//				marker = ( send_size - bytes_transferred < MaxNetworkBufferSize::MAX_BUFFER_SIZE_E );
//				tmp_send_size = ( !marker ) ? MaxNetworkBufferSize::MAX_BUFFER_SIZE_E : 
//					send_size - bytes_transferred;
//				SendPacketEx( this->m_send_buffer_ + bytes_transferred, tmp_send_size, sequence_number, marker, 
//					0, fragments, send_size_hdr, sizeof( size_t ) );
//				++sequence_number;
//				bytes_transferred += tmp_send_size;
//
//			}
//
//			//this->IncrementTimestamp( 1 );
//			this->IncrementTimestamp( 100 );
//		}
//		else
//		{
//			set_is_connected( true );
//		}
//	}
//
//	/*BinaryBuffer* RTPMediaCommunication::reconstruct_buffer( unsigned char* data, 
//		size_t data_len, size_t total_size, uint32_t timestamp, uint8_t fragment_number )
//
//	{
//		int oldest_buffer_index = -1;
//		size_t oldest_timestamp = timestamp;
//		BinaryBuffer* tmp_buffer = NULL;
//		for ( unsigned int i = 0; i < this->m_binary_buffers_.size(); ++i )
//		{
//			if ( this->m_binary_buffers_[i]->get_timestamp() == timestamp )
//			{
//				tmp_buffer = this->m_binary_buffers_[i];
//				break;
//			}
//			else if ( this->m_binary_buffers_[i]->get_timestamp() < oldest_timestamp )
//			{
//				tmp_buffer = this->m_binary_buffers_[i];
//				oldest_timestamp = tmp_buffer->get_timestamp();
//			}
//		}
//		if ( tmp_buffer != NULL )
//		{
//			tmp_buffer->set_timestamp_and_reset_info( timestamp );
//			tmp_buffer->increment_fragments();
//			tmp_buffer->set_binary_size( total_size );
//			if ( tmp_buffer->get_buffer_size() < total_size )
//			{
//				tmp_buffer->resize_buffer( total_size * 2 );
//			}
//			size_t offset = ( fragment_number - 1 ) * MaxNetworkBufferSize::MAX_BUFFER_SIZE_E;
//			memcpy( reinterpret_cast< char* >( tmp_buffer->get_buffer() + offset ), data, data_len );
//		}
//		return tmp_buffer;
//	}*/
//
//	void RTPMediaCommunication::ProcessRTPPacket( const RTPSourceData &srcdat,const RTPPacket &pack )
//	{	
//		size_t  bytes_transferred  = pack.GetPayloadLength();
//		unsigned char* packet_data = pack.GetPayloadData();
//		unsigned int timestamp = pack.GetTimestamp();
//		unsigned short sequence_number = pack.GetSequenceNumber();
//		unsigned short fragments = pack.GetExtensionID(); // I hijacked this to store the number of fragments
//		unsigned short fragment_number = pack.GetPayloadType();
//		//unsigned short fragment_number = pack.GetSequenceNumber();
//		bool marker = pack.HasMarker();
//		unsigned int send_size_hdr = *reinterpret_cast< unsigned int* >( pack.GetExtensionData() );
//		unsigned int ssrc = pack.GetSSRC();
//
//		//std::cout << "timestamp: " << timestamp << " packet size: " << bytes_transferred << " fragments: " << fragments << " fragment number: " << fragment_number << std::endl;
//
//		std::map<unsigned int, unsigned int>::iterator pos;
//		pos = m_receivers_timestamp_.find( ssrc );
//		if ( pos!= m_receivers_timestamp_.end() )
//		{
//			if (  timestamp < pos->second )
//			{
//				//std::cout << "I arrived, but late!" << timestamp << std::endl;
//				return;
//			}
//			else
//			{
//				pos->second = timestamp;
//			}
//		}
//		else
//		{
//			//m_receivers_timestamp_[pack.GetSSRC()] = timestamp;
//			m_receivers_timestamp_.insert( std::pair<unsigned int, unsigned int> (ssrc, timestamp) );
//
//		}
//
//		//std::cout << "second time --> timestamp: " << timestamp << " packet size: " << bytes_transferred << " fragments: " << fragments << " fragment number: " << fragment_number << std::endl;
//
//		//if ( timestamp < this->m_last_timestamp_displayed_ ) // we don't care about this frame anymore
//		//{
//		//	return;
//		//}
//
//		//BinaryBuffer* buf = reconstruct_buffer( packet_data, bytes_transferred, 
//		//	send_size_hdr, timestamp, fragment_number );
//
//		std::map<unsigned int, boost::shared_ptr< ReceiverBuffer > >::iterator receiver_pos = 
//			(m_receivers_buffer_.find( ssrc ));
//
//		boost::shared_ptr<ReceiverBuffer> receiver_buffer;
//
//		if ( receiver_pos != m_receivers_buffer_.end() ) //find right buffer
//		{
//			 receiver_buffer =  receiver_pos->second;
//		}
//		else
//		{
//			return;
//		}
//		
//
//		boost::shared_ptr< BinaryBuffer > buf = 
//			receiver_buffer->reconstruct_buffer( packet_data, bytes_transferred, 
//				send_size_hdr, timestamp, fragment_number );
//
//
//		if ( buf == NULL || buf->fragments_collected() < fragments )
//		{
//			return;
//		}
//
//		//std::cout << "Receive a complete frame has " << fragments << " fragments! \n\n"; 
//		//this->m_message_manager_->preprocess_buffer( &buf->get_buffer()[0] );
//		assert( false );
//
//
//		//this->m_last_timestamp_displayed_ = timestamp;
//	}
//
//	void RTPMediaCommunication::OnPollThreadStep()
//	{
//		BeginDataAccess();
//
//		// check incoming packets
//		if (GotoFirstSourceWithData())
//		{
//			do
//			{
//				RTPPacket *pack;
//				RTPSourceData *srcdat;
//
//				srcdat = GetCurrentSourceInfo();
//
//				while ((pack = GetNextPacket()) != NULL)
//				{
//					set_is_connected( true );
//					ProcessRTPPacket(*srcdat,*pack);
//					DeletePacket(pack);
//				}
//			} while (GotoNextSourceWithData());
//		}
//
//		EndDataAccess();
//	}
//
//	void RTPMediaCommunication::OnRTPPacket( RTPPacket *pack, const RTPTime &receivetime, const RTPAddress *senderaddress )
//	{
//		size_t  bytes_transferred  = pack->GetPayloadLength();
//		unsigned char* packet_data = pack->GetPayloadData();
//		unsigned int timestamp = pack->GetTimestamp();
//		unsigned short sequence_number = pack->GetSequenceNumber();
//		unsigned short fragments = pack->GetExtensionID(); // I hijacked this to store the number of fragments
//		unsigned short fragment_number = pack->GetPayloadType();
//		bool marker = pack->HasMarker();
//		unsigned int send_size_hdr = *reinterpret_cast< unsigned int* >( pack->GetExtensionData() );
//		unsigned int ssrc = pack->GetSSRC();
//
//		//std::cout << "On Arrive timestamp: " << timestamp << " packet size: " << bytes_transferred << " fragments: " << fragments << " fragment #: " << fragment_number << std::endl;
//
//	}
//
//	//void RTPMediaCommunication::OnRTPPacket( RTPPacket *pack, const RTPTime &receivetime, const RTPAddress *senderaddress )
//	//{
//	//	size_t  bytes_transferred  = pack->GetPayloadLength();
//	//	unsigned char* packet_data = pack->GetPayloadData();
//	//	unsigned int timestamp = pack->GetTimestamp();
//	//	unsigned short sequence_number = pack->GetSequenceNumber();
//	//	unsigned short fragments = pack->GetExtensionID(); // I hijacked this to store the number of fragments
//	//	unsigned short fragment_number = pack->GetPayloadType();
//	//	bool marker = pack->HasMarker();
//	//	unsigned int send_size_hdr = *reinterpret_cast< unsigned int* >( pack->GetExtensionData() );
//	//	unsigned int ssrc = pack->GetSSRC();
//
//	//	std::cout << "timestamp: " << timestamp << " packet size: " << bytes_transferred << " fragments: " << fragments << " fragment number: " << fragment_number << std::endl;
//
//	//	std::map<unsigned int, unsigned int>::iterator pos;
//	//	pos = m_receivers_timestamp_.find( ssrc );
//	//	if ( pos!= m_receivers_timestamp_.end() )
//	//	{
//	//		if (  timestamp < pos->second )
//	//		{
//	//			return;
//	//		}
//	//		else
//	//		{
//	//			pos->second = timestamp;
//	//		}
//	//	}
//	//	else
//	//	{
//	//		//m_receivers_timestamp_[pack.GetSSRC()] = timestamp;
//	//		m_receivers_timestamp_.insert( std::pair<unsigned int, unsigned int> (ssrc, timestamp) );
//
//	//	}
//
//	//	//if ( timestamp < this->m_last_timestamp_displayed_ ) // we don't care about this frame anymore
//	//	//{
//	//	//	return;
//	//	//}
//
//	//	//BinaryBuffer* buf = reconstruct_buffer( packet_data, bytes_transferred, 
//	//	//	send_size_hdr, timestamp, fragment_number );
//
//	//	std::map<unsigned int, boost::shared_ptr< ReceiverBuffer > >::iterator receiver_pos = 
//	//		(m_receivers_buffer_.find( ssrc ));
//
//	//	boost::shared_ptr<ReceiverBuffer> receiver_buffer;
//
//	//	if ( receiver_pos != m_receivers_buffer_.end() ) //find right buffer
//	//	{
//	//		receiver_buffer =  receiver_pos->second;
//	//	}
//	//	else
//	//	{
//	//		return;
//	//	}
//
//
//	//	boost::shared_ptr< BinaryBuffer > buf = 
//	//		receiver_buffer->reconstruct_buffer( packet_data, bytes_transferred, 
//	//		send_size_hdr, timestamp, fragment_number );
//
//
//	//	if ( buf == NULL || buf->fragments_collected() < fragments )
//	//	{
//	//		return;
//	//	}
//	//	this->m_message_manager_->preprocess_buffer( &buf->get_buffer()[0] );
//	//}
//
//} //end of namespace
