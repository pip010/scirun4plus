//#ifndef RTP_MEDIA_H
//#define RTP_MEDIA_H
//
//#include <map>
//#include <iostream>
//#include <fstream>
//
//#include "rtpsession.h"
//
//#include <boost/bind.hpp>
//#include <boost/cast.hpp>
//#include <boost/function.hpp>
//#include <boost/signals2.hpp>
//#include <boost/thread.hpp>
//
//#include "MediaCommunicationBase.h"
//
//#ifndef WIN32
//#include <netinet/in.h>
//#include <arpa/inet.h>
//#else
//#include <winsock2.h>
//#endif // WIN32
//
//#define MAX_PACKET_SIZE	((1024 * 64) - 1)
//#define MAX_BINARY_SIZE 6000000
//
////class RTPSession;
//class RTPSessionParams;
//class RTPUDPv4TransmissionParams;
//
//namespace BioMesh3d
//{
//
//	class BinaryBuffer
//	{
//	public:
//		BinaryBuffer() : m_fragments_( 0 ), m_timstamp_( 0 ), m_binary_size_( 0 ){}
//		~BinaryBuffer() {}
//
//		uint8_t fragments_collected() { return this->m_fragments_; }
//		void increment_fragments() { ++m_fragments_; }
//
//		uint32_t get_timestamp() { return this->m_timstamp_; }
//		void set_timestamp_and_reset_info( uint32_t timestamp ) 
//		{ 
//			if ( timestamp != this->m_timstamp_ )
//			{
//				this->m_fragments_ = 0;
//				this->m_binary_size_ = 0;
//				this->m_timstamp_ = timestamp; 
//			}
//		}
//
//		size_t get_binary_size() { return this->m_binary_size_; }
//		void set_binary_size( size_t size ) { this->m_binary_size_ = size; }
//
//		unsigned char* get_buffer() { return &this->m_buffer_[0]; }
//		void resize_buffer( size_t size ) 
//		{
//			this->m_buffer_.resize( size ); 
//		}
//		size_t get_buffer_size() { return this->m_buffer_.size(); }
//	private:
//		uint8_t m_fragments_;
//		uint32_t m_timstamp_;
//		size_t m_binary_size_;
//		std::vector< unsigned char > m_buffer_;
//	};
//
//	class ReceiverBuffer
//	{
//	public:
//		ReceiverBuffer() 
//		{
//			for ( unsigned int i = 0; i < 3; ++i )
//			{
//				boost::shared_ptr< BinaryBuffer > buf =
//					boost::shared_ptr< BinaryBuffer > ( new BinaryBuffer() );
//				this->m_receiver_buffers_.push_back( buf );
//			}
//		}
//
//		~ReceiverBuffer()
//		{
//			this->m_receiver_buffers_.clear();
//		}
//
//		boost::shared_ptr< BinaryBuffer > reconstruct_buffer( unsigned char* data, 
//			size_t data_len, size_t total_size, uint32_t timestamp, uint8_t fragment_number )
//
//		{
//			int oldest_buffer_index = -1;
//			size_t oldest_timestamp = timestamp;
//			boost::shared_ptr< BinaryBuffer > tmp_buffer;
//			for ( unsigned int i = 0; i < this->m_receiver_buffers_.size(); ++i )
//			{
//				if ( this->m_receiver_buffers_[i]->get_timestamp() == timestamp )
//				{
//					tmp_buffer = this->m_receiver_buffers_[i];
//					break;
//				}
//				else if ( this->m_receiver_buffers_[i]->get_timestamp() < oldest_timestamp )
//				{
//					tmp_buffer = this->m_receiver_buffers_[i];
//					oldest_timestamp = tmp_buffer->get_timestamp();
//				}
//			}
//			if ( tmp_buffer != NULL )
//			{
//				tmp_buffer->set_timestamp_and_reset_info( timestamp );
//				tmp_buffer->increment_fragments();
//				tmp_buffer->set_binary_size( total_size );
//				if ( tmp_buffer->get_buffer_size() < total_size )
//				{
//					tmp_buffer->resize_buffer( total_size * 2 );
//				}
//				size_t offset = ( fragment_number - 1 ) * MaxNetworkBufferSize::MAX_BUFFER_SIZE_E;
//				memcpy( reinterpret_cast< char* >( tmp_buffer->get_buffer() + offset ), data, data_len );
//			}
//			return tmp_buffer;
//		}
//
//	private:
//		std::vector< boost::shared_ptr< BinaryBuffer > > m_receiver_buffers_;
//
//	};
//
//	class MediaMemoryManager : public RTPMemoryManager
//	{
//	public:
//		MediaMemoryManager()
//		{
//
//		}
//		~MediaMemoryManager()
//		{
//
//		}
//
//		virtual void* AllocateBuffer( size_t numbytes, int memtype )
//		{
//			//return operator new( numbytes );
//			return operator new( 3000000 );
//		}
//
//		virtual void FreeBuffer( void *buffer )
//		{
//			delete buffer;
//			buffer = NULL;
//		}
//
//	private:
//	};
//
//	class RTPMediaCommunication : public RTPSession, public MediaCommunicationBase
//	{
//	public:
//		RTPMediaCommunication( const std::string& peer_ip_address, unsigned short port, 
//			boost::function< void ( MessageHandle ) > handle_message_functor );
//		~RTPMediaCommunication();
//
//		//void set_receive_callback( boost::function< void ( unsigned char* media_data, size_t bytes_transferred ) > func );
//
//		void run();
//		void send_message( MessageHandle msg );
//
//		void stop() {} //replaced by stop_service
//		void stop_service();
//
//		//virtual void OnRTPPacket( RTPPacket *pack, const RTPTime &receivetime, const RTPAddress *senderaddress );
//
//
//		void OnNewSource( RTPSourceData *dat );
//
//		void set_peer_port( unsigned short port ) { s_client_port_ = port; }
//
//		unsigned short get_client_port() const { return s_client_port_; }
//
//		boost::signals2::signal< void( ) > buffer_ready_sig_;
//		void send_buffer_ready_signal( ) { buffer_ready_sig_() ; }
//
//	protected:
//		virtual void OnRTPPacket( RTPPacket *pack, const RTPTime &receivetime, const RTPAddress *senderaddress );
//		virtual void ProcessRTPPacket(const RTPSourceData &srcdat,const RTPPacket &rtppack);
//		virtual void OnPollThreadStep();
//		virtual void OnBYEPacket( RTPSourceData *srcdat );
//		void OnRemoveSource( RTPSourceData *dat );
//
//	private:
//		//void send_ping();
//		//void start_timer();
//		//bool has_timed_out();
//		//void set_timed_out( bool time_out );
//		void delete_old_packets( RTPPacket* pack );
//		void init_service();
//		int report_error( int errCode );
//
//		std::string get_local_ip_address();
//
//		BinaryBuffer* reconstruct_buffer( unsigned char* data, size_t data_len, size_t total_size, uint32_t timestamp, uint8_t fragment_number );
//
//		std::string m_peer_ip_address_;
//
//		static unsigned short s_client_port_;
//		unsigned short m_port_;
//		
//		unsigned short m_client_rtp_port_;
//		unsigned short m_client_rctp_port_;
//
//		boost::thread* m_thread_;
//		boost::mutex m_mutex_;
//
//		//RTPSession* rtp_session_;
//		RTPSessionParams* sess_params_;
//		RTPUDPv4TransmissionParams* trans_params_;
//
//		unsigned int m_last_timestamp_displayed_;
//		std::vector< BinaryBuffer* > m_binary_buffers_;
//
//		//boost::mutex m_packet_mutex_;
//		//std::vector< RTPPacket* > m_packets_;
//		std::vector< RTPPacket* > m_old_packets_;
//		MediaMemoryManager* m_mem_mgr_;
//
//		bool m_send_ping_;
//		boost::mutex m_send_ping_mutex_;
//
//		bool m_time_out_;
//		boost::mutex m_time_out_mutex_;
//
//		std::map<unsigned int, unsigned int> m_receivers_timestamp_;
//		std::map<unsigned int, boost::shared_ptr< ReceiverBuffer > > m_receivers_buffer_;
//
//		std::string m_local_ip_address_;
//
//		std::ofstream port_writer;
//	};
//
//}// end namespace RTPMediaCommunication
//
//#endif 