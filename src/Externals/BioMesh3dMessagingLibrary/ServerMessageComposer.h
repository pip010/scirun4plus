#ifndef SERVER_MESSAGE_COMPOSER_H
#define SERVER_MESSAGE_COMPOSER_H

#include "Message.h"

namespace BioMesh3d
{

	class ServerMessageComposer
	{
	public:
		ServerMessageComposer();
		~ServerMessageComposer();

	protected:

		MessageHandle compose_set_hostname_message( const std::string& hostname );

		MessageHandle compose_stage_completed_message( 
			const std::string& user_name,
			const std::string& full_project_name,
			const unsigned short stage_number );

		MessageHandle compose_stage_started_message();

		MessageHandle compose_stage_stopped_message();

		MessageHandle compose_stage_terminated_message( 
			const std::string& full_project_name, 
			const unsigned short stage_number,
			const std::string& user_name );

		MessageHandle compose_visualization_terminated_message( 
			const std::string& full_project_name, const unsigned short stage_number, 
			const std::string& client_ip_address );

		MessageHandle compose_stage_error_message( 
			const std::string& full_project_name, 
			const unsigned short stage_number,
			const std::string& error_message,
			const std::string& user_name );

		MessageHandle compose_visualization_started_message( 
			const std::string& full_project_name,
			const std::string& client_ip_address,
			const unsigned short stage_number,
			const unsigned short port, const std::string& user_name,
			bool is_public_vis_session );
	private:
	};

}// end namespace Message
#endif
