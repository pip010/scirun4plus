#ifndef SERVER_MESSAGE_HANDLER_H
#define SERVER_MESSAGE_HANDLER_H

#include "MessageHandler.h"

namespace BioMesh3d
{
	class Message;

	class ServerMessageHandler : public MessageHandler
	{
	public:
		ServerMessageHandler();
		~ServerMessageHandler();

	protected:

		virtual void handle_set_hostname( const std::string& hostname ) { std::cout << "handle_set_hostname not connected" << std::endl; }
		virtual void handle_stage_complete( const std::string& full_project_name, 
			unsigned short stage_number, const std::string& user_name ) { std::cout << "handle_stage_complete not connected" << std::endl; }
		virtual void handle_stage_started() { std::cout << "handle_stage_started not connected" << std::endl; }
		virtual void handle_stage_stopped() { std::cout << "handle_stage_stopped not connected" << std::endl; }

		virtual void handle_stage_terminated(const std::string& full_project_name, 
			unsigned short stage_number, const std::string& user_name ) { std::cout << "handle_stage_terminated not connected" << std::endl; }
		virtual void handle_visualization_terminated(const std::string& full_project_name, 
			unsigned short stage_number, const std::string& client_ip_address ) { std::cout << "handle_visualization_terminated not connected" << std::endl; }
		virtual void handle_stage_error(const std::string& full_project_name, 
			unsigned short stage_number, const std::string& error_message,
			const std::string& user_name ) { std::cout << "handle_stage_error not connected" << std::endl; }
		virtual void handle_visualization_started( const std::string& full_project_name, 
					const std::string& client_ip_address, 
					unsigned short stage_number, 
			unsigned short port, const std::string& user_name, 
			bool is_public_vis_session ){ std::cout << "handle_visualization_started not connected" << std::endl; }
		
	private:
		void register_handlers();
		void unregister_handlers();

		void connect_set_hostname_handler( MessageHandle msg );
		void connect_stage_complete_handler( MessageHandle msg );
		void connect_stage_started_handler( MessageHandle msg );
		void connect_stage_stopped_handler( MessageHandle msg );
		void connect_stage_error_handler( MessageHandle msg );
		void connect_stage_terminated_handler( MessageHandle msg );
		void connect_visualization_terminated_handler( MessageHandle msg );
		void connect_visualization_started_handler( MessageHandle msg );

	};

}

#endif