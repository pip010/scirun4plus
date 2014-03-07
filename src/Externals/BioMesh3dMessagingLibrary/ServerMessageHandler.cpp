#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>

#include "ServerMessageHandler.h"
#include "Message.h"

namespace BioMesh3d
{

	ServerMessageHandler::ServerMessageHandler()
	{
		register_handlers();
	}

	ServerMessageHandler::~ServerMessageHandler()
	{
		unregister_handlers();
	}

	void ServerMessageHandler::register_handlers()
	{
		register_callback( MessageTypes::SERVER_SET_HOSTNAME_E, 
			boost::bind( &ServerMessageHandler::connect_set_hostname_handler, this, _1 ) );
		register_callback( MessageTypes::SERVER_STAGE_COMPLETE_E, 
			boost::bind( &ServerMessageHandler::connect_stage_complete_handler, this, _1 ) );
		register_callback( MessageTypes::SERVER_STAGE_STARTED_E, 
			boost::bind( &ServerMessageHandler::connect_stage_started_handler, this, _1 ) );
		register_callback( MessageTypes::SERVER_STAGE_STOPPED_E, 
			boost::bind( &ServerMessageHandler::connect_stage_stopped_handler, this, _1 ) );
		register_callback( MessageTypes::SERVER_STAGE_ERROR_E, 
			boost::bind( &ServerMessageHandler::connect_stage_error_handler, this, _1 ) );
		register_callback( MessageTypes::SERVER_STAGE_TERMINATED_E, 
			boost::bind( &ServerMessageHandler::connect_stage_terminated_handler, this, _1 ) );
		register_callback( MessageTypes::SERVER_VISUALIZATION_TERMINATED_E, 
			boost::bind( &ServerMessageHandler::connect_visualization_terminated_handler, this, _1 ) );
		register_callback( MessageTypes::SERVER_VISUALIZATION_STARTED_E, 
			boost::bind( &ServerMessageHandler::connect_visualization_started_handler, this, _1 ) );
	}

	void ServerMessageHandler::unregister_handlers()
	{
		unregister_callback( MessageTypes::SERVER_SET_HOSTNAME_E );
		unregister_callback( MessageTypes::SERVER_STAGE_COMPLETE_E );
		unregister_callback( MessageTypes::SERVER_STAGE_STARTED_E );
		unregister_callback( MessageTypes::SERVER_STAGE_STOPPED_E );
		unregister_callback( MessageTypes::SERVER_STAGE_ERROR_E );
		unregister_callback( MessageTypes::SERVER_STAGE_TERMINATED_E );
		unregister_callback( MessageTypes::SERVER_VISUALIZATION_TERMINATED_E );
	}

	void ServerMessageHandler::connect_set_hostname_handler( MessageHandle msg )
	{
		std::string hostname = "";
		if ( msg->retrieve_message_parameter( ParameterTypes::HOSTNAME_E, hostname ) )
		{
			handle_set_hostname( hostname );
		}
	}

	inline bool stage_message_helper( MessageHandle msg, std::string& full_project_name, 
		unsigned short& stage_number )
	{
		return msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name ) &&
			msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
	}

	void ServerMessageHandler::connect_stage_complete_handler( MessageHandle msg )
	{
		std::string full_project_name = "";
		unsigned short stage_number = 0;
		std::string user_name;
		if ( stage_message_helper( msg, full_project_name, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, user_name ) )
		{
			handle_stage_complete( full_project_name, stage_number, user_name );
		}
	}

	void ServerMessageHandler::connect_stage_started_handler( MessageHandle msg )
	{
		handle_stage_started();
	}

	void ServerMessageHandler::connect_stage_stopped_handler( MessageHandle msg )
	{
		handle_stage_stopped();
	}

	
	void ServerMessageHandler::connect_stage_terminated_handler( MessageHandle msg )
	{
		std::string full_project_name = "";
		unsigned short stage_number = 0;
		std::string user_name;
		if ( stage_message_helper( msg, full_project_name, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, user_name ) )
		{
			handle_stage_terminated( full_project_name, stage_number, user_name );
		}
	}
	
	void ServerMessageHandler::connect_stage_error_handler( MessageHandle msg )
	{
		std::string full_project_name, error_message = "";
		unsigned short stage_number = 0;
		std::string user_name;

		if ( stage_message_helper( msg, full_project_name, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::ERROR_MESSAGE_E, error_message ) &&
			 msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, user_name ) )
		{
			handle_stage_error( full_project_name, stage_number, error_message, user_name );
		}
	}

	void ServerMessageHandler::connect_visualization_terminated_handler( MessageHandle msg )
	{
		std::string full_project_name, client_ip_address = "";
		unsigned short stage_number = 0;
		if ( stage_message_helper( msg, full_project_name, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::CLIENT_IP_ADDRESS_E, client_ip_address ) )
		{
			handle_visualization_terminated( full_project_name, stage_number, client_ip_address );
		}
	}

	void ServerMessageHandler::connect_visualization_started_handler( MessageHandle msg )
	{
		std::string full_project_name, client_ip_address = "";
		unsigned short stage_number, port = 0;
		std::string user_name;
		bool is_public_vis_session;


		if ( stage_message_helper( msg, full_project_name, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::CLIENT_IP_ADDRESS_E, client_ip_address ) &&
			 msg->retrieve_message_parameter( ParameterTypes::PORT_E, port ) &&
			 msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, user_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session ) )
		{
			handle_visualization_started( full_project_name, 
				client_ip_address, stage_number, port,
				user_name, is_public_vis_session );
		}
	}

}// end namespace Message