#include "ServerMessageComposer.h"
#include "Message.h"

namespace BioMesh3d
{
	ServerMessageComposer::ServerMessageComposer()
	{
	}

	ServerMessageComposer::~ServerMessageComposer()
	{
	}

	MessageHandle ServerMessageComposer::compose_set_hostname_message( const std::string& hostname )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVER_SET_HOSTNAME_E );
		msg->add_message_parameter( ParameterTypes::HOSTNAME_E, hostname );
		return msg;
	}

	/*MessageHandle SCIClientMessageComposer::compose_stage_completed_message( 
		const std::string full_project_name, const unsigned short stage_number )*/
	MessageHandle ServerMessageComposer::compose_stage_completed_message( 
		const std::string& user_name,
		const std::string& full_project_name, 
		const unsigned short stage_number )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVER_STAGE_COMPLETE_E );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, user_name );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		return msg;
	}

	MessageHandle ServerMessageComposer::compose_stage_started_message()
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVER_STAGE_STARTED_E );
		return msg;
	}

	MessageHandle ServerMessageComposer::compose_stage_stopped_message()
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVER_STAGE_STOPPED_E );
		return msg;	
	}


	MessageHandle ServerMessageComposer::compose_stage_terminated_message( 
		const std::string& full_project_name, 
		const unsigned short stage_number,
		const std::string& user_name)
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVER_STAGE_TERMINATED_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, user_name );
		return msg;
	}

	MessageHandle ServerMessageComposer::compose_visualization_terminated_message( const std::string& full_project_name, 
		const unsigned short stage_number, const std::string& client_ip_address )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVER_VISUALIZATION_TERMINATED_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::CLIENT_IP_ADDRESS_E, client_ip_address );
		return msg;
	}

	MessageHandle ServerMessageComposer::compose_stage_error_message( 
		const std::string& full_project_name, const unsigned short stage_number,
		const std::string& error_message,
		const std::string& user_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVER_STAGE_ERROR_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::ERROR_MESSAGE_E, error_message );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, user_name );
		return msg;
	}

	MessageHandle ServerMessageComposer::compose_visualization_started_message( 
		const std::string& full_project_name, 
		const std::string& client_ip_address,
		const unsigned short stage_number,
		const unsigned short port,
		const std::string& user_name,
		bool is_public_vis_session )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVER_VISUALIZATION_STARTED_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::CLIENT_IP_ADDRESS_E, client_ip_address );
		msg->add_message_parameter( ParameterTypes::PORT_E, port );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, user_name );
		msg->add_message_parameter( 
			ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session );
		return msg;
	}

}// end namespace Message