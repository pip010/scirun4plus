#include "Message.h"
#include "BrokerMessageComposer.h"

namespace BioMesh3d
{
	BrokerMessageComposer::BrokerMessageComposer()
	{
	}

	BrokerMessageComposer::~BrokerMessageComposer()
	{
	}


	MessageHandle BrokerMessageComposer::compose_is_client_authenticated_message( 
		bool authenticated, const std::string& validation_message )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_IS_CLIENT_AUTHENTICATED_E );
		msg->add_message_parameter( ParameterTypes::AUTHENTICATED_E, authenticated );
		msg->add_message_parameter( ParameterTypes::VALIDATION_MESSAGE_E, validation_message );
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_client_login_fail_message( const std::string& error_description )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_CLIENT_LOGIN_FAIL_E );
		msg->add_message_parameter( 
			ParameterTypes::LOGIN_FAIL_DESCRIPTION_E, error_description );
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_initiate_file_upload_message( const std::string& full_project_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_INITIATE_FILE_UPLOAD_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		return msg;
	}

	BioMesh3d::MessageHandle BrokerMessageComposer::compose_mesh_status( 
		const std::string& label_map_name, const std::string& mesh_name, 
		const StageStates::StageState stage_state, const unsigned short stage_number, 
		const unsigned long date_processed, const std::string& error_message, 
		const std::string& mesh_description )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_MESH_STATUS_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::STAGE_STATE_E, ( int )stage_state );
		msg->add_message_parameter( ParameterTypes::DATE_PROCESSED_E, date_processed );
		msg->add_message_parameter( ParameterTypes::ERROR_MESSAGE_E, error_message );
		msg->add_message_parameter( ParameterTypes::DESCRIPTION_E, mesh_description );
		
		return msg;
	}

	BioMesh3d::MessageHandle BrokerMessageComposer::compose_model_config_info_message( 
		const std::string& label_map_name, const std::string& mesh_name, const 
		std::string& description, const std::string& mats, const std::string& mat_names, 
		const std::string& mat_radii, const std::string& refinement_levels, 
		const std::string& constant_sizing_value, const std::string& max_sizing_field, 
		const std::string& num_particle_iters, const std::string& tetgen_joined_vol_flags )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_MODEL_CONFIG_INFO_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::DESCRIPTION_E, description );
		msg->add_message_parameter( ParameterTypes::MODEL_MATS_E, mats );
		msg->add_message_parameter( ParameterTypes::MODEL_MAT_NAMES_E, mat_names );
		msg->add_message_parameter( ParameterTypes::MODEL_MAT_RADII_E, mat_radii );
		msg->add_message_parameter( ParameterTypes::MODEL_REFINEMENT_LEVELS_E, refinement_levels );
		msg->add_message_parameter( ParameterTypes::MODEL_CONSTANT_SIZING_VALUE_E, constant_sizing_value );		
		msg->add_message_parameter( ParameterTypes::MODEL_MAX_SIZING_FIELD_E, max_sizing_field );
		msg->add_message_parameter( ParameterTypes::MODEL_NUM_PARTICLE_ITERS_E, num_particle_iters );
		msg->add_message_parameter( ParameterTypes::MODEL_TETGEN_JOINED_VOL_FLAGS_E, tetgen_joined_vol_flags );
		return msg;
	}


	MessageHandle BrokerMessageComposer::compose_server_list_message( 
		const std::vector< std::string >& server_ip_address_list, const std::vector< std::string >& server_hostname_list )
	{
		std::string server_ip_list_string = "";
		std::string server_hostname_string = "";
		MessageHandle msg = MessageHandle( new Message() );
		assert ( server_ip_address_list.size() == server_hostname_list.size() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_SERVER_LIST_E );
		for ( size_t i = 0; i < server_ip_address_list.size(); ++i )
		{
			std::string comma = ( i < server_ip_address_list.size() - 1 ) ? "," : "";
			server_ip_list_string += server_ip_address_list[i] + comma;
			server_hostname_string += server_hostname_list[i] + comma;
		}
		msg->add_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_LIST_E, server_ip_list_string );
		msg->add_message_parameter( ParameterTypes::SERVER_HOSTNAME_LIST_E, server_hostname_string );
		return msg;
	}


	MessageHandle BrokerMessageComposer::compose_start_file_download_message( const std::string& full_project_name, 
		const std::string& file_name, size_t file_size )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_START_FILE_DOWNLOAD_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::FILE_NAME_E, file_name );
		msg->add_message_parameter( ParameterTypes::FILE_SIZE_E, file_size );
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_nrrd_file_names_message( std::vector< std::string > nrrd_file_names )
	{
		MessageHandle msg = MessageHandle( new Message() );
		std::string file_names_string = "";
		for ( unsigned int i = 0; i < nrrd_file_names.size(); ++i )
		{
			std::string comma = ( i < nrrd_file_names.size() - 1 ) ? "," : "";
			file_names_string += nrrd_file_names[i] + comma;
		}
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_file_transfer_complete_message( const std::string& full_project_name, unsigned char* file_data, size_t file_size )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_FILE_TRANSFER_COMPLETE_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::FILE_SIZE_E, file_size );
		msg->copy_binary_data( file_data, file_size );
		return msg;
	}

	BioMesh3d::MessageHandle BrokerMessageComposer::compose_transfer_file_fragment_message( 
		const std::string& label_map_name, const std::string& mesh_name, const std::string& file_name, 
		size_t file_size, unsigned char* file_data, unsigned short fragment_size )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_TRANSFER_FILE_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::FILE_NAME_E, file_name );
		msg->add_message_parameter( ParameterTypes::FILE_SIZE_E, file_size );
		msg->add_message_parameter( ParameterTypes::FRAGMENT_SIZE_E, fragment_size );
		msg->copy_binary_data( file_data, fragment_size );
		return msg;
	}

	BioMesh3d::MessageHandle BrokerMessageComposer::compose_request_run_stage_message( 
		const std::string& label_map_name, const std::string& mesh_name, 
		unsigned short start_stage_number, unsigned short end_stage_number )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_REQUEST_RUN_STAGE_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::START_STAGE_NUMBER_E, start_stage_number );
		msg->add_message_parameter( ParameterTypes::END_STAGE_NUMBER_E, end_stage_number );
		return msg;
	}

	BioMesh3d::MessageHandle BrokerMessageComposer::compose_request_stop_stage_message( 
		const std::string& label_map_name, const std::string& mesh_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_REQUEST_STOP_STAGE_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_request_terminate_stage_message( 
		const std::string& label_map_name,
		const std::string& mesh_name, 
		unsigned short stage_number )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_REQUEST_TERMINATE_STAGE_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_request_visualize_stage_message( 
		const std::string& project_path,
		const std::string& full_project_name, 
		unsigned short stage_number, const std::string& client_ip_address,
		const std::string& user_name,
		bool is_public_vis_session )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_REQUEST_VISUALIZE_STAGE_E );
		msg->add_message_parameter( ParameterTypes::PROJECT_PATH_E, project_path );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::CLIENT_IP_ADDRESS_E, client_ip_address );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, user_name );
		msg->add_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session );
		return msg;
	}

	BioMesh3d::MessageHandle BrokerMessageComposer::compose_visualization_started_message( 
		const std::string& label_map_name, const std::string& mesh_name, unsigned short group_id, 
		const std::string& server_ip_address, const unsigned short stage_number,
		const unsigned short port, const std::string& user_name, bool is_public_vis_session )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_VISUALIZATION_STARTED_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::GROUP_NAME_E, group_id );
		msg->add_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::PORT_E, port );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, user_name );
		msg->add_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session );

		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_visualization_resumed_message( 
		const std::string& full_project_name, 
		unsigned short group_id,
		const std::string& server_ip_address, 
		const unsigned short stage_number,
		const unsigned short port,
		const std::string& user_name,
		bool is_public_vis_session )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_VISUALIZATION_RESUMED_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::GROUP_NAME_E, group_id );
		msg->add_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::PORT_E, port );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, user_name );
		msg->add_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session );

		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_query_primary_project_message(
				const bool is_find )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( 
			MessageTypes::SERVICES_MANAGER_ANSWER_QUERY_PRIMARY_PROJECT_E );
		msg->add_message_parameter(
			ParameterTypes::PRIMARY_PROJECT_FIND_E, is_find );
		
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_delete_mesh_message( 
		const std::string& label_map_name, const std::string& mesh_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_DELETE_MESH_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );

		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_delete_label_map_message(
		const std::string& label_map_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::BROKER_DELETE_LABEL_MAP_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );

		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_answer_update_subproject_message(
		const bool is_update,  std::string description )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( 
			MessageTypes::SERVICES_MANAGER_ANSWER_UPDATE_SUBPROJECT_E );
		msg->add_message_parameter( 
			ParameterTypes::OPERATION_ANSWER_E, is_update );
		msg->add_message_parameter( 
			ParameterTypes::OPERATION_DESCRIPTION_E, description );

		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_label_map_message(
						const std::string& label_map_name, 
						const unsigned long date_created,
						const std::string& label_map_description )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_ANSWER_CREATE_PRIMARY_PROJECT_E );
		msg->add_message_parameter(	ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::DATE_PROCESSED_E, date_created );
		msg->add_message_parameter( ParameterTypes::DESCRIPTION_E, label_map_description );

		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_add_mesh_message( 
		const std::string& label_map_name, const std::string& mesh_name, 
		const unsigned short stage_number, const StageStates::StageState stage_state, 
		const unsigned long date_processed, const std::string& error_message, 
		const std::string& mesh_description )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( 
			MessageTypes::SERVICES_MANAGER_ADD_MESH_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::STAGE_STATE_E, ( unsigned int )stage_state );
		msg->add_message_parameter( ParameterTypes::DATE_PROCESSED_E, date_processed );
		msg->add_message_parameter( ParameterTypes::ERROR_MESSAGE_E, error_message );
		msg->add_message_parameter( ParameterTypes::DESCRIPTION_E, mesh_description );

		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_user_groups_message(
		const std::vector< std::string >& groups)
	{
		std::string groups_string = "";
		for ( size_t i = 0; i < groups.size(); ++i )
		{
			std::string comma = ( i == groups.size() - 1 ) ? "" : ",";
			groups_string += groups[i] + comma;
		}
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_ANSWER_QUERY_GROUPS_E );

		msg->add_message_parameter( ParameterTypes::USER_GROUPS_E, groups_string );

		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_set_project_path_message(
		const std::string& project_path )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type(
			MessageTypes::SERVICES_MANAGER_PROJECT_PATH_E );

		msg->add_message_parameter( 
			ParameterTypes::PROJECT_PATH_E, project_path );

		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_notify_vis_session_message( 
		const std::string& full_project_name, 
		const std::string& server_ip_address, 
		const unsigned short stage_number,
		const unsigned short port,
		const std::string& owner_name,
		bool is_public_vis_session,
		unsigned short group_id,
		const std::string& members,
		unsigned int sub_project_id)
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_VISUALIZATION_SESSION_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::PORT_E, port );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, owner_name );
		msg->add_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session );
		msg->add_message_parameter( ParameterTypes::GROUP_NAME_E, group_id );
		msg->add_message_parameter( ParameterTypes::SUB_PROJECT_ID_E, sub_project_id );
		msg->add_message_parameter( ParameterTypes::VIS_MEMBERS_E, members );
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_notify_leave_vis_session_message( 
		const std::string& full_project_name, 
		const std::string& server_ip_address, 
		const unsigned short stage_number,
		const unsigned short port,
		const std::string& owner_name,
		const std::string& member_name,
		bool is_public_vis_session,
		unsigned short group_id,
		unsigned int sub_project_id)
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_LEAVE_VIS_SESSION_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::PORT_E, port );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, owner_name );
		msg->add_message_parameter( ParameterTypes::MEMBERNAME_E, member_name );
		msg->add_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session );
		msg->add_message_parameter( ParameterTypes::GROUP_NAME_E, group_id );
		msg->add_message_parameter( ParameterTypes::SUB_PROJECT_ID_E, sub_project_id );
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_notify_join_vis_session_message( 
		const std::string& full_project_name, 
		const std::string& server_ip_address, 
		const unsigned short stage_number,
		const unsigned short port,
		const std::string& owner_name,
		bool is_public_vis_session,
		unsigned short group_id,
		unsigned int sub_project_id)
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_JOIN_VIS_SESSION_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::PORT_E, port );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, owner_name );
		msg->add_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session );
		msg->add_message_parameter( ParameterTypes::GROUP_NAME_E, group_id );
		msg->add_message_parameter( ParameterTypes::SUB_PROJECT_ID_E, sub_project_id );
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_ping_message()
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::PING_E );
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_alert_message( const std::string& alert_string )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_ALERT_MESSAGE_E );
		msg->add_message_parameter( ParameterTypes::ALERT_E, alert_string );
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_confirm_file_transfer_cancel_message()
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_FILE_TRANSFER_CANCELED_E );
		return msg;
	}

	MessageHandle BrokerMessageComposer::compose_rsa_pem_public_key_message( 
		const std::string& rsa_pem_public_key_string )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_RSA_PEM_PUBLIC_KEY_E );
		msg->add_message_parameter( ParameterTypes::RSA_PEM_PUBLIC_KEY_E, 
			rsa_pem_public_key_string );
		return msg;
	}

	BioMesh3d::MessageHandle BrokerMessageComposer::compose_noop_message()
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::NOOP_E );
		return msg;
	}

	BioMesh3d::MessageHandle BrokerMessageComposer::compose_users_info( const std::vector< std::string >& usernames, 
		const std::vector< std::string >& first_names, const std::vector< std::string >& last_names, 
		const std::vector< std::vector< std::string > >& groups )
	{
		std::string parameters = "";
		for ( size_t i = 0; i < usernames.size(); ++i )
		{
			parameters += usernames[i] + '\r' + first_names[i] + '\r' + last_names[i] + '\r';
			for ( size_t j = 0; j < groups.at( i ).size(); ++j )
			{
				char comma_or_end = ( j < groups.at( i ).size() - 1 ) ? ',' : '\n';
				parameters += groups.at( i ).at( j ) + comma_or_end;
			}
		}
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::SERVICES_MANAGER_USERS_INFO_E );
		msg->add_message_parameter( ParameterTypes::USERS_INFO_E,  parameters );
		return msg;
	}

}