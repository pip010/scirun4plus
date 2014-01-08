#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "Message.h"
#include "BrokerMessageHandler.h"

namespace BioMesh3d
{

	BrokerMessageHandler::BrokerMessageHandler()
	{
		register_handlers();
	}

	BrokerMessageHandler::~BrokerMessageHandler()
	{
		unregister_handlers();
	}

	void BrokerMessageHandler::register_handlers()
	{
		register_callback( MessageTypes::SERVICES_MANAGER_IS_CLIENT_AUTHENTICATED_E, 
			boost::bind( &BrokerMessageHandler::connect_is_client_authenticated_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_INITIATE_FILE_UPLOAD_E, 
			boost::bind( &BrokerMessageHandler::connect_initiate_file_upload_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_MESH_STATUS_E, 
			boost::bind( &BrokerMessageHandler::connect_mesh_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_SERVER_LIST_E, 
			boost::bind( &BrokerMessageHandler::connect_server_list_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_START_FILE_DOWNLOAD_E, 
			boost::bind( &BrokerMessageHandler::connect_start_file_download_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_REQUEST_RUN_STAGE_E, 
			boost::bind( &BrokerMessageHandler::connect_request_run_stage_message_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_REQUEST_STOP_STAGE_E, 
			boost::bind( &BrokerMessageHandler::connect_request_stop_stage_message_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_REQUEST_TERMINATE_STAGE_E, 
			boost::bind( &BrokerMessageHandler::connect_request_terminate_stage_message_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_REQUEST_VISUALIZE_STAGE_E, 
			boost::bind( &BrokerMessageHandler::connect_request_visualize_stage_message_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_MODEL_CONFIG_INFO_E, 
			boost::bind( &BrokerMessageHandler::connect_model_config_info_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_FILE_TRANSFER_COMPLETE_E, 
			boost::bind( &BrokerMessageHandler::connect_file_transfer_complete_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_TRANSFER_FILE_E, 
			boost::bind( &BrokerMessageHandler::connect_transfer_file_fragment_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_VISUALIZATION_STARTED_E, 
			boost::bind( &BrokerMessageHandler::connect_visualization_started_message_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_ANSWER_QUERY_PRIMARY_PROJECT_E, 
			boost::bind( &BrokerMessageHandler::connect_answer_primary_project_message_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_ANSWER_QUERY_GROUPS_E, 
			boost::bind( &BrokerMessageHandler::connect_answer_user_groups_message_handler, this, _1 ) );
		register_callback( MessageTypes::BROKER_DELETE_LABEL_MAP_E, 
			boost::bind( &BrokerMessageHandler::connect_delete_label_map_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_DELETE_MESH_E, 
			boost::bind( &BrokerMessageHandler::connect_delete_mesh_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_ADD_MESH_E, 
			boost::bind( &BrokerMessageHandler::connect_add_mesh_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_ANSWER_CREATE_PRIMARY_PROJECT_E, 
			boost::bind( &BrokerMessageHandler::connect_label_map_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_PROJECT_PATH_E, 
			boost::bind( &BrokerMessageHandler::connect_set_project_path_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_CLIENT_LOGIN_FAIL_E, 
			boost::bind( &BrokerMessageHandler::connect_login_fail_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_VISUALIZATION_SESSION_E, 
			boost::bind( &BrokerMessageHandler::connect_vis_notify_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_LEAVE_VIS_SESSION_E, 
			boost::bind( &BrokerMessageHandler::connect_notify_leave_vis_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_JOIN_VIS_SESSION_E, 
			boost::bind( &BrokerMessageHandler::connect_notify_join_vis_handler, this, _1 ) );
		register_callback( MessageTypes::PING_E, 
			boost::bind( &BrokerMessageHandler::connect_ping_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_VISUALIZATION_RESUMED_E, 
			boost::bind( &BrokerMessageHandler::connect_visualization_resumed_message_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_ALERT_MESSAGE_E,
			boost::bind( &BrokerMessageHandler::connect_alert_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_FILE_TRANSFER_CANCELED_E,
			boost::bind( &BrokerMessageHandler::connect_confirm_file_transfer_cancel_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_RSA_PEM_PUBLIC_KEY_E,
			boost::bind( &BrokerMessageHandler::connect_rsa_pem_public_key_handler, this, _1 ) );
		register_callback( MessageTypes::SERVICES_MANAGER_USERS_INFO_E,
			boost::bind( &BrokerMessageHandler::connect_users_info_handler, this, _1 ) );
	}

	void BrokerMessageHandler::unregister_handlers()
	{
		unregister_callback( MessageTypes::SERVICES_MANAGER_IS_CLIENT_AUTHENTICATED_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_INITIATE_FILE_UPLOAD_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_MESH_STATUS_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_SERVER_LIST_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_START_FILE_DOWNLOAD_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_REQUEST_RUN_STAGE_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_REQUEST_TERMINATE_STAGE_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_REQUEST_VISUALIZE_STAGE_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_ANSWER_QUERY_PRIMARY_PROJECT_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_ANSWER_QUERY_GROUPS_E );
		unregister_callback( MessageTypes::BROKER_DELETE_LABEL_MAP_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_DELETE_MESH_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_ADD_MESH_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_ANSWER_CREATE_PRIMARY_PROJECT_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_PROJECT_PATH_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_CLIENT_LOGIN_FAIL_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_VISUALIZATION_SESSION_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_LEAVE_VIS_SESSION_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_JOIN_VIS_SESSION_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_VISUALIZATION_RESUMED_E );
		unregister_callback( MessageTypes::PING_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_ALERT_MESSAGE_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_FILE_TRANSFER_CANCELED_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_RSA_PEM_PUBLIC_KEY_E );
		unregister_callback( MessageTypes::SERVICES_MANAGER_USERS_INFO_E );

	}


	void BrokerMessageHandler::connect_is_client_authenticated_handler( MessageHandle msg )
	{
		bool authenticated = true;
		std::string validation_message = "";
		if ( msg->retrieve_message_parameter( ParameterTypes::AUTHENTICATED_E, authenticated ) && 
			 msg->retrieve_message_parameter( ParameterTypes::VALIDATION_MESSAGE_E, validation_message ))
		{
			handle_is_client_authenticated( authenticated, validation_message );
		}
	}


	void BrokerMessageHandler::connect_initiate_file_upload_handler( MessageHandle msg )
	{
		std::string full_project_name = "";
		bool is_valid = msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		if ( is_valid )
		{
			handle_initiate_file_upload( full_project_name );
		}
	}

	void BrokerMessageHandler::connect_mesh_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name, error_message;
		unsigned short stage_number;
		int stage_state_int;
		StageStates::StageState stage_state;
		unsigned long date_processed;
		std::string mesh_description;

		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) && 
			 msg->retrieve_message_parameter( ParameterTypes::DATE_PROCESSED_E, date_processed )&& 
			 msg->retrieve_message_parameter( ParameterTypes::ERROR_MESSAGE_E, error_message ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_STATE_E, stage_state_int ) &&
			 msg->retrieve_message_parameter( ParameterTypes::DESCRIPTION_E, mesh_description ) )
		{
			stage_state = ( StageStates::StageState )stage_state_int;
			handle_mesh_status( label_map_name, mesh_name, stage_state, 
				stage_number, date_processed, error_message, mesh_description );
		}
	}

	void BrokerMessageHandler::connect_server_list_handler( MessageHandle msg )
	{
		std::string server_ip_address_list_string = "";
		std::string server_hostname_list_string = "";
		std::vector< std::string > server_ip_address_list;
		std::vector< std::string > server_hostname_list;
		if ( msg->retrieve_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_LIST_E, server_ip_address_list_string ) &&
			 msg->retrieve_message_parameter( ParameterTypes::SERVER_HOSTNAME_LIST_E, server_hostname_list_string ) )
		{
			boost::algorithm::split( server_ip_address_list, server_ip_address_list_string, boost::is_any_of(",") );
			boost::algorithm::split( server_hostname_list, server_hostname_list_string, boost::is_any_of(",") );
			handle_server_list( server_ip_address_list, server_hostname_list );
		}
	}

	void BrokerMessageHandler::connect_start_file_download_handler( MessageHandle msg )
	{
		std::string full_project_name, file_name = "";
		unsigned int file_size;
		if ( msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FILE_NAME_E, file_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FILE_SIZE_E, file_size ) )
		{
			handle_start_file_download( full_project_name, file_name, file_size );
		}
	}

	void BrokerMessageHandler::connect_file_transfer_complete_handler( MessageHandle msg )
	{
		std::string full_project_name = "";
		size_t file_size = 0;
		unsigned char* file_data = msg->get_binary_buf();
		if ( msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FILE_SIZE_E, file_size ) &&
			 file_data != NULL )
		{
			handle_file_transfer_complete( full_project_name, file_data, file_size );
		}
	}

	void BrokerMessageHandler::connect_transfer_file_fragment_handler( MessageHandle msg )
	{
		std::string file_path, label_map_name, mesh_name;
		size_t file_size = 0;
		unsigned short fragment_size = 0;
		unsigned char* file_data = msg->get_binary_buf();

		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FILE_NAME_E, file_path ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FILE_SIZE_E, file_size ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FRAGMENT_SIZE_E, fragment_size ) &&
			 file_data != NULL )
		{
			handle_transfer_file_fragment( label_map_name, mesh_name, file_path, file_size, file_data, fragment_size );
		}
	}

	void BrokerMessageHandler::connect_model_config_info_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name, description, mats, mat_names, 
			mat_radii, refinement_levels, constant_sizing_value, max_sizing_field, 
			num_particle_iters, tetgen_joined_vol_flags = "";
		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) &&
			msg->retrieve_message_parameter( ParameterTypes::DESCRIPTION_E, description ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MATS_E, mats ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MAT_NAMES_E, mat_names ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MAT_RADII_E, mat_radii ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_REFINEMENT_LEVELS_E, refinement_levels ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_CONSTANT_SIZING_VALUE_E, constant_sizing_value ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MAX_SIZING_FIELD_E, max_sizing_field ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_NUM_PARTICLE_ITERS_E, num_particle_iters ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_TETGEN_JOINED_VOL_FLAGS_E, tetgen_joined_vol_flags ) )
		{
			handle_model_config_info( label_map_name, mesh_name, description, mats, mat_names, mat_radii, 
				refinement_levels, constant_sizing_value, max_sizing_field, num_particle_iters, tetgen_joined_vol_flags );
		}
	}

	inline bool stage_message_helper( MessageHandle msg, std::string& full_project_name, 
		unsigned short& stage_number )
	{
		return msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name ) &&
			msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
	}
	
	void BrokerMessageHandler::connect_request_run_stage_message_handler( MessageHandle msg )
	{
		std::string label_map, mesh;
		unsigned short start_stage_number, end_stage_number;
		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh ) &&
			 msg->retrieve_message_parameter( ParameterTypes::START_STAGE_NUMBER_E, start_stage_number ),
			 msg->retrieve_message_parameter( ParameterTypes::END_STAGE_NUMBER_E, end_stage_number ) )
		{
			handle_request_run_stage( label_map, mesh, start_stage_number, end_stage_number );
		}
	}

	void BrokerMessageHandler::connect_request_stop_stage_message_handler( MessageHandle msg )
	{
		std::string label_map, mesh;
		unsigned short start_stage_number, end_stage_number;
		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh ) )
		{
			handle_request_stop_stage( label_map, mesh );
		}
	}

	void BrokerMessageHandler::connect_request_terminate_stage_message_handler( MessageHandle msg )
	{
		std::string label_map, mesh;
		unsigned short stage_number;
		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number ) )
		{
			handle_request_terminate_stage( label_map, mesh, stage_number );
		}
	}

	void BrokerMessageHandler::connect_request_visualize_stage_message_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name, client_ip_address = "";
		std::string user_name;
		unsigned short stage_number;
		bool is_public_vis_session;


		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::CLIENT_IP_ADDRESS_E, client_ip_address ) &&
			 msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, user_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session ) )
		{
			handle_request_visualize_stage_message( label_map_name,
				mesh_name, stage_number, client_ip_address,
				user_name, is_public_vis_session );
		}
	}

	void BrokerMessageHandler::connect_visualization_started_message_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name, server_ip_address = "";
		unsigned short stage_number, port = 0;
		std::string user_name;
		bool is_public_vis_session;
		unsigned short group_id;

		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::PORT_E, port ) &&
			 msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, user_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session ) &&
			 msg->retrieve_message_parameter( ParameterTypes::GROUP_NAME_E, group_id ) )
		{
			handle_visualization_started_message( 
				label_map_name, mesh_name, group_id, server_ip_address,
				stage_number, port, user_name, is_public_vis_session );
		}
	}

	void BrokerMessageHandler::connect_answer_primary_project_message_handler( MessageHandle msg )
	{
		bool is_find;
		if ( msg->retrieve_message_parameter( ParameterTypes::PRIMARY_PROJECT_FIND_E, is_find ) )
		{
			handle_answer_primary_project_message( is_find );
		}
	}

	void BrokerMessageHandler::connect_answer_user_groups_message_handler( MessageHandle msg )
	{
		std::string groups_string = "";
		bool is_valid = msg->retrieve_message_parameter( ParameterTypes::USER_GROUPS_E, groups_string );
		if ( is_valid )
		{
			std::vector< std::string > groups;
			boost::algorithm::split( groups, groups_string, boost::is_any_of( "," ) );
			handle_answer_user_groups_message( groups );
		}
	}

	void BrokerMessageHandler::connect_delete_mesh_handler ( MessageHandle msg )
	{
		std::string label_map_name, mesh_name;

		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ))
		{
			handle_delete_mesh_message( label_map_name, mesh_name );
		}
	}

	void BrokerMessageHandler::connect_delete_label_map_handler( 
		MessageHandle msg )
	{
		std::string label_map_name;

		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) )
		{
			handle_delete_label_map_message( label_map_name );
		}

	}

	void BrokerMessageHandler::connect_add_mesh_handler( 
								MessageHandle msg )
	{
		std::string label_map_name, mesh_name, error_message, mesh_description;
		unsigned int stage_number, stage_state;
		unsigned long date_processed;

		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_STATE_E, stage_state ) &&
			 msg->retrieve_message_parameter( ParameterTypes::ERROR_MESSAGE_E, error_message ) &&
			 msg->retrieve_message_parameter( ParameterTypes::DESCRIPTION_E, mesh_description ) &&
			 msg->retrieve_message_parameter( ParameterTypes::DATE_PROCESSED_E, date_processed ) )
		{
			StageStates::StageState stage_state_tmp = ( StageStates::StageState )stage_state;
			handle_add_mesh_message( label_map_name, mesh_name, stage_number, 
				stage_state_tmp, date_processed, error_message, mesh_description );
		}

	}


	void BrokerMessageHandler::connect_label_map_handler( 
		MessageHandle msg )
	{
		std::string label_map_name;
		unsigned long date_processed;
		std::string label_map_description;

		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::DATE_PROCESSED_E, date_processed ) &&
			 msg->retrieve_message_parameter( ParameterTypes::DESCRIPTION_E, label_map_description ) )
		{
			handle_label_map_message( label_map_name, date_processed, label_map_description  );
		}

	}

	void BrokerMessageHandler::connect_set_project_path_handler( MessageHandle msg )
	{
		std::string project_path;
		if ( msg->retrieve_message_parameter( ParameterTypes::PROJECT_PATH_E, project_path ) )
		{
			handle_set_project_path_message( project_path );
		}
	}

	void BrokerMessageHandler::connect_login_fail_handler( MessageHandle msg )
	{
		std::string error_description;
		if ( msg->retrieve_message_parameter( 
			ParameterTypes::LOGIN_FAIL_DESCRIPTION_E, error_description ) )
		{
			handle_client_login_fail_message( error_description );
		}
	}

	void BrokerMessageHandler::connect_vis_notify_handler( MessageHandle msg )
	{
		std::string full_project_name, server_ip_address = "";
		unsigned short stage_number, port = 0;
		std::string owner_name;
		bool is_public_vis_session;

		unsigned short group_id;
		std::string members;
		unsigned int sub_project_id;

		if ( msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::PORT_E, port ) &&
			 msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, owner_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session ) &
			 msg->retrieve_message_parameter( ParameterTypes::GROUP_NAME_E, group_id ) &&
			 msg->retrieve_message_parameter( ParameterTypes::VIS_MEMBERS_E, members ) &&
			 msg->retrieve_message_parameter( ParameterTypes::SUB_PROJECT_ID_E, sub_project_id ) )
		{
			handle_vis_notify_message( 
				full_project_name, server_ip_address,
				stage_number, port, owner_name, is_public_vis_session,
				group_id, members, sub_project_id );
		}
	}

	void BrokerMessageHandler::connect_notify_leave_vis_handler( MessageHandle msg )
	{
		std::string full_project_name, server_ip_address = "";
		unsigned short stage_number, port = 0;
		std::string owner_name;
		std::string member_name;
		bool is_public_vis_session;

		unsigned short group_id;
		unsigned int sub_project_id;

		if ( msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::PORT_E, port ) &&
			 msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, owner_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MEMBERNAME_E, member_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session ) &&
			 msg->retrieve_message_parameter( ParameterTypes::GROUP_NAME_E, group_id ) &&
			 msg->retrieve_message_parameter( ParameterTypes::SUB_PROJECT_ID_E, sub_project_id ) )
		{
			handle_notify_leave_vis_message( 
				full_project_name, server_ip_address,
				stage_number, port, owner_name, member_name, is_public_vis_session,
				group_id, sub_project_id );
		}
	}

	void BrokerMessageHandler::connect_notify_join_vis_handler( MessageHandle msg )
	{
		std::string full_project_name, server_ip_address = "";
		unsigned short stage_number, port = 0;
		std::string member_name;
		bool is_public_vis_session;

		unsigned short group_id;
		unsigned int sub_project_id;

		if ( msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::PORT_E, port ) &&
			 msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, member_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session ) &&
			 msg->retrieve_message_parameter( ParameterTypes::GROUP_NAME_E, group_id ) &&
			 msg->retrieve_message_parameter( ParameterTypes::SUB_PROJECT_ID_E, sub_project_id ) )
		{
			handle_notify_join_vis_message( 
				full_project_name, server_ip_address,
				stage_number, port, member_name, is_public_vis_session,
				group_id, sub_project_id );
		}
	}

	void BrokerMessageHandler::connect_visualization_resumed_message_handler( MessageHandle msg )
	{
		std::string full_project_name, server_ip_address = "";
		unsigned short stage_number, port = 0;
		std::string user_name;
		bool is_public_vis_session;
		unsigned short group_id;

		if ( msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::PORT_E, port ) &&
			 msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, user_name ) &&
			 msg->retrieve_message_parameter(ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session ) &&
			 msg->retrieve_message_parameter(ParameterTypes::GROUP_NAME_E, group_id ) )
		{
			handle_visualization_resumed_message( 
				full_project_name, group_id, server_ip_address,
				stage_number, port, user_name, is_public_vis_session );
		}
	}

	void BrokerMessageHandler::connect_ping_handler( MessageHandle msg )
	{
		handle_ping_message();
	}
	
	void BrokerMessageHandler::connect_alert_handler( MessageHandle msg )
	{
		std::string alert_string;
		if ( msg->retrieve_message_parameter( ParameterTypes::ALERT_E, alert_string ) )
		{
			handle_alert_message( alert_string );
		}
	}

	void BrokerMessageHandler::connect_confirm_file_transfer_cancel_handler( MessageHandle msg )
	{
		handle_confirm_file_transfer_cancel_message();
	}

	void BrokerMessageHandler::connect_rsa_pem_public_key_handler( MessageHandle msg )
	{
		std::string rsa_pem_public_key_string = "";
		if ( msg->retrieve_message_parameter( ParameterTypes::RSA_PEM_PUBLIC_KEY_E, rsa_pem_public_key_string ) )
		{
			handle_rsa_pem_public_key_message( rsa_pem_public_key_string );
		}
	}

	void BrokerMessageHandler::connect_users_info_handler( MessageHandle msg )
	{
		std::vector< std::string > usernames, first_names, last_names;
		std::vector< std::vector< std::string > > groups;
		std::string parameters = "";
		if ( msg->retrieve_message_parameter( ParameterTypes::USERS_INFO_E, parameters ) )
		{
			std::vector< std::string > users_strings;
			boost::algorithm::split( users_strings, parameters, boost::is_any_of( "\n" ) );
			groups.resize( users_strings.size() - 1 ); // we subtract the last one because it's always an empty string
			for ( size_t i = 0; i < users_strings.size() - 1; ++i )
			{
				std::vector< std::string > user_info;
				boost::algorithm::split( user_info, users_strings[i], boost::is_any_of( "\r" ) );
				usernames.push_back( user_info[0] );
				first_names.push_back( user_info[1] );
				last_names.push_back( user_info[2] );
				std::vector< std::string > user_groups;
				boost::algorithm::split( user_groups, user_info[3], boost::is_any_of( "," ) );
				groups[i] = user_groups;
			}
			handle_users_info( usernames, first_names, last_names, groups );
		}
	}

}// end namespace Message