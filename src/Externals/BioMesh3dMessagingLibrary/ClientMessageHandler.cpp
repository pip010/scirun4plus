#include <boost/bind.hpp>

#include "Message.h"
#include "ClientMessageHandler.h"

namespace BioMesh3d
{

	ClientMessageHandler::ClientMessageHandler()
	{
		register_handlers();
	}

	ClientMessageHandler::~ClientMessageHandler()
	{
		unregister_handlers();
	}

	void ClientMessageHandler::register_handlers()
	{

		register_callback( MessageTypes::CLIENT_KEYBOARD_E, 
			boost::bind( &ClientMessageHandler::connect_keyboard_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_MOUSE_E, 
			boost::bind( &ClientMessageHandler::connect_mouse_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_CLIPPING_E, 
			boost::bind( &ClientMessageHandler::connect_clipping_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_RESIZE_E, 
			boost::bind( &ClientMessageHandler::connect_resize_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_QUALITY_ADJUST_E,
			boost::bind( &ClientMessageHandler::connect_quality_adjust_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_EXIT_SCIRUN_E,
			boost::bind( &ClientMessageHandler::connect_exit_scirun_handler, this, _1 ) );
		register_callback( MessageTypes::PING_E,
			boost::bind( &ClientMessageHandler::connect_ping_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_QUALITY_OF_SERVICE_E,
			boost::bind( &ClientMessageHandler::connect_quality_of_service_handler, this, _1 ) );

		register_callback( MessageTypes::CLIENT_AUTHENTICATE_CLIENT_E,
			boost::bind( &ClientMessageHandler::connect_authenticate_client_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_ADD_USER_E,
			boost::bind( &ClientMessageHandler::connect_add_user_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_ADD_LABEL_MAP_E,
			boost::bind( &ClientMessageHandler::connect_add_label_map_message_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_ADD_MESH_E,
			boost::bind( &ClientMessageHandler::connect_add_mesh_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_UPDATE_MODEL_CONFIG_INFO_E,
			boost::bind( &ClientMessageHandler::connect_update_model_config_info_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_MODEL_CONFIG_INFO_E,
			boost::bind( &ClientMessageHandler::connect_request_model_config_info_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_NRRD_LIST_E,
			boost::bind( &ClientMessageHandler::connect_request_nrrd_list_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_RUN_STAGE_E,
			boost::bind( &ClientMessageHandler::connect_request_run_stage_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_TERMINATE_STAGE_E,
			boost::bind( &ClientMessageHandler::connect_request_terminate_stage_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_VISUALIZE_STAGE_E,
			boost::bind( &ClientMessageHandler::connect_request_visualize_stage_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_UPLOAD_FILE_E,
			boost::bind( &ClientMessageHandler::connect_upload_nrrd_create_project_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_TRANSFER_FILE_E,
			boost::bind( &ClientMessageHandler::connect_transfer_file_fragment_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_DOWNLOAD_FILE_E,
			boost::bind( &ClientMessageHandler::connect_request_download_file_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_CLEAR_ERROR_MESSAGE_E,
			boost::bind( &ClientMessageHandler::connect_clear_error_message_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_LATENCY_E,
			boost::bind( &ClientMessageHandler::connect_latency_handler, this, _1 ) );

		register_callback( MessageTypes::CLIENT_REQUEST_QUERY_PRIMARY_PROJECT_E,
			boost::bind( &ClientMessageHandler::connect_query_primary_project_message_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_CREATE_PRIMARY_PROJECT_UPLOAD_FILE_E,
			boost::bind( &ClientMessageHandler::connect_upload_nrrd_create_primary_project_message_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_DELETE_LABEL_MAP_E,
			boost::bind( &ClientMessageHandler::connect_request_delete_label_map_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_DELETE_MESH_E,
			boost::bind( &ClientMessageHandler::connect_request_delete_mesh_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_JOIN_SESSION_E,
			boost::bind( &ClientMessageHandler::connect_request_join_vis_session_message_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_LEAVE_SESSION_E,
			boost::bind( &ClientMessageHandler::connect_request_leave_vis_session_message_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_REQUEST_RESUME_SESSION_E,
			boost::bind( &ClientMessageHandler::connect_request_resume_vis_session_message_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_CANCEL_FILE_TRANSFER_E,
			boost::bind( &ClientMessageHandler::connect_handle_cancel_file_transfer_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_MODIFY_USER_GROUPS_E,
			boost::bind( &ClientMessageHandler::connect_modify_user_groups_handler, this, _1 ) );
		register_callback( MessageTypes::CLIENT_ADD_USER_GROUP_E,
			boost::bind( &ClientMessageHandler::connect_add_user_group_handler, this, _1 ) );
	}

	void ClientMessageHandler::unregister_handlers()
	{
		unregister_callback( MessageTypes::CLIENT_KEYBOARD_E );
		unregister_callback( MessageTypes::CLIENT_MOUSE_E );
		unregister_callback( MessageTypes::CLIENT_CLIPPING_E );
		unregister_callback( MessageTypes::CLIENT_RESIZE_E );
		unregister_callback( MessageTypes::CLIENT_QUALITY_ADJUST_E );
		unregister_callback( MessageTypes::CLIENT_EXIT_SCIRUN_E );
		unregister_callback( MessageTypes::PING_E );
		unregister_callback( MessageTypes::CLIENT_QUALITY_OF_SERVICE_E );

		unregister_callback( MessageTypes::CLIENT_AUTHENTICATE_CLIENT_E );
		unregister_callback( MessageTypes::CLIENT_ADD_USER_E );
		unregister_callback( MessageTypes::CLIENT_ADD_LABEL_MAP_E );
		unregister_callback( MessageTypes::CLIENT_ADD_MESH_E );
		unregister_callback( MessageTypes::CLIENT_UPDATE_MODEL_CONFIG_INFO_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_MODEL_CONFIG_INFO_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_NRRD_LIST_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_RUN_STAGE_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_TERMINATE_STAGE_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_VISUALIZE_STAGE_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_UPLOAD_FILE_E );
		unregister_callback( MessageTypes::CLIENT_TRANSFER_FILE_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_QUERY_PRIMARY_PROJECT_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_JOIN_SESSION_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_LEAVE_SESSION_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_RESUME_SESSION_E );
		unregister_callback( MessageTypes::CLIENT_CANCEL_FILE_TRANSFER_E );
		unregister_callback( MessageTypes::CLIENT_LATENCY_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_DELETE_LABEL_MAP_E );
		unregister_callback( MessageTypes::CLIENT_REQUEST_DELETE_MESH_E );
		unregister_callback( MessageTypes::CLIENT_MODIFY_USER_GROUPS_E );
		unregister_callback( MessageTypes::CLIENT_ADD_USER_GROUP_E );
	}

	void ClientMessageHandler::connect_keyboard_handler( MessageHandle msg )
	{
		// This function directly accesses the header for more speedy retrieval of our data
		char key_pressed;
        //char input = 0;
		bool ctrl, alt, shift;
        //bool is_valid = true;

		if ( msg->retrieve_message_parameter( ParameterTypes::KEY_PRESSED_E, key_pressed ) &&
			 msg->retrieve_message_parameter( ParameterTypes::CTRL_PRESSED_E, ctrl ) &&
			 msg->retrieve_message_parameter( ParameterTypes::ALT_PRESSED_E, alt ) &&
			 msg->retrieve_message_parameter( ParameterTypes::SHIFT_PRESSED_E, shift ) )
		{
			handle_keyboard( key_pressed, ctrl, alt, shift );
		}
	}

	void ClientMessageHandler::connect_mouse_handler( MessageHandle msg )
	{
		// This function directly accesses the header for more speedy retrieval of our data
		short mouse_x, mouse_y = 0;
		//char input = 0;
		bool mouse_right, mouse_middle, mouse_left, ctrl, alt, shift;
        //bool is_valid = true;
		unsigned long timestamp;

		if ( msg->retrieve_message_parameter( ParameterTypes::MOUSE_X_E, mouse_x ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MOUSE_Y_E, mouse_y ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MOUSE_RIGHT_E, mouse_right ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MOUSE_MIDDLE_E, mouse_middle ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MOUSE_LEFT_E, mouse_left ) &&
			 msg->retrieve_message_parameter( ParameterTypes::CTRL_PRESSED_E, ctrl ) &&
			 msg->retrieve_message_parameter( ParameterTypes::ALT_PRESSED_E, alt ) &&
			 msg->retrieve_message_parameter( ParameterTypes::SHIFT_PRESSED_E, shift ) &&
			 msg->retrieve_message_parameter( ParameterTypes::TIMESTAMP_E, timestamp ) )
		{
			handle_mouse( mouse_x, mouse_y, mouse_left, mouse_middle, mouse_right, ctrl, alt, 
				shift, timestamp );
		}
	}

	void ClientMessageHandler::connect_clipping_handler( MessageHandle msg )
	{
		unsigned short plane;
		bool enable;
		bool show_frame;
		bool reverse_normal;
		float x, y, z, d;

		if ( msg->retrieve_message_parameter( ParameterTypes::CLIPPING_PLANE_E, plane ) &&
			msg->retrieve_message_parameter( ParameterTypes::CLIPPING_ENABLE_PLANE_E, enable ) &&
			msg->retrieve_message_parameter( ParameterTypes::CLIPPING_SHOW_FRAME_E, show_frame ) &&
			msg->retrieve_message_parameter( ParameterTypes::CLIPPING_REVERSE_NORMAL_E,  reverse_normal ) &&
			msg->retrieve_message_parameter( ParameterTypes::CLIPPING_X_E, x ) &&
			msg->retrieve_message_parameter( ParameterTypes::CLIPPING_Y_E, y ) &&
			msg->retrieve_message_parameter( ParameterTypes::CLIPPING_Z_E, z ) &&
			msg->retrieve_message_parameter( ParameterTypes::CLIPPING_D_E, d ) )
		{
			this->handle_clipping( plane, enable, show_frame, reverse_normal, x, y, z, d );
		}
	}


	void ClientMessageHandler::connect_resize_handler( MessageHandle msg )
	{
		unsigned short width, height = 0;
		if ( msg->retrieve_message_parameter( ParameterTypes::WIDTH_E, width ) &&
			 msg->retrieve_message_parameter( ParameterTypes::HEIGHT_E, height ) )
		{
			handle_resize( width, height );
		}
	}

	void ClientMessageHandler::connect_quality_adjust_handler( MessageHandle msg )
	{
		unsigned int quality = 0;
		unsigned short gop = 0;

		
		if ( msg->retrieve_message_parameter( ParameterTypes::GOP_E, gop ) &&
			 msg->retrieve_message_parameter( ParameterTypes::QUALITY_E, quality ) )
		{
			handle_quality_adjust( quality, gop );
		}
	}

	void ClientMessageHandler::connect_exit_scirun_handler( MessageHandle msg )
	{
		handle_exit_scirun();
	}

	void ClientMessageHandler::connect_ping_handler( MessageHandle msg )
	{
		handle_ping();
	}

	void ClientMessageHandler::connect_quality_of_service_handler( MessageHandle msg )
	{
		unsigned short messages_received = 0;
		if ( msg->retrieve_message_parameter( ParameterTypes::MESSAGES_RECEIVED_E, messages_received ) )
		{
			handle_quality_of_service( messages_received );
		}
	}

	void ClientMessageHandler::connect_latency_handler( MessageHandle msg )
	{
		unsigned long timestamp = 0;
		if ( msg->retrieve_message_parameter( ParameterTypes::TIMESTAMP_E, timestamp ) )
		{
			handle_latency( timestamp );
		}
	}

	void ClientMessageHandler::connect_authenticate_client_handler( MessageHandle msg )
	{
		//size_t binary_size = msg->get_binary_size();
		unsigned char* binary_data = msg->get_binary_buf();
		std::vector< char > username, password;
		std::string username_size_string( ( char* )binary_data );
		size_t username_size = boost::lexical_cast< size_t >( username_size_string );
		username.resize( username_size );
		password.resize( msg->get_binary_size() - ( username_size_string.size() + 1 + username_size ) );
		memcpy( ( unsigned char* )&username[0], &binary_data[username_size_string.size() + 1], username_size );
		memcpy( ( unsigned char* )&password[0], &binary_data[username_size_string.size() + 1 + username_size], username_size );
		handle_authenticate_client( username, password );
	}

	void ClientMessageHandler::connect_add_user_handler( MessageHandle msg )
	{
		std::string first_name, last_name = "";
//		size_t binary_size = msg->get_binary_size();
		unsigned char* binary_data = msg->get_binary_buf();
		std::vector< char > username, password;
		std::string username_size_string( ( char* )binary_data );
		size_t username_size = boost::lexical_cast< size_t >( username_size_string );
		username.resize( username_size );
		password.resize( msg->get_binary_size() - ( username_size_string.size() + 1 + username_size ) );
		memcpy( ( unsigned char* )&username[0], &binary_data[username_size_string.size() + 1], username_size );
		memcpy( ( unsigned char* )&password[0], &binary_data[username_size_string.size() + 1 + username_size], username_size );
		if ( msg->retrieve_message_parameter( ParameterTypes::FIRST_NAME_E, first_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::LAST_NAME_E, last_name ) )
		{
			handle_add_user( username, password, first_name, last_name );
		}
	}

	void ClientMessageHandler::connect_add_label_map_message_handler( MessageHandle msg )
	{
		std::string label_map_name, label_map_description, group_name = "";
		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::GROUP_NAME_E, group_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::DESCRIPTION_E, label_map_description ) )
		{
			handle_add_label_map( label_map_name, group_name, label_map_description );
		}
	}

	void ClientMessageHandler::connect_add_mesh_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name = "";

		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) )
		{
			handle_add_mesh( label_map_name, mesh_name );
		}
	}
	
	void ClientMessageHandler::connect_request_model_config_info_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name = "";
		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) )
		{
			handle_request_model_config_info( label_map_name, mesh_name );
		}
	}
	
	void ClientMessageHandler::connect_request_nrrd_list_handler( MessageHandle msg )
	{
		handle_request_nrrd_list();
	}

	inline bool stage_message_helper( MessageHandle msg, std::string& label_map_name,
		std::string& mesh_name,	unsigned short& stage_number, std::string& server_ip_address )
	{
		return msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) &&
			msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number ) &&
			msg->retrieve_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address );
	}
	
	void ClientMessageHandler::connect_request_run_stage_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name, server_ip_address = "";
		unsigned short stage_number = 0;
		if ( stage_message_helper( msg, label_map_name, mesh_name, stage_number, server_ip_address ) )
		{
			handle_request_run_stage( label_map_name, mesh_name, server_ip_address, stage_number );
		}
	}
	
	void ClientMessageHandler::connect_request_terminate_stage_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name = "";
		unsigned short stage_number = 0;
		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number ) )
		{
			handle_request_terminate_stage( label_map_name, mesh_name, stage_number );
		}
	}
	
	void ClientMessageHandler::connect_request_visualize_stage_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name, server_ip_address = "";
		unsigned short stage_number = 0;
		bool is_public_vis_session;

		if( stage_message_helper( msg, label_map_name, mesh_name, stage_number, server_ip_address ) &&
			msg->retrieve_message_parameter(ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session ) )
		{
			handle_request_visualize_stage( label_map_name, mesh_name, stage_number, 
				server_ip_address, is_public_vis_session );
		}
	}
	
	void ClientMessageHandler::connect_update_model_config_info_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name, description, mats, mat_names, mat_radii, refinement_levels, 
			constant_sizing_value, max_sizing_field, num_particle_iters, tetgen_joined_vol_flags = "";
		unsigned short stage_number;
		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) &&
			msg->retrieve_message_parameter( ParameterTypes::DESCRIPTION_E, description ) &&
			 msg->retrieve_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MATS_E, mats ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MAT_NAMES_E, mat_names ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MAT_RADII_E, mat_radii ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_REFINEMENT_LEVELS_E, refinement_levels ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_CONSTANT_SIZING_VALUE_E, constant_sizing_value ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MAX_SIZING_FIELD_E, max_sizing_field ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_NUM_PARTICLE_ITERS_E, num_particle_iters ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_TETGEN_JOINED_VOL_FLAGS_E, tetgen_joined_vol_flags ) )
		{
			handle_update_model_config_info( label_map_name, mesh_name, description, stage_number, mats, 
				mat_names, mat_radii, refinement_levels, constant_sizing_value, max_sizing_field, num_particle_iters, 
				tetgen_joined_vol_flags );
		}
	}
	
	void ClientMessageHandler::connect_upload_nrrd_create_project_handler( MessageHandle msg )
	{
		std::string full_project_name, file_name, mats, mat_names, 
			mat_radii, refinement_levels, constant_sizing_value, max_sizing_field, 
			num_particle_iters, tetgen_joined_vol_flags = "";
		size_t file_size = 0;
		unsigned char* file_data = msg->get_binary_buf();
		if ( msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_PATH_E, full_project_name ) && 
			 msg->retrieve_message_parameter( ParameterTypes::FILE_NAME_E, file_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FILE_SIZE_E, file_size ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MATS_E, mats ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MAT_NAMES_E, mat_names ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MAT_RADII_E, mat_radii ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_REFINEMENT_LEVELS_E, refinement_levels ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_CONSTANT_SIZING_VALUE_E, constant_sizing_value ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_MAX_SIZING_FIELD_E, max_sizing_field ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_NUM_PARTICLE_ITERS_E, num_particle_iters ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MODEL_TETGEN_JOINED_VOL_FLAGS_E, tetgen_joined_vol_flags ) &&
			 file_data != NULL )
		{
			handle_upload_nrrd_create_project( full_project_name, file_name, file_size, file_data, mats, mat_names, mat_radii, 
				refinement_levels, constant_sizing_value, max_sizing_field, num_particle_iters, tetgen_joined_vol_flags );
		}
	}

	void ClientMessageHandler::connect_transfer_file_fragment_handler( MessageHandle msg )
	{
		std::string file_name, label_map_name, mesh_name;
		size_t file_size = 0;
		unsigned short fragment_size = 0;
		unsigned char* file_data = msg->get_binary_buf();
		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) &&
			msg->retrieve_message_parameter( ParameterTypes::FILE_NAME_E, file_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FILE_SIZE_E, file_size ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FRAGMENT_SIZE_E, fragment_size ) &&
			 file_data != NULL )
		{
			handle_transfer_file_fragment( label_map_name, mesh_name, file_name, file_size, file_data, fragment_size );
		}
	}

	void ClientMessageHandler::connect_request_download_file_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name, file_path = "";
		size_t bytes_transferred = 0;
		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FILE_NAME_E, file_path ) &&
			 msg->retrieve_message_parameter( ParameterTypes::BYTES_TRANSFERRED_E, bytes_transferred ) )
		{
			handle_request_file_fragment( label_map_name, mesh_name, file_path, bytes_transferred );
		}
	}

	void ClientMessageHandler::connect_clear_error_message_handler( MessageHandle msg )
	{
		std::string full_project_name = "";
		if ( msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_PATH_E, full_project_name ) )
		{
			handle_clear_error_message( full_project_name );
		}
	}

	void ClientMessageHandler::connect_query_primary_project_message_handler( MessageHandle msg )
	{
		std::string full_project_name = "";
		if ( msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_PATH_E, full_project_name ) )
		{
			handle_query_primary_project_message( full_project_name );
		}
		
	}

	void ClientMessageHandler::connect_upload_nrrd_create_primary_project_message_handler( MessageHandle msg )
	{
		std::string owner = "";
		std::string group = "";
		std::string full_project_name = "";
		std::string file_name = "";
		size_t file_size = 0;

		unsigned char* file_data = msg->get_binary_buf();
		if ( msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, owner ) && 
			 msg->retrieve_message_parameter( ParameterTypes::USER_GROUPS_E, group ) && 
			 msg->retrieve_message_parameter( ParameterTypes::FULL_PROJECT_PATH_E, full_project_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FILE_NAME_E, file_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::FILE_SIZE_E, file_size ) &&
			 file_data != NULL )
		{
			handle_upload_nrrd_create_primary_project( owner, group, full_project_name,
				file_name, file_size, file_data );
		}
	}

	void ClientMessageHandler::connect_request_delete_mesh_handler( MessageHandle msg )
	{
		std::string label_map_name, mesh_name = "";

		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) &&
			 msg->retrieve_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name ) )
		{
			handle_request_delete_mesh_message( label_map_name, mesh_name );
		}
	}

	void ClientMessageHandler::connect_request_delete_label_map_handler( MessageHandle msg )
	{
		std::string label_map_name = "";
		
		if ( msg->retrieve_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name ) )
		{
			handle_request_delete_label_map_message( label_map_name );
		}
	}

	void ClientMessageHandler::connect_request_join_vis_session_message_handler( MessageHandle msg )
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
			handle_request_join_vis_session_message( 
				full_project_name, server_ip_address,
				stage_number, port, member_name, is_public_vis_session,
				group_id, sub_project_id );
		}

	}

	void ClientMessageHandler::connect_request_resume_vis_session_message_handler( MessageHandle msg )
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
			handle_request_resume_vis_session_message( 
				full_project_name, server_ip_address,
				stage_number, port, owner_name, member_name, is_public_vis_session,
				group_id, sub_project_id );
		}

	}

	void ClientMessageHandler::connect_request_leave_vis_session_message_handler( MessageHandle msg )
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
			handle_request_leave_vis_session_message( 
				full_project_name, server_ip_address,
				stage_number, port, member_name, is_public_vis_session,
				group_id, sub_project_id );
		}


	}

	void ClientMessageHandler::connect_handle_cancel_file_transfer_handler( MessageHandle msg )
	{
		handle_cancel_file_transfer_message();
	}

	void ClientMessageHandler::connect_modify_user_groups_handler( MessageHandle msg )
	{
		std::string username, groups_string = "";
		std::vector< std::string > groups;
		if ( msg->retrieve_message_parameter( ParameterTypes::USERNAME_E, username ) &&
			msg->retrieve_message_parameter( ParameterTypes::USER_GROUPS_E, groups_string ) )
		{
			boost::algorithm::split( groups, groups_string, boost::is_any_of( "," ) );
			handle_modify_user_groups( username, groups );
		}
	}

	void ClientMessageHandler::connect_add_user_group_handler( MessageHandle msg )
	{
		std::string group_name = "";
		if ( msg->retrieve_message_parameter( ParameterTypes::GROUP_NAME_E, group_name ) )
		{
			handle_add_user_group( group_name );
		}
	}

}// end namespace Message