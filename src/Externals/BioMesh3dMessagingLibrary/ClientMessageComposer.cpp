#include "ClientMessageComposer.h"
#include "Message.h"

namespace BioMesh3d
{
	ClientMessageComposer::ClientMessageComposer()
	{
	}

	ClientMessageComposer::~ClientMessageComposer()
	{
		//delete msg;
	}

	MessageHandle ClientMessageComposer::compose_keyboard_message( const char key_pressed, 
		bool ctrl, bool alt, bool shift )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_KEYBOARD_E );
		msg->add_message_parameter( ParameterTypes::KEY_PRESSED_E, key_pressed );
		msg->add_message_parameter( ParameterTypes::CTRL_PRESSED_E, ctrl );
		msg->add_message_parameter( ParameterTypes::ALT_PRESSED_E, alt );
		msg->add_message_parameter( ParameterTypes::SHIFT_PRESSED_E, shift );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_mouse_message( short _x,
		short _y, bool left_click, bool middle_click, bool right_click,
		bool ctrl, bool alt, bool shift, unsigned long timestamp )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_MOUSE_E );
		msg->add_message_parameter( ParameterTypes::MOUSE_X_E, _x );
		msg->add_message_parameter( ParameterTypes::MOUSE_Y_E, _y );
		msg->add_message_parameter( ParameterTypes::MOUSE_LEFT_E, left_click );
		msg->add_message_parameter( ParameterTypes::MOUSE_MIDDLE_E, middle_click );
		msg->add_message_parameter( ParameterTypes::MOUSE_RIGHT_E, right_click );
		msg->add_message_parameter( ParameterTypes::CTRL_PRESSED_E, ctrl );
		msg->add_message_parameter( ParameterTypes::ALT_PRESSED_E, alt );
		msg->add_message_parameter( ParameterTypes::SHIFT_PRESSED_E, shift );
		msg->add_message_parameter( ParameterTypes::TIMESTAMP_E, timestamp );
		return msg;
	}

	BioMesh3d::MessageHandle ClientMessageComposer::compose_clipping_message( 
		unsigned short plane, bool enable, bool show_frame, 
		bool reverse_normal, float x, float y, float z, float d )
	{
		MessageHandle msg = MessageHandle( new Message );
		msg->set_message_type( MessageTypes::CLIENT_CLIPPING_E );
		msg->add_message_parameter( ParameterTypes::CLIPPING_PLANE_E, plane );
		msg->add_message_parameter( ParameterTypes::CLIPPING_ENABLE_PLANE_E, enable );
		msg->add_message_parameter( ParameterTypes::CLIPPING_SHOW_FRAME_E, show_frame );
		msg->add_message_parameter( ParameterTypes::CLIPPING_REVERSE_NORMAL_E,  reverse_normal );
		msg->add_message_parameter( ParameterTypes::CLIPPING_X_E, x );
		msg->add_message_parameter( ParameterTypes::CLIPPING_Y_E, y );
		msg->add_message_parameter( ParameterTypes::CLIPPING_Z_E, z );
		msg->add_message_parameter( ParameterTypes::CLIPPING_D_E, d );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_resize_message( const unsigned short width, 
		const unsigned short height )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_RESIZE_E );
		msg->add_message_parameter( ParameterTypes::WIDTH_E, width );
		msg->add_message_parameter( ParameterTypes::HEIGHT_E, height );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_quality_adjust_message( const unsigned int quality,
		const unsigned short gop )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_QUALITY_ADJUST_E );
		msg->add_message_parameter( ParameterTypes::QUALITY_E, quality );
		msg->add_message_parameter( ParameterTypes::GOP_E, gop );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_quality_of_service_message( 
		unsigned short messages_received )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_QUALITY_OF_SERVICE_E );
		msg->add_message_parameter( ParameterTypes::MESSAGES_RECEIVED_E, 
			messages_received );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_latency_message( unsigned long timestamp )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_LATENCY_E );
		msg->add_message_parameter( ParameterTypes::TIMESTAMP_E, timestamp );
		return msg;
	}




	MessageHandle ClientMessageComposer::compose_authenticate_client_message( const std::vector< char >& username, const std::vector< char >& password )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_AUTHENTICATE_CLIENT_E );
		std::vector< char > binary_data;
		std::string username_size_string = boost::lexical_cast< std::string >( username.size() );
		binary_data.resize( username_size_string.size() + username.size() + password.size() + 1 );
		memcpy( ( unsigned char* )&binary_data[0], username_size_string.c_str(), username.size() );
		memcpy( ( unsigned char* )&binary_data[username_size_string.size() + 1], ( unsigned char* )&username[0], username.size() );
		memcpy( ( unsigned char* )&binary_data[username_size_string.size() + username.size() + 1], ( unsigned char* )&password[0], password.size() );
		msg->copy_binary_data( ( unsigned char* )&binary_data[0], binary_data.size() );
		return msg;
	}

	BioMesh3d::MessageHandle ClientMessageComposer::compose_add_user_message( const std::vector< char >& username, const std::vector< char >& password 
		, const std::string& first_name, const std::string& last_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_ADD_USER_E );
		std::vector< char > binary_data;
		std::string username_size_string = boost::lexical_cast< std::string >( username.size() );
		binary_data.resize( username_size_string.size() + username.size() + password.size() + 1 );
		memcpy( ( unsigned char* )&binary_data[0], username_size_string.c_str(), username.size() );
		memcpy( ( unsigned char* )&binary_data[username_size_string.size() + 1], ( unsigned char* )&username[0], username.size() );
		memcpy( ( unsigned char* )&binary_data[username_size_string.size() + username.size() + 1], ( unsigned char* )&password[0], password.size() );
		msg->copy_binary_data( ( unsigned char* )&binary_data[0], binary_data.size() );
		msg->add_message_parameter( ParameterTypes::FIRST_NAME_E, first_name );
		msg->add_message_parameter( ParameterTypes::LAST_NAME_E, last_name );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_add_label_map_message( 
		const std::string& primary_project_name, 
		const std::string& group_name,
		const std::string& primary_project_description )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_ADD_LABEL_MAP_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, primary_project_name );
		msg->add_message_parameter( ParameterTypes::GROUP_NAME_E, group_name );
		msg->add_message_parameter( ParameterTypes::DESCRIPTION_E,
			primary_project_description );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_add_mesh_message( 
		const std::string& label_map_name, const std::string& mesh_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_ADD_MESH_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		return msg;
	}

	BioMesh3d::MessageHandle ClientMessageComposer::compose_request_model_config_info_message( 
		const std::string& label_map__name, const std::string& mesh_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_MODEL_CONFIG_INFO_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map__name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_request_nrrd_list_message()
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_NRRD_LIST_E );
		return msg;
	}

	BioMesh3d::MessageHandle ClientMessageComposer::compose_request_run_stage_message( 
		const std::string& label_map_name, const std::string& mesh_name, 
		const std::string& server_ip_address, const unsigned short stage_number )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_RUN_STAGE_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address );
		return msg;
	}

	BioMesh3d::MessageHandle ClientMessageComposer::compose_request_terminate_stage_message( const std::string& label_map_name, const std::string& mesh_name, unsigned short stage_number )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_TERMINATE_STAGE_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		return msg;
	}

	BioMesh3d::MessageHandle ClientMessageComposer::compose_request_visualize_stage_message( const std::string& label_map_name, const std::string& mesh_name, unsigned short stage_number, const std::string server_ip_address, bool is_public_vis_session )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_VISUALIZE_STAGE_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address );
		msg->add_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session );
		return msg;
	}

	BioMesh3d::MessageHandle ClientMessageComposer::compose_update_model_config_info_message( 
		const std::string label_map_name, const std::string mesh_name,
		const std::string& description, const unsigned short stage_number, 
		const std::string mats, const std::string mat_names, const std::string mat_radii, 
		const std::string refinement_levels, const std::string constant_sizing_value, 
		const std::string max_sizing_field, 
		const std::string num_particle_iters, const std::string tetgen_joined_vol_flags )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_UPDATE_MODEL_CONFIG_INFO_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::DESCRIPTION_E, description );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
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

	MessageHandle ClientMessageComposer::compose_upload_nrrd_create_project_message( const std::string& full_project_path, 
		const std::string& file_name, size_t file_size, unsigned char* file_data, const std::string& mats, const std::string& mat_names,
		const std::string& mat_radii, const std::string& refinement_levels,
		const std::string& constant_sizing_value,
		const std::string& max_sizing_field, const std::string& num_particle_iters, 
		const std::string& tetgen_joined_vol_flags )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_UPLOAD_FILE_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_PATH_E, full_project_path );
		msg->add_message_parameter( ParameterTypes::FILE_NAME_E, file_name );
		msg->add_message_parameter( ParameterTypes::FILE_SIZE_E, file_size );
		msg->add_message_parameter( ParameterTypes::MODEL_MATS_E, mats );
		msg->add_message_parameter( ParameterTypes::MODEL_MAT_NAMES_E, mat_names );
		msg->add_message_parameter( ParameterTypes::MODEL_MAT_RADII_E, mat_radii );
		msg->add_message_parameter( ParameterTypes::MODEL_REFINEMENT_LEVELS_E, constant_sizing_value );
		msg->add_message_parameter( ParameterTypes::MODEL_CONSTANT_SIZING_VALUE_E, max_sizing_field );		
		msg->add_message_parameter( ParameterTypes::MODEL_MAX_SIZING_FIELD_E, max_sizing_field );
		msg->add_message_parameter( ParameterTypes::MODEL_NUM_PARTICLE_ITERS_E, num_particle_iters );
		msg->add_message_parameter( ParameterTypes::MODEL_TETGEN_JOINED_VOL_FLAGS_E, tetgen_joined_vol_flags );
		msg->copy_binary_data( file_data, file_size );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_upload_nrrd_create_primary_project_message(
		const std::string& owner,
		const std::string& group,
		const std::string& full_project_path, 
		const std::string& file_name, size_t file_size,
		unsigned char* file_data )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type(
			MessageTypes::CLIENT_REQUEST_CREATE_PRIMARY_PROJECT_UPLOAD_FILE_E );
		
		msg->add_message_parameter(
			ParameterTypes::USERNAME_E, owner );

		msg->add_message_parameter(
			ParameterTypes::USER_GROUPS_E, group );

		msg->add_message_parameter( 
			ParameterTypes::FULL_PROJECT_PATH_E, full_project_path );
		msg->add_message_parameter( ParameterTypes::FILE_NAME_E, file_name );
		msg->add_message_parameter( ParameterTypes::FILE_SIZE_E, file_size );

		msg->copy_binary_data( file_data, file_size );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_query_primary_project_message( const std::string& primary_project_name)
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_QUERY_PRIMARY_PROJECT_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_PATH_E, primary_project_name );

		return msg;
	}

	MessageHandle ClientMessageComposer::compose_request_delete_mesh_message( 
		const std::string& label_map_name, const std::string& mesh_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_DELETE_MESH_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );

		return msg;
	}

	MessageHandle ClientMessageComposer::compose_request_delete_label_map_message( 
		const std::string& label_map_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_DELETE_LABEL_MAP_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );

		return msg;
	}

	BioMesh3d::MessageHandle ClientMessageComposer::compose_transfer_file_fragment_message( 
		const std::string& label_map_name, const std::string& mesh_name, 
		const std::string& file_name, size_t file_size, unsigned char* file_data, 
		unsigned short fragment_size )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_TRANSFER_FILE_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::FILE_NAME_E, file_name );
		msg->add_message_parameter( ParameterTypes::FILE_SIZE_E, file_size );
		msg->add_message_parameter( ParameterTypes::FRAGMENT_SIZE_E, fragment_size );
		msg->copy_binary_data( file_data, fragment_size );
		return msg;
	}

	BioMesh3d::MessageHandle ClientMessageComposer::compose_request_file_message( const std::string& label_map_name, 
		const std::string& mesh_name, const std::string& file_name, size_t bytes_transferred )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_DOWNLOAD_FILE_E );
		msg->add_message_parameter( ParameterTypes::LABEL_MAP_NAME_E, label_map_name );
		msg->add_message_parameter( ParameterTypes::MESH_NAME_E, mesh_name );
		msg->add_message_parameter( ParameterTypes::FILE_NAME_E, file_name );
		msg->add_message_parameter( ParameterTypes::BYTES_TRANSFERRED_E, bytes_transferred );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_clear_project_error( const std::string& full_project_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_CLEAR_ERROR_MESSAGE_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_PATH_E, full_project_name );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_exit_scirun_message()
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_EXIT_SCIRUN_E );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_ping_message()
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::PING_E );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_leave_vis_session_message( 
		const std::string& full_project_name, 
		const std::string& server_ip_address, 
		const unsigned short stage_number,
		const unsigned short port,
		const std::string& member_name,
		bool is_public_vis_session,
		unsigned short group_id,
		unsigned int sub_project_id)
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_LEAVE_SESSION_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::PORT_E, port );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, member_name );
		msg->add_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session );
		msg->add_message_parameter( ParameterTypes::GROUP_NAME_E, group_id );
		msg->add_message_parameter( ParameterTypes::SUB_PROJECT_ID_E, sub_project_id );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_join_vis_session_message( 
		const std::string& full_project_name, 
		const std::string& server_ip_address, 
		const unsigned short stage_number,
		const unsigned short port,
		const std::string& member_name,
		bool is_public_vis_session,
		unsigned short group_id,
		unsigned int sub_project_id)
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_JOIN_SESSION_E );
		msg->add_message_parameter( ParameterTypes::FULL_PROJECT_NAME_E, full_project_name );
		msg->add_message_parameter( ParameterTypes::SERVER_IP_ADDRESS_E, server_ip_address );
		msg->add_message_parameter( ParameterTypes::STAGE_NUMBER_E, stage_number );
		msg->add_message_parameter( ParameterTypes::PORT_E, port );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, member_name );
		msg->add_message_parameter( ParameterTypes::IS_PUBLIC_VIS_SESSION_E, is_public_vis_session );
		msg->add_message_parameter( ParameterTypes::GROUP_NAME_E, group_id );
		msg->add_message_parameter( ParameterTypes::SUB_PROJECT_ID_E, sub_project_id );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_resume_vis_session_message( 
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
		msg->set_message_type( MessageTypes::CLIENT_REQUEST_RESUME_SESSION_E );
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

	MessageHandle ClientMessageComposer::compose_cancel_file_transfer_message()
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_CANCEL_FILE_TRANSFER_E );
		return msg;
	}

	MessageHandle ClientMessageComposer::compose_modify_user_groups( const std::string& username, const std::vector< std::string >& groups )
	{
		std::string groups_string = "";
		for ( size_t i = 0; i < groups.size(); ++i )
		{
			std::string comma = ( i < groups.size() - 1 ) ? "," : "";
			groups_string += groups[i] + comma;
		}
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_MODIFY_USER_GROUPS_E );
		msg->add_message_parameter( ParameterTypes::USERNAME_E, username );
		msg->add_message_parameter( ParameterTypes::USER_GROUPS_E, groups_string );
		return msg;
	}

	BioMesh3d::MessageHandle ClientMessageComposer::compose_add_user_group( const std::string& group_name )
	{
		MessageHandle msg = MessageHandle( new Message() );
		msg->set_message_type( MessageTypes::CLIENT_ADD_USER_GROUP_E );
		msg->add_message_parameter( ParameterTypes::GROUP_NAME_E, group_name );
		return msg;
	}


}// end namespace Message