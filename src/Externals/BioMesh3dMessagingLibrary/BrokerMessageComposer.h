#ifndef SERVER_MANAGER_MESSAGE_COMPOSER
#define SERVER_MANAGER_MESSAGE_COMPOSER

#include <vector>

#include "Message.h"

namespace BioMesh3d
{

	class BrokerMessageComposer
	{
	public:
		BrokerMessageComposer();
		~BrokerMessageComposer();
	protected:
		MessageHandle compose_is_client_authenticated_message( bool authenticated, 
			const std::string& validation_message="" );
		MessageHandle compose_client_login_fail_message( const std::string& error_description );
		MessageHandle compose_initiate_file_upload_message( const std::string& full_project_name );
		MessageHandle compose_mesh_status( const std::string& label_map_name, 
			const std::string& mesh_name, const StageStates::StageState stage_state, 
			const unsigned short stage_number, const unsigned long date_processed, 
			const std::string& error_message, const std::string& mesh_description );

		MessageHandle compose_server_list_message( const std::vector< std::string >& server_ip_address_list, 
			const std::vector< std::string >& server_hostname_list );
		MessageHandle compose_start_file_download_message( const std::string& full_project_name, 
			const std::string& file_name, size_t file_size );
		MessageHandle compose_nrrd_file_names_message( std::vector< std::string > nrrd_file_names );
		MessageHandle compose_file_transfer_complete_message( const std::string& full_project_name,
			unsigned char* file_data, size_t file_size );
		MessageHandle compose_transfer_file_fragment_message( const std::string& label_map_name, 
			const std::string& mesh_name, const std::string& file_name, size_t file_size, 
			unsigned char* file_data, unsigned short fragment_size );
		MessageHandle compose_model_config_info_message( const std::string& label_map_name, 
			const std::string& mesh_name, const std::string& description, const std::string& mats, 
			const std::string& mat_names, const std::string& mat_radii, 
			const std::string& refinement_levels, 
			const std::string& constant_sizing_value, const std::string& max_sizing_field, 
			const std::string& num_particle_iters, const std::string& tetgen_joined_vol_flags );

		MessageHandle compose_request_run_stage_message( const std::string& label_map_name, 
			const std::string& mesh_name, unsigned short start_stage_number, 
			unsigned short end_stage_number );

		MessageHandle compose_request_stop_stage_message( const std::string& label_map_name, 
			const std::string& mesh_name );

		MessageHandle compose_request_terminate_stage_message( 
			const std::string& label_map_name,
			const std::string& mesh_name, 
			unsigned short stage_number );

		MessageHandle compose_request_visualize_stage_message( 
			const std::string& project_path,
			const std::string& full_project_name,
			unsigned short stage_number,
			const std::string& client_ip_address,
			const std::string& user_name, bool is_public_vis_session );

		MessageHandle compose_visualization_started_message( const std::string& label_map_name, 
			const std::string& mesh_name, unsigned short group_id,
			const std::string& server_ip_address, const unsigned short stage_number, 
			const unsigned short port, const std::string& user_name, bool is_public_vis_session );

		MessageHandle compose_visualization_resumed_message( 
			const std::string& full_project_name, 
			unsigned short group_id,
			const std::string& server_ip_address, 
			const unsigned short stage_number,
			const unsigned short port,
			const std::string& user_name,
			bool is_public_vis_session );

		MessageHandle compose_query_primary_project_message ( const bool is_find );
		MessageHandle compose_delete_label_map_message( const std::string& label_map_name );
		MessageHandle compose_delete_mesh_message( const std::string& label_map_name, 
			const std::string& mesh_name );
		MessageHandle compose_answer_update_subproject_message( const bool is_update,
			std::string description );
		MessageHandle compose_add_mesh_message( const std::string& label_map_name,
			const std::string& mesh_name, const unsigned short stage_number, 
			const StageStates::StageState stage_state, const unsigned long date_processed, 
			const std::string& error_message, const std::string& mesh_description );

		MessageHandle compose_label_map_message( 
			const std::string& label_map_name, 
			const unsigned long date_processed,
			const std::string& label_map_description );

		MessageHandle compose_user_groups_message( const std::vector< std::string >& groups);
		MessageHandle compose_set_project_path_message( const std::string& project_path );

		MessageHandle compose_notify_vis_session_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address, 
			const unsigned short stage_number,
			const unsigned short port,
			const std::string& owner_name,
			bool is_public_vis_session,
			unsigned short group_id,
			const std::string& members,
			unsigned int sub_project_id = 0);

		MessageHandle compose_notify_leave_vis_session_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address, 
			const unsigned short stage_number,
			const unsigned short port,
			const std::string& owner_name,
			const std::string& member_name,
			bool is_public_vis_session,
			unsigned short group_id,
			unsigned int sub_project_id);

		MessageHandle compose_notify_join_vis_session_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address, 
			const unsigned short stage_number,
			const unsigned short port,
			const std::string& member_name,
			bool is_public_vis_session,
			unsigned short group_id,
			unsigned int sub_project_id);

		MessageHandle compose_ping_message();
		MessageHandle compose_alert_message( const std::string& alert_string );
		MessageHandle compose_confirm_file_transfer_cancel_message();
		MessageHandle compose_rsa_pem_public_key_message( const std::string& rsa_pem_public_key_string );

		MessageHandle compose_noop_message();
		MessageHandle compose_users_info( const std::vector< std::string >& usernames,
			const std::vector< std::string >& first_names, const std::vector< std::string >& last_names,
			const std::vector< std::vector< std::string > >& groups );

	private:
	};
}

#endif