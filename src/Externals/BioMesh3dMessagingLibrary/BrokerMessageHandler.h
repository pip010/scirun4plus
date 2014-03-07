#ifndef SERVER_MESSAGE_HANDLER_H
#define SERVER_MESSAGE_HANDLER_H

#include "MessageHandler.h"

namespace BioMesh3d
{

	class BrokerMessageHandler : public MessageHandler
	{
	public:
		BrokerMessageHandler();
		~BrokerMessageHandler();

	protected:
		virtual void handle_is_client_authenticated( bool authenticated, const std::string& validation_message ) { std::cout << "handle_is_client_authenticated not connected" << std::endl; }
		virtual void handle_initiate_file_upload( const std::string& full_project_name ) { std::cout << "handle_initiate_file_upload not connected" << std::endl; }
		virtual void handle_mesh_status( const std::string& label_map_name, const std::string& mesh_name, 
			const StageStates::StageState stage_state, const unsigned short stage_number, 
			const unsigned long date_processed, const std::string& error_message, 
			const std::string& mesh_description )
		{ std::cout << "handle_project_info not connected" << std::endl; }
		virtual void handle_server_list( const std::vector< std::string >& server_ip_address_list, const std::vector< std::string >& server_hostname_list ) { std::cout << "handle_server_list not connected" << std::endl; }
		virtual void handle_start_file_download( const std::string& full_project_name, 
			const std::string& file_name, unsigned int file_size ) { std::cout << "handle_start_file_download not connected" << std::endl; }
		virtual void handle_file_transfer_complete( const std::string& full_project_name, unsigned char* file_data, size_t file_size ) { std::cout << "handle_file_transfer_complete not connected" << std::endl; }
		virtual void handle_transfer_file_fragment( const std::string& label_map_name, const std::string& mesh_name, const std::string& file_name, 
			size_t file_size, unsigned char* file_data, unsigned short fragment_size ) { std::cout << "handle_transfer_file_fragment not connected" << std::endl; }
		virtual void handle_model_config_info( const std::string& label_map_name, 
			const std::string& mesh_name, const std::string& description, const std::string& mats, 
			const std::string& mat_names, const std::string& mat_radii, const std::string& refinement_levels, 
			const std::string& constant_sizing_value, const std::string& max_sizing_field, 
			const std::string& num_particle_iters, 
			const std::string& tetgen_joined_vol_flags ) { std::cout << "handle_model_config_info not connected" << std::endl; }

		virtual void handle_request_run_stage( 
			const std::string& label_map_name,
			const std::string& mesh_name, 
			unsigned short start_stage_number,
			unsigned short end_stage_number ) { std::cout << "handle_request_run_stage not connected" << std::endl; }

		virtual void handle_request_stop_stage( 
			const std::string& label_map_name,
			const std::string& mesh_name ) { std::cout << "handle_request_stop_stage not connected" << std::endl; }

		virtual void handle_request_terminate_stage(
			const std::string& label_map_name,
			const std::string& mesh_name, 
			unsigned short stage_number ) { std::cout << "handle_request_terminate_stage not connected" << std::endl; }

		virtual void handle_request_visualize_stage_message( 
			const std::string& label_map_name, 
			const std::string& mesh_name, 
			unsigned short stage_number,
			const std::string& client_ip_address,
			const std::string& user_name, bool is_public_vis_session ) { std::cout << "handle_request_visualize_stage_message not connected" << std::endl; }

		virtual void handle_visualization_started_message( 
			const std::string& label_map_name, 
			const std::string& mesh_name,
			unsigned short group_id,
			const std::string& server_ip_address, unsigned short stage_number,
			unsigned short port, const std::string& user_name,
			bool is_public_vis_session ) { std::cout << "handle_visualization_started_message not connected" << std::endl; }

		virtual void handle_answer_primary_project_message ( bool is_find ){ std::cout << "handle_answer_primary_project_message not connected" << std::endl; }
		virtual void handle_answer_user_groups_message ( const std::vector< std::string >& groups ){ std::cout << "handle_answer_groups_message not connected" << std::endl; }
		virtual void handle_delete_label_map_message( const std::string& label_map_name ){ std::cout << "handle_answer_delete_project_message not connected" << std::endl; }
		virtual void handle_delete_mesh_message( const std::string& label_map_name, const std::string& mesh_name ){ std::cout << "handle_answer_delete_project_message not connected" << std::endl; }
		virtual void handle_label_map_message( 
			const std::string& label_map_name,
			unsigned long date_processed, const std::string& label_map_description ){ std::cout << "handle_answer_create_primary_project_message not connected" << std::endl; }
		virtual void handle_add_mesh_message( const std::string& label_map_name,
			const std::string& mesh_name, const unsigned short stage_number, 
			const StageStates::StageState stage_state, const unsigned long date_processed, 
			const std::string& error_message, const std::string& mesh_description ){ std::cout << "handle_answer_create_subproject_message not connected" << std::endl; }
		virtual void handle_set_project_path_message( const std::string& project_path ) { std::cout << "handle_set_project_path_message not connected" << std::endl; }
		virtual void handle_client_login_fail_message( const std::string& error_description){ std::cout << "handle_client_login_fail_message not connected" << std::endl; }
		virtual void handle_vis_notify_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address,
			unsigned short stage_number, unsigned short port,const std::string& owner_name,
			bool is_public_vis_session, unsigned short group_id,
			const std::string& members, unsigned int sub_project_id ){ std::cout << "handle_vis_notify_message not connected" << std::endl; }

		virtual void handle_notify_leave_vis_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address,
			unsigned short stage_number, unsigned short port, 
			const std::string& owner_name,
			const std::string& member_name,
			bool is_public_vis_session, unsigned short group_id,
			unsigned int sub_project_id ){ std::cout << "handle_notify_leave_vis_message not connected" << std::endl; }

		virtual void handle_notify_join_vis_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address,
			unsigned short stage_number, unsigned short port, 
			const std::string& member_name,
			bool is_public_vis_session, unsigned short group_id,
			unsigned int sub_project_id ){ std::cout << "handle_notify_join_vis_message not connected" << std::endl; }

		virtual void handle_visualization_resumed_message( 
			const std::string& full_project_name, 
			unsigned short group_id,
			const std::string& server_ip_address, unsigned short stage_number,
			unsigned short port, const std::string& user_name,
			bool is_public_vis_session ) { std::cout << "handle_visualization_resumed_message not connected" << std::endl; }

		virtual void handle_ping_message() { std::cout << "handle_ping_message not connected" << std::endl; }
		virtual void handle_alert_message( const std::string& alert_message ) { std::cout << "handle_alert_message not connected" << std::endl; }
		virtual void handle_confirm_file_transfer_cancel_message() { std::cout << "handle_confirm_file_transfer_cancel_message not connected" << std::endl; }
		virtual void handle_rsa_pem_public_key_message( const std::string& rsa_pem_public_key ) { std::cout << "handle_rsa_pem_public_key_message not connected" << std::endl; }
		virtual void handle_users_info( const std::vector< std::string >& usernames,
			const std::vector< std::string >& first_names, const std::vector< std::string >& last_names,
			const std::vector< std::vector< std::string > >& groups ) { std::cout << "handle_users_info not connected" << std::endl; }

		private:
		void register_handlers();
		void unregister_handlers();

		void connect_is_client_authenticated_handler( MessageHandle msg );
		void connect_initiate_file_upload_handler( MessageHandle msg );
		void connect_mesh_handler( MessageHandle msg );
		void connect_server_list_handler( MessageHandle msg );
		void connect_start_file_download_handler( MessageHandle msg );
		void connect_file_transfer_complete_handler( MessageHandle msg );
		void connect_transfer_file_fragment_handler( MessageHandle msg );
		void connect_model_config_info_handler( MessageHandle msg );
		void connect_delete_label_map_handler ( MessageHandle msg );
		void connect_delete_mesh_handler ( MessageHandle msg );
		
		void connect_request_run_stage_message_handler( MessageHandle msg );
		void connect_request_stop_stage_message_handler( MessageHandle msg );
		void connect_request_terminate_stage_message_handler( MessageHandle msg );
		void connect_request_visualize_stage_message_handler( MessageHandle msg );
		void connect_visualization_started_message_handler( MessageHandle msg );
		void connect_answer_primary_project_message_handler( MessageHandle msg);
		void connect_answer_user_groups_message_handler( MessageHandle msg);
		void connect_add_mesh_handler( MessageHandle msg );
		void connect_label_map_handler( MessageHandle msg );
		void connect_set_project_path_handler( MessageHandle msg );
		void connect_login_fail_handler( MessageHandle msg );
		void connect_vis_notify_handler( MessageHandle msg );
		void connect_notify_leave_vis_handler( MessageHandle msg );
		void connect_notify_join_vis_handler( MessageHandle msg );
		void connect_visualization_resumed_message_handler( MessageHandle msg );
		void connect_ping_handler( MessageHandle msg );
		void connect_alert_handler( MessageHandle msg );
		void connect_confirm_file_transfer_cancel_handler( MessageHandle msg );
		void connect_rsa_pem_public_key_handler( MessageHandle msg );
		void connect_users_info_handler( MessageHandle msg );
	};
}

#endif