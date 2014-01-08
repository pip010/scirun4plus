#ifndef CLIENT_MESSAGE_H
#define CLIENT_MESSAGE_H

#include "MessageHandler.h"

namespace BioMesh3d
{

	class ClientMessageHandler : public MessageHandler
	{
	public:
		ClientMessageHandler();
		~ClientMessageHandler();

	protected:
		virtual void handle_keyboard( char key_pressed, 
			bool ctrl, bool alt, bool shift ) { std::cout << "handle_keyboard not connected" << std::endl; }
		virtual void handle_mouse( short _x,
			short _y, bool left_click, bool middle_click, bool right_click,
			bool ctrl, bool alt, bool shift, unsigned long timestamp ) { std::cout << "handle_mouse not connected" << std::endl; }
		
		virtual void handle_clipping( unsigned short plane, bool enable, 
			bool show_frame, bool reverse_normal, float x, float y, float z, float d )
		{ 
			std::cout << "handle_clipping not connected" << std::endl; 
		}
		
		virtual void handle_resize( unsigned short width, 
			unsigned short height ) { std::cout << "handle_resize not connected" << std::endl; }
		virtual void handle_quality_adjust( unsigned int quality,
			unsigned short gop ) { std::cout << "handle_quality_adjust not connected" << std::endl; }
		virtual void handle_exit_scirun() { std::cout << "handle_exit_scirun not connected" << std::endl; }
		virtual void handle_ping() { std::cout << "handle_ping not connected" << std::endl; }
		virtual void handle_quality_of_service( unsigned short messages_received ) { std::cout << "handle_quality_of_service not connected" << std::endl; }
		virtual void handle_latency( unsigned long timestamp ) { std::cout << "handle_latency not connected" << std::endl; }


		virtual void handle_authenticate_client( const std::vector< char >& username, 
			const std::vector< char >& password ) { std::cout << "handle_authenticate_client not connected" << std::endl; }
		virtual void handle_add_user( const std::vector< char >& username, 
			const std::vector< char >& password, const std::string& first_name, const std::string& last_name ) { std::cout << "handle_add_user not connected" << std::endl; }
		virtual void handle_add_label_map( const std::string& label_map_name, 
			const std::string& group_name,
			const std::string& label_map_description ) { std::cout << "handle_add_primary_project not connected" << std::endl; }
		virtual void handle_add_mesh( const std::string& label_map_name, const std::string& mesh_name )
		{ std::cout << "handle_add_mesh not connected" << std::endl; }
		virtual void handle_request_model_config_info( const std::string& label_map_name, const std::string& mesh_name ) { std::cout << "handle_request_model_config_info not connected" << std::endl; }
		virtual void handle_request_nrrd_list() { std::cout << "handle_request_nrrd_list not connected" << std::endl; }
		virtual void handle_request_run_stage( const std::string& label_map_name, const std::string& mesh_name, const std::string& server_ip_address, unsigned short stage_number ) { std::cout << "handle_request_run_stage not connected" << std::endl; }
		virtual void handle_request_terminate_stage( const std::string& label_map_name, const std::string& mesh_name, unsigned short stage_number ) { std::cout << "handle_request_terminate_stage not connected" << std::endl; }
		virtual void handle_request_visualize_stage( const std::string& label_map_name, const std::string& mesh_name, unsigned short stage_number, const std::string& server_ip_address, bool is_public_vis_session ) { std::cout << "handle_request_visualize_stage not connected" << std::endl; }
		virtual void handle_update_model_config_info( const std::string& label_map_name, 
			const std::string& mesh_name,
			const std::string& description, const unsigned short stage_number, 
			const std::string& mats, const std::string& mat_names, const std::string& mat_radii, 
			const std::string& refinement_levels, const std::string& constant_sizing_value, const std::string& max_sizing_field, 
			const std::string& num_particle_iters, const std::string& tetgen_joined_vol_flags ) { std::cout << "handle_update_model_config_info not connected" << std::endl; }
		virtual void handle_upload_nrrd_create_project( const std::string& full_project_name, const std::string& file_name,
			size_t file_size, unsigned char* file_data, const std::string& mats, const std::string& mat_names, 
			const std::string& mat_radii, const std::string& refinement_levels, const std::string& constant_sizing_value, const std::string& max_sizing_field, 
			const std::string& num_particle_iters, const std::string& tetgen_joined_vol_flags ) { std::cout << "handle_upload_nrrd_create_project not connected" << std::endl; }

		virtual void handle_transfer_file_fragment( const std::string& label_map_name, const std::string& mesh_name, const std::string& file_name, size_t file_size, unsigned char* file_data, unsigned short fragment_size ) { std::cout << "handle_transfer_file_fragment not connected" << std::endl; }

		virtual void handle_request_file_fragment( const std::string& label_map_name, const std::string& mesh_name, const std::string& file_path, size_t bytes_transferred ) { std::cout << "handle_request_file not connected" << std::endl; }
		virtual void handle_clear_error_message( const std::string& full_project_name ) { std::cout << "handle_clear_error_message not connected" << std::endl; }

		virtual void handle_query_primary_project_message( const std::string& full_project_name) { std::cout << "handle_query_primary_project_message not connected" << std::endl; }
		virtual void handle_request_delete_label_map_message( const std::string& label_map_name ) { std::cout << "handle_request_delete_label_map_message not connected" << std::endl; }
		virtual void handle_request_delete_mesh_message( const std::string& label_map_name,
			const std::string& mesh_name ) { std::cout << "handle_request_delete_mesh_message not connected" << std::endl; }
		virtual void handle_upload_nrrd_create_primary_project( 
			const std::string& owner, const std::string& group,
			const std::string& full_project_name, const std::string& file_name, 
			size_t file_size, unsigned char* file_data){ std::cout << "handle_upload_nrrd_create_primary_project not connected" << std::endl; }

		virtual void handle_request_join_vis_session_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address,
			unsigned short stage_number, unsigned short port,const std::string& owner_name,
			bool is_public_vis_session, unsigned short group_id,
			unsigned int sub_project_id ){ std::cout << "handle_request_join_vis_session_message not connected" << std::endl; }

		virtual void handle_request_resume_vis_session_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address,
			unsigned short stage_number, unsigned short port,const std::string& owner_name,
			const std::string& member_name,
			bool is_public_vis_session, unsigned short group_id,
			unsigned int sub_project_id ){ std::cout << "handle_request_resume_vis_session_message not connected" << std::endl; }

		virtual void handle_request_leave_vis_session_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address,
			unsigned short stage_number, unsigned short port,const std::string& owner_name,
			bool is_public_vis_session, unsigned short group_id,
			unsigned int sub_project_id ){ std::cout << "handle_request_leave_vis_session_message not connected" << std::endl; }

		virtual void handle_cancel_file_transfer_message() { std::cout << "handle_cancel_file_transfer_message not connected" << std::endl; }
		virtual void handle_modify_user_groups( const std::string& username, const std::vector< std::string >& groups ) { std::cout << "handle_modify_user_groups not connected" << std::endl; }
		virtual void handle_add_user_group( const std::string& group_name ) { std::cout << "handle_add_user_group not connected" << std::endl; }

	private:
		void register_handlers();
		void unregister_handlers();

		void connect_keyboard_handler( MessageHandle msg );
		void connect_mouse_handler( MessageHandle msg );
		void connect_clipping_handler( MessageHandle msg );
		void connect_resize_handler( MessageHandle msg );
		void connect_quality_adjust_handler( MessageHandle msg );
		void connect_exit_scirun_handler( MessageHandle msg );
		void connect_ping_handler( MessageHandle msg );
		void connect_quality_of_service_handler( MessageHandle msg );
		void connect_latency_handler( MessageHandle msg );

		void connect_authenticate_client_handler( MessageHandle msg );
		void connect_add_user_handler( MessageHandle msg );
		void connect_add_label_map_message_handler( MessageHandle msg );
		void connect_add_mesh_handler( MessageHandle msg );
		void connect_request_model_config_info_handler( MessageHandle msg );
		void connect_request_nrrd_list_handler( MessageHandle msg );
		void connect_request_run_stage_handler( MessageHandle msg );
		void connect_request_terminate_stage_handler( MessageHandle msg );
		void connect_request_visualize_stage_handler( MessageHandle msg );
		void connect_update_model_config_info_handler( MessageHandle msg );
		void connect_upload_nrrd_create_project_handler( MessageHandle msg );
		void connect_transfer_file_fragment_handler( MessageHandle msg );
		void connect_request_download_file_handler( MessageHandle msg );
		void connect_clear_error_message_handler( MessageHandle msg );
		void connect_query_primary_project_message_handler( MessageHandle msg );
		void connect_request_delete_mesh_handler( MessageHandle msg );
		void connect_request_delete_label_map_handler( MessageHandle msg );
		void connect_upload_nrrd_create_primary_project_message_handler( MessageHandle msg );
		void connect_request_join_vis_session_message_handler( MessageHandle msg );
		void connect_request_resume_vis_session_message_handler( MessageHandle msg );
		void connect_request_leave_vis_session_message_handler( MessageHandle msg );
		void connect_handle_cancel_file_transfer_handler( MessageHandle msg );
		void connect_modify_user_groups_handler( MessageHandle msg );
		void connect_add_user_group_handler( MessageHandle msg );
	};
}

#endif