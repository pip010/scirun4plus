#ifndef CLIENT_MESSAGE_COMPOSER_H
#define CLIENT_MESSAGE_COMPOSER_H

#include "Message.h"

namespace BioMesh3d
{

	class ClientMessageComposer
	{
	public:
		ClientMessageComposer();
		~ClientMessageComposer();

	protected:

		MessageHandle compose_keyboard_message( const char key_pressed, 
			bool ctrl, bool alt, bool shift );
		MessageHandle compose_mouse_message( short _x,
			short _y, bool left_click, bool middle_click, bool right_click,
			bool ctrl, bool alt, bool shift, unsigned long timestamp );
		MessageHandle compose_clipping_message( unsigned short plane, bool enable, 
			bool show_frame, bool reverse_normal, float x, float y, float z, float d );
		MessageHandle compose_resize_message( const unsigned short width, 
			const unsigned short height );
		MessageHandle compose_quality_adjust_message( const unsigned int quality,
			const unsigned short gop );
		MessageHandle compose_exit_scirun_message();
		MessageHandle compose_ping_message();
		MessageHandle compose_quality_of_service_message( unsigned short messages_received );
		MessageHandle compose_latency_message( unsigned long timestamp );


		MessageHandle compose_authenticate_client_message( const std::vector< char >& username, const std::vector< char >& password );
		MessageHandle compose_add_user_message( const std::vector< char >& username, const std::vector< char >& password, const std::string& first_name, const std::string& last_name );
		MessageHandle compose_add_label_map_message( 
			const std::string& primary_project_name, 
			const std::string& group_name,
			const std::string& description );

		MessageHandle compose_add_mesh_message( const std::string& label_map_name, 
			const std::string& mesh_name );

		MessageHandle compose_request_model_config_info_message( const std::string& label_map__name, const std::string& mesh_name );
		MessageHandle compose_request_nrrd_list_message();
		MessageHandle compose_request_run_stage_message( const std::string& label_map_name, 
			const std::string& mesh_name, const std::string& server_ip_address, 
			const unsigned short stage_number );
		MessageHandle compose_request_terminate_stage_message( const std::string& label_map_name, 
			const std::string& mesh_name, unsigned short stage_number );
		MessageHandle compose_request_visualize_stage_message( const std::string& label_map_name, 
			const std::string& mesh_name, unsigned short stage_number, 
			const std::string server_ip_address, bool is_public_vis_session );
		MessageHandle compose_update_model_config_info_message( const std::string label_map_name, 
			const std::string mesh_name, const std::string& description, 
			const unsigned short stage_number, const std::string mats, const std::string mat_names,
			const std::string mat_radii, const std::string refinement_levels, 
			const std::string constant_sizing_value, 
			const std::string max_sizing_field, const std::string num_particle_iters, 
			const std::string tetgen_joined_vol_flags );
		MessageHandle compose_upload_nrrd_create_project_message( const std::string& full_project_path, const std::string& file_name, 
			size_t file_size, unsigned char* file_data, const std::string& mats, const std::string& mat_names, 
			const std::string& mat_radii, const std::string& refinement_levels, 
			const std::string& constant_sizing_value, const std::string& max_sizing_field, 
			const std::string& num_particle_iters, const std::string& tetgen_joined_vol_flags );
		MessageHandle compose_transfer_file_fragment_message( const std::string& label_map_name, const std::string& mesh_name, const std::string& file_name, size_t file_size, unsigned char* file_data, unsigned short fragment_size );
		MessageHandle compose_request_file_message( const std::string& label_map_name, const std::string& mesh_name, const std::string& file_path, size_t bytes_transferred );


		MessageHandle compose_query_primary_project_message( 
			const std::string& primary_project_name);

		MessageHandle compose_upload_nrrd_create_primary_project_message( 
			const std::string& owner, const std::string& group,
			const std::string& full_project_path, const std::string& file_name,
			size_t file_size,unsigned char* file_data );

		MessageHandle compose_request_delete_label_map_message( 
			const std::string& label_map_name );
		MessageHandle compose_request_delete_mesh_message( 
			const std::string& label_map_name, const std::string& mesh_name );
		MessageHandle compose_clear_project_error( const std::string& full_project_path );
		MessageHandle compose_leave_vis_session_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address, 
			const unsigned short stage_number,
			const unsigned short port,
			const std::string& member_name,
			bool is_public_vis_session,
			unsigned short group_id,
			unsigned int sub_project_id);

		MessageHandle compose_join_vis_session_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address, 
			const unsigned short stage_number,
			const unsigned short port,
			const std::string& member_name,
			bool is_public_vis_session,
			unsigned short group_id,
			unsigned int sub_project_id);

		MessageHandle compose_resume_vis_session_message( 
			const std::string& full_project_name, 
			const std::string& server_ip_address, 
			const unsigned short stage_number,
			const unsigned short port,
			const std::string& owner_name,
			const std::string& member_name,
			bool is_public_vis_session,
			unsigned short group_id,
			unsigned int sub_project_id);

		MessageHandle compose_cancel_file_transfer_message();
		MessageHandle compose_modify_user_groups( const std::string& username, const std::vector< std::string >& groups );
		MessageHandle compose_add_user_group( const std::string& group_name );
	private:
		
	};

}// end namespace Message
#endif