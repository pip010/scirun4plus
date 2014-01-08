/*
 For more information, please see: http://software.sci.utah.edu

 The MIT License

 Copyright (c) 2011 Scientific Computing and Imaging Institute,
 University of Utah.


 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the "Software"),
 to deal in the Software without restriction, including without limitation
 the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
 */


#ifndef MESSAGE_H
#define MESSAGE_H

#include <boost/algorithm/string.hpp>
#include <boost/asio.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/noncopyable.hpp>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>

#include <iostream>
#include <map>

#include <Core/Utils/Buffer.h>

#define BINARY_DATA_SIZE 1000000
#define MESSAGE_HEADER_SIZE 6

namespace BioMesh3d
{

class Message;
class TCPCommunicator;
class RTPMediaCommunication;
typedef boost::shared_ptr< Message > MessageHandle;

namespace StageStates
{
	enum StageState
	{
		UNINITIALIZED_E,
		COMPLETED_E,
		RUNNING_E,
		ERROR_E
	};
}// end namespace StageStates

namespace HeaderOffsets
{
	enum Offset
	{
		SIZE_E = 0,
		TYPE_E = 4,
	};
}// end namespace HeaderOffsets

namespace MessageTypes
{
	enum MessageType
	{
		INVALID_E = 0,
		PING_E = 1,
		// Messages created by the Server
		SERVER_SET_HOSTNAME_E = 100,
		SERVER_STAGE_COMPLETE_E = 101,
		SERVER_STAGE_STARTED_E = 102,
		SERVER_STAGE_STOPPED_E = 103,
		SERVER_STAGE_ERROR_E = 104,
		SERVER_STAGE_TERMINATED_E = 105,
		SERVER_VISUALIZATION_TERMINATED_E = 106,
		SERVER_VISUALIZATION_STARTED_E = 107,
		// Messages created by the ServerManager
		SERVICES_MANAGER_STAGE_COMPLETE_E = 200,
		SERVICES_MANAGER_STAGE_ERROR_E = 201,
		SERVICES_MANAGER_STAGE_TERMINATED_E = 202,
		SERVICES_MANAGER_REQUEST_RUN_STAGE_E = 203,
		SERVICES_MANAGER_REQUEST_STOP_STAGE_E = 204,
		SERVICES_MANAGER_REQUEST_TERMINATE_STAGE_E = 205,
		SERVICES_MANAGER_REQUEST_VISUALIZE_STAGE_E = 206,
		SERVICES_MANAGER_INITIATE_FILE_UPLOAD_E = 207,
		SERVICES_MANAGER_IS_CLIENT_AUTHENTICATED_E = 208,
		SERVICES_MANAGER_MESH_STATUS_E = 209,
		SERVICES_MANAGER_SERVER_LIST_E = 210,
		SERVICES_MANAGER_START_FILE_DOWNLOAD_E = 211,
		SERVICES_MANAGER_VISUALIZATION_STARTED_E = 212,
		SERVICES_MANAGER_FILE_TRANSFER_COMPLETE_E = 213,
		SERVICES_MANAGER_MODEL_CONFIG_INFO_E = 214,
		SERVICES_MANAGER_TRANSFER_FILE_E = 215,
		SERVICES_MANAGER_ANSWER_QUERY_PRIMARY_PROJECT_E = 216,
		SERVICES_MANAGER_ANSWER_CREATE_PRIMARY_PROJECT_E = 217,
		SERVICES_MANAGER_PROJECT_PATH_E = 218,
		SERVICES_MANAGER_DELETE_LABEL_MAP_E = 219,
		SERVICES_MANAGER_DELETE_MESH_E = 220,
		SERVICES_MANAGER_ADD_MESH_E = 221,
		SERVICES_MANAGER_ANSWER_UPDATE_SUBPROJECT_E = 222,
		SERVICES_MANAGER_ANSWER_QUERY_GROUPS_E = 223,
		SERVICES_MANAGER_CLIENT_LOGIN_FAIL_E = 224,
		SERVICES_MANAGER_VISUALIZATION_SESSION_E = 225,
		SERVICES_MANAGER_LEAVE_VIS_SESSION_E = 226,
		SERVICES_MANAGER_JOIN_VIS_SESSION_E = 227,
		SERVICES_MANAGER_VISUALIZATION_RESUMED_E = 228,
		SERVICES_MANAGER_ALERT_MESSAGE_E = 229,
		SERVICES_MANAGER_FILE_TRANSFER_CANCELED_E = 230,
		SERVICES_MANAGER_RSA_PEM_PUBLIC_KEY_E = 231,
		SERVICES_MANAGER_USERS_INFO_E = 232,
		BROKER_DELETE_LABEL_MAP_E = 233,

		// Messages created by SCIRun
		SCIRUN_IMAGE_E = 300,
		// TODO: Remove the next four
		SCIRUN_MPEG4_E = 301,
		SCIRUN_VP8_E = 302,
		SCIRUN_JPEG_E = 303,
		SCIRUN_PNG_E = 304,
		SCIRUN_QUALITY_OF_SERVICE_E = 305,
		SCIRUN_NOTIFY_SERVER_OF_PORT_E = 306,
		SCIRUN_LATENCY_E = 307,

		// Messages created by the Client
		CLIENT_AUTHENTICATE_CLIENT_E = 400,
		CLIENT_ADD_USER_E = 401,
		CLIENT_ADD_LABEL_MAP_E = 402,
		CLIENT_ADD_MESH_E = 403,
		CLIENT_REQUEST_MODEL_CONFIG_INFO_E = 404,
		CLIENT_REQUEST_NRRD_LIST_E = 405,
		CLIENT_REQUEST_RUN_STAGE_E = 406,
		CLIENT_REQUEST_TERMINATE_STAGE_E = 407,
		CLIENT_REQUEST_VISUALIZE_STAGE_E = 408,
		CLIENT_TRANSFER_FILE_E = 409,
		CLIENT_REQUEST_UPLOAD_FILE_E = 410,
		CLIENT_REQUEST_DOWNLOAD_FILE_E = 411,

		CLIENT_REQUEST_QUERY_PRIMARY_PROJECT_E = 500,
		CLIENT_REQUEST_CREATE_PRIMARY_PROJECT_UPLOAD_FILE_E = 501,
		CLIENT_REQUEST_DELETE_LABEL_MAP_E = 502,
		CLIENT_REQUEST_DELETE_MESH_E = 503,
		CLIENT_REQUEST_LEAVE_SESSION_E = 504,
		CLIENT_REQUEST_JOIN_SESSION_E = 505,
		CLIENT_REQUEST_RESUME_SESSION_E = 506,
		CLIENT_CANCEL_FILE_TRANSFER_E = 507,

		CLIENT_CLEAR_ERROR_MESSAGE_E = 600,
		CLIENT_UPDATE_MODEL_CONFIG_INFO_E = 601,
		CLIENT_KEYBOARD_E = 602,
		CLIENT_MOUSE_E = 603,
		CLIENT_CLIPPING_E = 604,
		CLIENT_RESIZE_E = 605,
		CLIENT_QUALITY_ADJUST_E = 606,
		CLIENT_EXIT_SCIRUN_E = 607,
		CLIENT_QUALITY_OF_SERVICE_E = 608,
		CLIENT_LATENCY_E = 609,
		CLIENT_MODIFY_USER_GROUPS_E = 610,
		CLIENT_ADD_USER_GROUP_E = 611,

		NOOP_E = 700,

		// DONT PUT ANYTHING PAST MESSAGE_TYPE_SIZE_E
		MESSAGE_TYPE_SIZE_E = 1000 // the max value in this enum
	};
} // end namespace MessageTypes

namespace ParameterTypes
{
	enum ParamType
	{
		INVALID_E = 0,
		FULL_PROJECT_NAME_E = 1,
		LABEL_MAP_NAME_E = 2,
		MESH_NAME_E = 3,
		FULL_PROJECT_PATH_E = 4,
		STAGE_NUMBER_E = 5,
		STAGE_STATE_E = 6,
		STAGES_COMPLETE_E = 7,
		STAGE_RUNNING_E = 8,
		STAGE_ERROR_E = 9,
		DATE_PROCESSED_E = 10,
		ERROR_MESSAGE_E = 11,
		NRRD_FILE_NAMES_E = 12,
		HOSTNAME_E = 13,
		USERNAME_E = 14,
		PASSWORD_E = 15,
		FIRST_NAME_E = 16,
		LAST_NAME_E = 17,
		VALIDATION_MESSAGE_E = 18,
		PORT_E = 19,
		CLIENT_IP_ADDRESS_E = 20,
		FILE_NAME_E = 21,
		FILE_SIZE_E = 22,
		FRAGMENT_NUMBER_E = 23,
		FRAGMENT_SIZE_E = 24,
		BYTES_TRANSFERRED_E = 25,
		BINARY_DATA_E = 26,
		KEY_PRESSED_E = 27,
		CTRL_PRESSED_E = 28,
		ALT_PRESSED_E = 29,
		SHIFT_PRESSED_E = 30,
		MOUSE_LEFT_E = 31,
		MOUSE_MIDDLE_E = 32,
		MOUSE_RIGHT_E = 33,
		MOUSE_X_E = 34,
		MOUSE_Y_E = 35,
		WIDTH_E = 36,
		HEIGHT_E = 37,
		CLIPPING_FLOAT_E = 38,
		CLIPPING_DIMENSION_E = 39,
		CLIPPING_PLANE_E = 40,
		GOP_E = 41,
		QUALITY_E = 42,
		SERVER_IP_ADDRESS_E = 43,
		MODEL_MATS_E = 44,
		MODEL_MAT_NAMES_E = 45,
		MODEL_MAT_RADII_E = 46,
		MODEL_REFINEMENT_LEVELS_E = 47,
		MODEL_MAX_SIZING_FIELD_E = 48,
		MODEL_NUM_PARTICLE_ITERS_E = 49,
		MODEL_TETGEN_JOINED_VOL_FLAGS_E = 50,
		AUTHENTICATED_E = 51,
		PROJECT_PATH_E = 52,
		SUB_PROJECT_ID_E = 53,

		PRIMARY_PROJECT_FIND_E = 54,
		USER_GROUPS_E = 55,
		GROUP_NAME_E = 56,
		VIS_MEMBERS_E = 57,
		DESCRIPTION_E = 58,
		OPERATION_ANSWER_E = 59,
		OPERATION_DESCRIPTION_E = 60,
		IS_PUBLIC_VIS_SESSION_E = 61,
		LOGIN_FAIL_DESCRIPTION_E = 62,
		MEMBERNAME_E = 63,
		MESSAGES_RECEIVED_E = 64,
		RSA_PEM_PUBLIC_KEY_E = 65,
		USERS_INFO_E = 66,
		TIMESTAMP_E = 67,

		START_STAGE_NUMBER_E = 68,
		END_STAGE_NUMBER_E = 69,
		
		ALERT_E = 70,
		SERVER_IP_ADDRESS_LIST_E = 71,
		SERVER_HOSTNAME_LIST_E = 72,

		CODEC_SERIAL_E = 73,
		CODEC_TYPE_E = 74,

		CLIPPING_ENABLE_PLANE_E = 75,
		CLIPPING_SHOW_FRAME_E = 76,
		CLIPPING_REVERSE_NORMAL_E = 77,
		CLIPPING_X_E = 78,
		CLIPPING_Y_E = 79,
		CLIPPING_Z_E = 80,
		CLIPPING_D_E = 81,
		
		// I realize this is out of order, but I didn't 
		// want to renumber everything.
		MODEL_CONSTANT_SIZING_VALUE_E = 82,
		
		// DONT PUT ANYTHING PAST PARAM_TYPES_SIZE_E
		PARAM_TYPE_SIZE_E = 500
	};
} // end namespace MessageParameterTypes



class Message
{
public:

	Message();
	~Message();

	// TODO: Hmmm.....
	void reset_message();

	template < class T >
	bool add_message_parameter( const int key, const T& value )
	{
		boost::mutex::scoped_lock lock( m_message_values_mutex_ );
		bool is_valid = true;
		try
		{
			std::string value_string = boost::lexical_cast< std::string >( value );
			boost::algorithm::replace_all( value_string, "&", "%26" );
			boost::algorithm::replace_all( value_string, "=", "%3D" );
			this->m_message_values_[key] = value_string;
		}
		catch ( boost::bad_lexical_cast& exc )
		{
			std::cout << "Could not add a parameter to the message: " << exc.what() << std::endl;
			is_valid = false;
		}
		catch ( std::exception& exc )
		{
			std::cout << "Standard error: " << exc.what() << std::endl;
			is_valid = false;
		}
		is_valid = ( get_message_type() != MessageTypes::INVALID_E ) ? is_valid : false; // if we aren't valid then return false
		return is_valid;
	}


	// if the key exists then set value with the correct value and
	// return true. if the key exists then set the value to an empty
	// string and return false.
	template < class T >
	bool retrieve_message_parameter( const int key, T& value )
	{
		boost::mutex::scoped_lock lock( m_message_values_mutex_ );
		bool is_valid = false;
		try
		{
			if ( this->m_message_values_.find( key ) != this->m_message_values_.end() )
			{
				std::string value_string = this->m_message_values_[ key ];
				boost::algorithm::replace_all( value_string, "%26", "&" );
				boost::algorithm::replace_all( value_string, "%3D", "=" );
				try
				{
					value = boost::lexical_cast< T >( value_string );
				}
				catch ( std::exception& exc )
				{
					std::cout << "Lexical cast error: " << exc.what() << std::endl;
				}
				is_valid = true;
			}
		}
		catch ( std::exception& exc )
		{
			std::cout << "Could not retrieve parameter due to bad lexical cast: " << exc.what() << std::endl;
		}
		is_valid = ( get_message_type() != MessageTypes::INVALID_E ) ? is_valid : false; // if we aren't valid then return false
		return is_valid;	
	}

	// TODO: Replace with get_buffer
	boost::uint8_t* get_binary_buf();
	boost::uint32_t get_binary_size();
	
	Core::BufferHandle get_binary_buffer();
	void set_binary_buffer( Core::BufferHandle buffer );
	
	// TODO: Remove this one
	void set_binary_size( boost::uint32_t bin_size );
	
	// TODO: Should not really do this
	void copy_binary_data( boost::uint8_t* data, size_t data_size );

	MessageTypes::MessageType get_message_type();
	void set_message_type( MessageTypes::MessageType message_type );
	std::string serialize_parameters();
	
	
private:

	MessageTypes::MessageType m_message_type_;
	
	// TODO: Test without this mutex
	boost::mutex m_message_values_mutex_;
	std::map< int, std::string > m_message_values_;
	
	// Replace with BufferHandle
	Core::BufferHandle buffer_;
};

	
} // end namespace Message

#endif
