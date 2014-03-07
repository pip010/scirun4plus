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

#include <boost/bind.hpp>

#include "Message.h"
#include "SCIRunMessageHandler.h"

#include <Core/Codec/CompressedImage3.h>

namespace BioMesh3d
{

	SCIRunMessageHandler::SCIRunMessageHandler()
	{
		register_handlers();
	}

	SCIRunMessageHandler::~SCIRunMessageHandler()
	{
		unregister_handlers();
	}

	void SCIRunMessageHandler::register_handlers()
	{
		register_callback( MessageTypes::SCIRUN_IMAGE_E,
			boost::bind( &SCIRunMessageHandler::connect_compressed_image_handler, this, _1 ) );

		register_callback( MessageTypes::SCIRUN_QUALITY_OF_SERVICE_E,
			boost::bind( &SCIRunMessageHandler::connect_quality_of_service_handler, this, _1 ) );
		register_callback( MessageTypes::SCIRUN_NOTIFY_SERVER_OF_PORT_E,
			boost::bind( &SCIRunMessageHandler::connect_notify_server_of_port, this, _1 ) );
		register_callback( MessageTypes::SCIRUN_LATENCY_E,
			boost::bind( &SCIRunMessageHandler::connect_latency_handler, this, _1 ) );
	}

	void SCIRunMessageHandler::unregister_handlers()
	{
		unregister_callback( MessageTypes::SCIRUN_IMAGE_E );
		unregister_callback( MessageTypes::SCIRUN_QUALITY_OF_SERVICE_E );
		unregister_callback( MessageTypes::SCIRUN_NOTIFY_SERVER_OF_PORT_E );
		unregister_callback( MessageTypes::SCIRUN_LATENCY_E );
	}

	void SCIRunMessageHandler::connect_compressed_image_handler( MessageHandle msg )
	{
		size_t width = 0;
		size_t height = 0;
		Core::CodecTypeID codec_type = Core::UNKNOWN_CODEC_TYPE_C;
		Core::CodecSerialID codec_serial = 0;
		
		if ( !( msg->retrieve_message_parameter( ParameterTypes::WIDTH_E, width ) &&
			  msg->retrieve_message_parameter( ParameterTypes::HEIGHT_E, height ) &&
			  msg->retrieve_message_parameter( ParameterTypes::CODEC_TYPE_E, codec_type ) &&
			  msg->retrieve_message_parameter( ParameterTypes::CODEC_SERIAL_E, codec_serial ) ) )
		{
			std::cerr << "Could not read message" << std::endl;
		}
		
		Core::CompressedImage3Handle cimage = Core::CompressedImage3::New( 
			width, height, codec_type, codec_serial, msg->get_binary_buffer() );

		this->handle_compressed_image( cimage );
	}

	void SCIRunMessageHandler::connect_quality_of_service_handler( MessageHandle msg )
	{
		unsigned short messages_received = 0;
		bool is_valid = msg->retrieve_message_parameter( 
			ParameterTypes::MESSAGES_RECEIVED_E, messages_received );
		if ( is_valid )
		{
			handle_quality_of_service( messages_received );
		}
	}

	void SCIRunMessageHandler::connect_latency_handler( MessageHandle msg )
	{
		unsigned long timestamp = 0;
		bool is_valid = msg->retrieve_message_parameter( ParameterTypes::TIMESTAMP_E, timestamp );
		if ( is_valid )
		{
			handle_latency( timestamp );
		}
	}

	void SCIRunMessageHandler::connect_notify_server_of_port( MessageHandle msg )
	{
		unsigned short port = 0;
		
		std::cerr << "Server received SCIRun port." << std::endl;
		
		bool is_valid = msg->retrieve_message_parameter( ParameterTypes::PORT_E, port );
		if ( is_valid )
		{
			std::cerr << "SCIRun port is: " << port << std::endl;
			handle_notify_server_of_port( port );
		}
	}

}// end namespace Message