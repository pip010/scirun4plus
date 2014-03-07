#include "SCIRunMessageComposer.h"
#include "Message.h"

namespace BioMesh3d
{

MessageHandle SCIRunMessageComposer::compose_image_message( Core::CompressedImage3Handle image )
{
	MessageHandle msg = MessageHandle( new Message() );
	msg->set_message_type( MessageTypes::SCIRUN_IMAGE_E );
	msg->add_message_parameter( ParameterTypes::WIDTH_E, image->get_width() );
	msg->add_message_parameter( ParameterTypes::HEIGHT_E, image->get_height() );
	msg->add_message_parameter( ParameterTypes::CODEC_TYPE_E, image->get_codec_type() );
	msg->add_message_parameter( ParameterTypes::CODEC_SERIAL_E, image->get_codec_serial_id() );
	msg->set_binary_buffer( image->get_buffer() );
	return msg;
}

MessageHandle SCIRunMessageComposer::compose_quality_of_service_message( 
	unsigned short messages_received )
{
	MessageHandle msg = MessageHandle( new Message() );
	msg->set_message_type( MessageTypes::SCIRUN_QUALITY_OF_SERVICE_E );
	msg->add_message_parameter( ParameterTypes::MESSAGES_RECEIVED_E, 
		messages_received );
	return msg;
}

MessageHandle SCIRunMessageComposer::compose_latency_message( unsigned long timestamp )
{
	MessageHandle msg = MessageHandle( new Message() );
	msg->set_message_type( MessageTypes::SCIRUN_LATENCY_E );
	msg->add_message_parameter( ParameterTypes::TIMESTAMP_E, timestamp );
	return msg;
}

MessageHandle SCIRunMessageComposer::compose_notify_server_of_port_message( unsigned short port )
{
	MessageHandle msg = MessageHandle( new Message() );
	msg->set_message_type( MessageTypes::SCIRUN_NOTIFY_SERVER_OF_PORT_E );
	msg->add_message_parameter( ParameterTypes::PORT_E, port );
	return msg;
}

} // end namespace Message
