#ifndef SCIRUN_MESSAGE_COMPOSER_H
#define SCIRUN_MESSAGE_COMPOSER_H

#include "Message.h"
#include <Core/Codec/CompressedImage3.h>

namespace BioMesh3d
{

class SCIRunMessageComposer
{
protected:
	MessageHandle compose_image_message( Core::CompressedImage3Handle image );
	MessageHandle compose_quality_of_service_message( unsigned short messages_received );
	MessageHandle compose_latency_message( unsigned long timestamp );
	MessageHandle compose_notify_server_of_port_message( unsigned short port );
};

} // end namespace BioMesh3d

#endif
