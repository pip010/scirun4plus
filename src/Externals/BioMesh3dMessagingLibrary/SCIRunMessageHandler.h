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


#ifndef SCIRUN_MESSAGE_H
#define SCIRUN_MESSAGE_H

#include "MessageHandler.h"

#include <Core/Codec/CompressedImage3.h>

namespace BioMesh3d
{

class SCIRunMessageHandler : public MessageHandler
{
public:
	SCIRunMessageHandler();
	~SCIRunMessageHandler();

protected:

	virtual void handle_compressed_image( Core::CompressedImage3Handle compressed_image ) {}
	virtual void handle_quality_of_service( unsigned short messages_received ){}
	virtual void handle_latency( unsigned long timestamp ){}
	virtual void handle_notify_server_of_port( unsigned short port ){}

private:
	void register_handlers();
	void unregister_handlers();

	void connect_compressed_image_handler( MessageHandle msg );
	void connect_quality_of_service_handler( MessageHandle msg );
	void connect_latency_handler( MessageHandle msg );

	void connect_notify_server_of_port( MessageHandle msg );
};

}

#endif