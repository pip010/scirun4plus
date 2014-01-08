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

// Core includes
#include <Core/Utils/Log.h>
#include <Core/Utils/StdVectorBuffer.h>

namespace Core
{

StdVectorBuffer::StdVectorBuffer()
{
}

StdVectorBuffer::~StdVectorBuffer()
{
}

void* StdVectorBuffer::get_data_ptr()
{
	return reinterpret_cast<void *>( &( this->buffer_[ 0 ] ) );
}

size_t StdVectorBuffer::get_data_size() const
{
	return this->buffer_.size();
}

bool StdVectorBuffer::resize( size_t size )
{
	try {
		this->buffer_.resize( size );
	}
	catch (...)
	{
		CORE_LOG_ERROR( "Failed to resize StdVectorBuffer." );
		return false;
	}
	
	return true;
}

BufferHandle StdVectorBuffer::New( size_t size )
{

	if ( size == 0) size = 1;
	
	StdVectorBuffer* buffer = 0;
	
	try 
	{		
		buffer = new StdVectorBuffer;
	}
	catch ( ... ) 
	{
		CORE_LOG_ERROR( "Could not allocate StdVectorBuffer." );
		return BufferHandle();
	}
	
	// if we cannot resize the vector, log an error and return
	// an empty handle.
	if (! buffer->resize( size ) )
	{
		return BufferHandle();
	}
	return BufferHandle( buffer );
}
	
} // namespace Core

