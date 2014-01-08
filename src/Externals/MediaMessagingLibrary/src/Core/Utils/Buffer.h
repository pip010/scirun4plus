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

#ifndef CORE_UTILS_BUFFER_H
#define CORE_UTILS_BUFFER_H

// Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

namespace Core
{

class Buffer;
typedef boost::shared_ptr<Buffer> BufferHandle;

class Buffer : public boost::noncopyable
{

public:
	virtual ~Buffer();

	// GET_DATA_PTR:
	// Get the underlying data pointer of this data block
	virtual void* get_data_ptr() = 0;
	
	// GET_DATA_SIZE:
	// Get the size of the underlying data block
	virtual size_t get_data_size() const = 0;
	
	// RESIZE
	// Deallocate previously allocated memory and allocate to new size.
	// return true upon success
	virtual bool resize( size_t size ) = 0;
public:

	// NEW
	// Function for creating a new data block with a certain size
	// NOTE: The new object will be a StdVectorBuffer
	static BufferHandle New( size_t size );
	
};

} // namespace Core

#endif