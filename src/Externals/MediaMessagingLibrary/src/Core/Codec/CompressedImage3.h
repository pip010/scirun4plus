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

#ifndef CORE_CODEC_COMPRESSEDIMAGE3_H
#define CORE_CODEC_COMPRESSEDIMAGE3_H

#include <boost/shared_ptr.hpp>

#include <Core/Codec/Codec.h>
#include <Core/Codec/Image3.h>

namespace Core {
	
class CompressedImage3;
typedef boost::shared_ptr<CompressedImage3> CompressedImage3Handle;

class CompressedImage3: public Image3
{

protected:
	// users should not invoke the constructor directly -- they should use New below
	CompressedImage3( size_t width, size_t height, BufferHandle buf );

	CompressedImage3( size_t width, size_t height, CodecTypeID codec_type,
		CodecSerialID codec_serial_id, BufferHandle buffer );
public:
	virtual ~CompressedImage3();

	// accessors
	bool get_is_keyframe() const;
	void set_is_keyframe( bool is_keyframe );
	CodecTypeID get_codec_type() const;
	void set_codec_type( CodecTypeID codec_type );
	CodecSerialID get_codec_serial_id () const;
	void set_codec_serial_id( CodecSerialID codec_serial_id );
	
private:

	bool is_keyframe_;
	CodecTypeID codec_type_;
	CodecSerialID codec_serial_id_;
	
public:
	
	static CompressedImage3Handle New( size_t width, size_t height, BufferHandle buf );

	static CompressedImage3Handle New( size_t width, size_t height, CodecTypeID codec_type,
		CodecSerialID codec_serial_id, BufferHandle buf );
}; // class CompressedImage3

} // namespace Core

#endif // CORE_CODEC_COMPRESSEDIMAGE3_H
