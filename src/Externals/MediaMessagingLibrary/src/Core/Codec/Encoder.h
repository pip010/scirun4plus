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

#ifndef CORE_CODEC_ENCODER_H
#define CORE_CODEC_ENCODER_H

#include <boost/shared_ptr.hpp>

#include <Core/Codec/Codec.h>
#include <Core/Codec/CompressedImage3.h>

#include <string>

namespace Core
{
	
class Encoder;
class EncoderPrivate;

typedef boost::shared_ptr<Encoder> EncoderHandle; 
typedef boost::shared_ptr<EncoderPrivate> EncoderPrivateHandle; 

// Is not threadsafe	
class Encoder : public Codec
{
public:
	
	Encoder( CodecTypeID );
	virtual ~Encoder();

	CodecSerialID get_codec_serial_id() const;
	
	// Encodes the raw image represented by input into the compressed image
	// represented by output.
	virtual bool encode_image( Image3Handle input, CompressedImage3Handle& output ) = 0;
	
	// Set the bitrate per second of the encoder.  This depends on the target framerate.
	void set_bitrate( unsigned int bitrate );
	
	unsigned int get_bitrate() const;

	// Set the interval between keyframes.
	void set_keyframe_interval( unsigned int keyframe_interval );
	
	unsigned int get_keyframe_interval() const;
	
	// -- internals --
protected:
	void set_parameters_changed();
	void reset_parameters_changed();
	bool get_parameters_changed() const;

	// called by the derived class to indicate that
	// a new encoder has been configured.
	void update_codec_serial_id();
	
private:
	
	EncoderPrivateHandle private_;
};

}// end namespace Core

#endif // CORE_CODEC_ENCODER_H
