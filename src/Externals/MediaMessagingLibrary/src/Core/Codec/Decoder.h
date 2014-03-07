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

#ifndef CORE_CODEC_DECODER_H
#define CORE_CODEC_DECODER_H

#include <boost/shared_ptr.hpp>

#include <Core/Codec/Codec.h>
#include <Core/Codec/CompressedImage3.h>

#include <string>

namespace Core
{
	
class Decoder;
typedef boost::shared_ptr<Decoder> DecoderHandle; 

class DecoderPrivate;
typedef boost::shared_ptr<DecoderPrivate> DecoderPrivateHandle;

// Is not threadsafe	
class Decoder : public Codec
{
public:
	
	Decoder( CodecTypeID );
	virtual ~Decoder();

	CodecSerialID get_codec_serial_id() const;

	// Encodes the raw image represented by input into the compressed image
	// represented by output.
	virtual bool decode_image( CompressedImage3Handle input, Image3Handle& output ) = 0;
	
protected:
	DecoderPrivateHandle private_;
	void set_codec_serial_id( CodecSerialID );
	
};

}// end namespace Core

#endif // CORE_CODEC_DECODER_H
