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
#ifndef CORE_CODEC_MPEGENCODER_H
#define CORE_CODEC_MPEGENCODER_H

#include <Core/Codec/Encoder.h>

#include <boost/shared_ptr.hpp>

namespace Core
{

class MPEGEncoder;
typedef boost::shared_ptr<MPEGEncoder> MPEGEncoderHandle;

class MPEGEncoderPrivate;
typedef boost::shared_ptr<MPEGEncoderPrivate> MPEGEncoderPrivateHandle;

const CodecTypeID MPEG_CODEC_TYPE_C="MPEG Codec";
	
class MPEGEncoder : public Encoder
{
	friend class MPEGEncoderPrivate;

public:
	MPEGEncoder();
	virtual ~MPEGEncoder();

	virtual bool encode_image( Image3Handle input, CompressedImage3Handle& output );
	
	// -- internals --
private:
	MPEGEncoderPrivateHandle private_;

};

} // end namespace Core

#endif // CORE_CODEC_MPEGENCODER_H
