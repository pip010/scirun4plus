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

#include <boost/shared_ptr.hpp>

#include <Core/Codec/Encoder.h>
#include <Core/Utils/AtomicCounter.h> 

namespace Core
{
	
class EncoderPrivate
{
public:
	EncoderPrivate() : 
		codec_serial_id_( Counter++ ),
		bitrate_( 0 ), 
		keyframe_interval_( 1 ), 
		parameters_changed_( false )
	{ }
	
	CodecSerialID codec_serial_id_;
	unsigned int bitrate_;
	unsigned int keyframe_interval_;
	bool parameters_changed_;

	static AtomicCounter Counter;
};

AtomicCounter EncoderPrivate::Counter;

Encoder::Encoder( CodecTypeID codec_type ) :
	Codec(codec_type),
	private_( new EncoderPrivate )
{
}

Encoder::~Encoder()
{	
}
	
CodecSerialID Encoder::get_codec_serial_id() const
{
	lock_type lock( this->get_mutex() );
	return this->private_->codec_serial_id_;
}

// Set the bitrate per second of the encoder.  This depends on the target framerate.
void Encoder::set_bitrate( unsigned int bitrate )
{
	lock_type lock( this->get_mutex() );

	if ( bitrate != this->private_->bitrate_) 
	{
		this->private_->bitrate_ = bitrate;
		this->set_parameters_changed();
	}
}

unsigned int Encoder::get_bitrate() const
{
	lock_type lock( this->get_mutex() );
	return this->private_->bitrate_;
}

// Set the interval between keyframes.
void Encoder::set_keyframe_interval( unsigned int keyframe_interval )
{
	lock_type lock( this->get_mutex() );
	if ( keyframe_interval != this->private_->keyframe_interval_ ) 
	{
		this->private_->keyframe_interval_ = keyframe_interval;
		this->set_parameters_changed();
	}
}

unsigned int Encoder::get_keyframe_interval() const
{
	lock_type lock( this->get_mutex() );
	return this->private_->keyframe_interval_;
}

void Encoder::set_parameters_changed()
{
	lock_type lock( this->get_mutex() );
	this->private_->parameters_changed_ = true;
}

void Encoder::reset_parameters_changed()
{
	lock_type lock( this->get_mutex() );
	this->private_->parameters_changed_ = false;
}

bool Encoder::get_parameters_changed() const
{
	lock_type lock( this->get_mutex() );
	return this->private_->parameters_changed_;
}

void Encoder::update_codec_serial_id()
{
	lock_type lock( this->get_mutex() );
	this->private_->codec_serial_id_ = this->private_->Counter++;
}


}// end namespace Core
