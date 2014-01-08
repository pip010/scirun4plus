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



#include <Core/Codec/FFMPEGDecoder.h>
#include <Core/Codec/MPEGEncoder.h>
#include <Core/Codec/VP8Encoder.h>
#include <Core/Utils/Buffer.h>
#include <Core/Utils/Log.h>

// FFMPEG includes
extern "C" {
	
#include <libavutil/pixfmt.h>
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libswscale/swscale.h>

}

#include <boost/shared_ptr.hpp>

#include <algorithm>

namespace Core
{
	
class FFMPEGDecoderPrivate
{

	// Member variables
public:
	// Pointer back to the main class
	FFMPEGDecoder* decoder_;
	
	bool is_initialized_;
	
	AVCodec* codec_;
	AVCodecContext* codec_cxt_;
	AVFrame* frame_yuv_;
	AVFrame* frame_rgb_;
	SwsContext* img_convert_ctx_;
	BufferHandle yuv_buf_;
	CodecID ffmpeg_codec_id_;
	
	// Private functions
public:	
	FFMPEGDecoderPrivate( CodecTypeID type ) :
		decoder_( 0 ),
		is_initialized_( false ),
		codec_( 0 ),
		codec_cxt_( 0 ),
		frame_yuv_( 0 ),
		frame_rgb_( 0 ),
		img_convert_ctx_( 0 )
	{
		// translate our Codec types into the 
		// enumerations defined by FFMPEG
		if ( type == MPEG_CODEC_TYPE_C )
		{
			ffmpeg_codec_id_ = CODEC_ID_MPEG4;
		} else // our default case is Google's VP8
		{ 
			ffmpeg_codec_id_ = CODEC_ID_VP8;
		}
		
		avcodec_init();
		av_register_all();
	}
	
	~FFMPEGDecoderPrivate()
	{
		clear_codec();
	}
	
	
	bool initialize_codec( CodecSerialID serial_id,
						   unsigned int image_width, unsigned int image_height );
	
	void clear_codec();
	
};

bool FFMPEGDecoderPrivate::initialize_codec( CodecSerialID codec_serial_id, unsigned int image_width, unsigned int image_height )
{
	this->clear_codec();
		
	// Make a new context		
	this->codec_cxt_ = avcodec_alloc_context();
	if ( this->codec_cxt_ == 0 )
	{
		// We could not allocate the context
		CORE_LOG_ERROR( "Failed to allocate FFMPEG context." );
		return false;
	}

	// Allocate the RGB image.  Because the buffer will be input as a parameter to
	// the encode_image function, we do not need to pass in a buffer
	this->frame_rgb_ = avcodec_alloc_frame();
	if ( !this->frame_rgb_ )
	{
		CORE_LOG_ERROR( "Could not allocate RGB frame for FFMPEG codec." );
		return false;
	}
	
	// Allocate the YUV image.
	this->yuv_buf_ = Buffer::New( avpicture_get_size( PIX_FMT_YUV420P, 
													  image_width, image_height ) );
	if ( !this->yuv_buf_ )
	{
		CORE_LOG_ERROR( "Could not allocate YUV buffer object for FFMPEG codec." );
		return false;
	}
	
	// Allocate the YUV picture.
	this->frame_yuv_ = avcodec_alloc_frame();
	
	if ( !this->frame_yuv_ )
	{
		CORE_LOG_ERROR( "Could not allocate YUV frame for FFMPEG codec." );
		return false;
	}
	
	avpicture_fill( reinterpret_cast<AVPicture*>( this->frame_yuv_ ), 
				    static_cast<boost::uint8_t*>( this->yuv_buf_->get_data_ptr() ),
				    PIX_FMT_YUV420P, image_width, image_height );
	
	
	this->codec_ = avcodec_find_decoder( this->ffmpeg_codec_id_ );
	if ( this->codec_ == 0 )
	{
		// We could not find the decoder
		CORE_LOG_ERROR( "Could not find FFMPEG codec." );
		return false;
	}
	
	if ( avcodec_open( this->codec_cxt_, this->codec_ ) < 0 ) 
	{
		// Cannot open Codec
		CORE_LOG_ERROR( "Could not open FFMPEG codec." );
		return false;
	}
	
	this->img_convert_ctx_ = sws_getContext( image_width, image_height, 
											 PIX_FMT_YUV420P,
											 image_width, image_height, 
											 PIX_FMT_RGB24, SWS_BICUBIC, 
											 NULL, NULL, NULL );
	
	this->decoder_->set_codec_serial_id( codec_serial_id );
	
	this->is_initialized_ = true;
	return true;
}
	

	
void FFMPEGDecoderPrivate::clear_codec()
{
	if ( !this->is_initialized_ ) return;
	
	this->yuv_buf_.reset();
	av_free( this->frame_rgb_ );
	this->frame_rgb_ = 0;
	av_free( this->frame_yuv_ );
	this->frame_yuv_ = 0;
	avcodec_close( this->codec_cxt_ );
	this->codec_cxt_ = 0;
	av_free( this->codec_cxt_ );
	
	sws_freeContext( this->img_convert_ctx_ ); 
	this->img_convert_ctx_ = 0;
	
	this->is_initialized_ = false;
}

//------ FFMPEGDecoder class

FFMPEGDecoder::FFMPEGDecoder( CodecTypeID type ) :
	Decoder( type ), 
	private_( new FFMPEGDecoderPrivate( type ) )
{
	// Fill in the pointer back to this class
	this->private_->decoder_ = this;
	
}
	
FFMPEGDecoder::~FFMPEGDecoder()
{
}

bool FFMPEGDecoder::decode_image( CompressedImage3Handle input, Image3Handle& output )
{
	lock_type lock( this->get_mutex() );

	// At the outset we free the output object.
	output.reset();

	if ( !this->private_->is_initialized_ || 
		 this->get_codec_serial_id() != input->get_codec_serial_id() ) 
	{
		if (! this->private_->initialize_codec( input->get_codec_serial_id(), input->get_width(), input->get_height() ) )
		{
			CORE_LOG_ERROR( "Failed to encode image." );
			return false;
		}
	}
	
	// unsigned int bytes_remaining = out_img_buf.size();
	unsigned int bytes_remaining = input->get_buffer()->get_data_size();
	int bytes_decoded = 0;
	int frame_finished = 0; // Zero if no frame could be decompressed, otherwise it is non-zero
	boost::uint8_t* inbuf = reinterpret_cast<boost::uint8_t*>(input->get_buffer()->get_data_ptr());
	
	while ( bytes_remaining > 0 )
	{
		/* Decode the next chunk of data */
		bytes_decoded = avcodec_decode_video( this->private_->codec_cxt_, this->private_->frame_yuv_,
											 &frame_finished, inbuf, bytes_remaining);
		
		if ( bytes_decoded < 0 )
		{
			output.reset();
			CORE_LOG_ERROR( "Failed to decode video with avcodec_decode_video" );
			return false;
		}
		
		bytes_remaining -= bytes_decoded;
		inbuf += bytes_decoded;
	}
	
	// Allocate the output RGB image
	BufferHandle output_buffer = Buffer::New( avpicture_get_size( PIX_FMT_RGB24, 
		this->private_->codec_cxt_->width, this->private_->codec_cxt_->height ) );
	
	if (! output_buffer )
	{
		CORE_LOG_ERROR( "Failed to allocate output buffer." );
		return false;
	}	
	
	// Put the input rgb image into the AV Picture data structure.
	avpicture_fill( reinterpret_cast<AVPicture*>( this->private_->frame_rgb_ ), 
				    static_cast<boost::uint8_t*>( output_buffer->get_data_ptr() ),
				    PIX_FMT_RGB24, this->private_->codec_cxt_->width,
					this->private_->codec_cxt_->height );
	
	
	if ( frame_finished == 0 )
	{
		output.reset();
		CORE_LOG_ERROR( "Failed to finish frame in decode_image" );
		return false;
	}
	
	// Convert the YUV image to RGB.
	sws_scale( this->private_->img_convert_ctx_, this->private_->frame_yuv_->data, 
			   this->private_->frame_yuv_->linesize, 0, 
			   this->private_->codec_cxt_->height, this->private_->frame_rgb_->data, 
			   this->private_->frame_rgb_->linesize );

	output = Image3::New( this->private_->codec_cxt_->width, this->private_->codec_cxt_->height, 
		output_buffer );
	if ( !output )
	{
		CORE_LOG_ERROR( "Failed to create output image." );
		return false;
	}

	return true;
}

} // end namespace Core
