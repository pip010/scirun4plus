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

#include <Core/Codec/MPEGEncoder.h>
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
	
class MPEGEncoderPrivate
{

	// Member variables
public:
	// Pointer back to the main class
	MPEGEncoder* encoder_;
	
	bool is_initialized_;

	// The size of the image for which we initialize the codec
	unsigned int image_width_;
	unsigned int image_height_;
	
	AVCodec* codec_;
	AVCodecContext* codec_cxt_;
	AVFrame* frame_yuv_;
	AVFrame* frame_rgb_;
	SwsContext* img_convert_ctx_;
	BufferHandle yuv_buf_;
	
	// Private functions
public:	
	MPEGEncoderPrivate() :
		encoder_( 0 ),
		is_initialized_( false ),
		image_width_( 0 ),
		image_height_( 0 ),
		codec_( 0 ),
		codec_cxt_( 0 ),
		frame_yuv_( 0 ),
		frame_rgb_( 0 ),
		img_convert_ctx_( 0 )
	{
		avcodec_init();
		av_register_all();
	}
	
	~MPEGEncoderPrivate()
	{
		clear_codec();
	}
	
	
	bool initialize_codec( unsigned int image_width, unsigned int image_height );
	
	void clear_codec();
	
};

bool MPEGEncoderPrivate::initialize_codec( unsigned int image_width, unsigned int image_height )
{
	this->clear_codec();
	
	// Set the size for the new image series
	this->image_width_ = image_width;
	this->image_height_ = image_height;
	
	// Make a new context		
	this->codec_cxt_ = avcodec_alloc_context();
	if ( this->codec_cxt_ == 0 )
	{
		// We could not allocate the context
		CORE_LOG_ERROR( "Failed to allocate FFMPEG context." );
		return false;
	}
	
	// this is kind of a ghetto way to do this but it works and isn't too much code
	this->codec_cxt_->codec_id = CODEC_ID_MPEG4;
	this->codec_cxt_->codec_type = CODEC_TYPE_VIDEO;
	this->codec_cxt_->bit_rate = this->encoder_->get_bitrate();
	this->codec_cxt_->gop_size = this->encoder_->get_keyframe_interval();
	this->codec_cxt_->pix_fmt = PIX_FMT_YUV420P;
	
	// this is all framerate stuff, which doesn't seem to affect what we are doing
	this->codec_cxt_->time_base.num = 1;
	this->codec_cxt_->time_base.den = 25;
	
	// Image size this codec is processing
	this->codec_cxt_->width = this->image_width_;
	this->codec_cxt_->height = this->image_height_;
	
	// Allocate the RGB image.  Because the buffer will be input as a parameter to
	// the encode_image function, we do not need to pass in a buffer
	this->frame_rgb_ = avcodec_alloc_frame();
	if ( !this->frame_rgb_ )
	{
		CORE_LOG_ERROR( "Could not allocate RGB frame for MPEG4 codec." );
		return false;
	}
	
	// Allocate the YUV image.
	this->yuv_buf_ = Buffer::New( avpicture_get_size( this->codec_cxt_->pix_fmt, 
													  this->image_width_, this->image_height_ ) );
	if ( !this->yuv_buf_ )
	{
		CORE_LOG_ERROR( "Could not allocate YUV buffer object for MPEG4 codec." );
		return false;
	}
	// Allocate the YUV picture.
	this->frame_yuv_ = avcodec_alloc_frame();
	
	if ( !this->frame_yuv_ )
	{
		CORE_LOG_ERROR( "Could not allocate YUV frame for MPEG4 codec." );
		return false;
	}
	
	avpicture_fill( reinterpret_cast<AVPicture*>( this->frame_yuv_ ), 
				    static_cast<boost::uint8_t*>( this->yuv_buf_->get_data_ptr() ),
				    PIX_FMT_YUV420P, this->image_width_, this->image_height_ );
	
	
	this->codec_ = avcodec_find_encoder( this->codec_cxt_->codec_id );
	if ( this->codec_ == 0 )
	{
		// We could not find the encoder
		CORE_LOG_ERROR( "Could not find FFMPEG codec for MPEG4." );
		return false;
	}
	
	if ( avcodec_open( this->codec_cxt_, this->codec_ ) < 0 ) 
	{
		// Cannot open Codec
		CORE_LOG_ERROR( "Could not open FFMPEG codec for MPEG4." );
		return false;
	}
	
	this->img_convert_ctx_ = sws_getContext( this->image_width_, this->image_height_, 
											 PIX_FMT_RGB24, 
											 this->codec_cxt_->width, 
											 this->codec_cxt_->height, 
											 PIX_FMT_YUV420P, SWS_BICUBIC, 
											 NULL, NULL, NULL );
	
	this->encoder_->update_codec_serial_id();
	this->encoder_->reset_parameters_changed();

	this->is_initialized_ = true;
	
	return true;
}
	
void MPEGEncoderPrivate::clear_codec()
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

//------ MPEGEncoder class

MPEGEncoder::MPEGEncoder() :
	Encoder( MPEG_CODEC_TYPE_C ), 
	private_( new MPEGEncoderPrivate )
{
	// Fill in the pointer back to this class
	this->private_->encoder_ = this;
	
}
	
MPEGEncoder::~MPEGEncoder() {}

bool MPEGEncoder::encode_image( Image3Handle input, CompressedImage3Handle& output )
{
	lock_type lock( this->get_mutex() );

	// At the outset we free the output object.
	output.reset();

	if ( !this->private_->is_initialized_ || 
		 input->get_width() != this->private_->image_width_ || 
		 input->get_height() != this-> private_->image_height_ || 
		 this->get_parameters_changed() ) 
	{
		if (! this->private_->initialize_codec( input->get_width(), input->get_height() ) )
		{
			CORE_LOG_ERROR( "Failed to encode image." );
			return false;
		}
	}
	
	// Put the input rgb image into the AV Picture data structure.
	avpicture_fill( reinterpret_cast<AVPicture*>( this->private_->frame_rgb_ ), 
				    static_cast<boost::uint8_t*>( input->get_buffer()->get_data_ptr() ),
				    PIX_FMT_RGB24, input->get_width(), input->get_height() );

	// Convert the RGB image to YUV.
	sws_scale( this->private_->img_convert_ctx_, this->private_->frame_rgb_->data, 
			   this->private_->frame_rgb_->linesize, 0, 
			   this->private_->codec_cxt_->height, this->private_->frame_yuv_->data, 
			   this->private_->frame_yuv_->linesize );


	// A guess at the maximum size of the compressed image
	size_t output_buffer_size = std::max( static_cast<size_t>( 4096 ), 
										  3 * input->get_width() * input->get_height() );
	BufferHandle output_buffer = Buffer::New( output_buffer_size );

	if (! output_buffer )
	{
		CORE_LOG_ERROR( "Failed to allocate output buffer." );
		return false;
	}
	
	output = CompressedImage3::New( input->get_width(), input->get_height(), 
								    output_buffer );
	if (! output )
	{
		CORE_LOG_ERROR( "Failed to create output image." );
		return false;
	}

	int encoded_size = avcodec_encode_video( this->private_->codec_cxt_, 
											 static_cast<boost::uint8_t*>( output_buffer->get_data_ptr() ), 
											 output_buffer->get_data_size(),
											 this->private_->frame_yuv_ );	
																				 
	// if encode fails, one possibility is that the allocated buffer was not large enough
	// try again with one twice as large
	if ( encoded_size < 0 )
	{
		output_buffer = Buffer::New( 2*output_buffer_size );
		if (! output_buffer )
		{
			output.reset();
			CORE_LOG_ERROR( "Failed to allocate output buffer." );
			return false;
		}
		output = CompressedImage3::New( input->get_width(), input->get_height(),
									    output_buffer );
		
		encoded_size = avcodec_encode_video( this->private_->codec_cxt_, 
											 static_cast<boost::uint8_t*>( output_buffer->get_data_ptr() ), 
											 output_buffer->get_data_size(),
											 this->private_->frame_yuv_ );	
	}
	
	if ( encoded_size < 0 )
	{
		output.reset();
		CORE_LOG_ERROR( "Failed to encode video with avcodec_encode_video." );
		return false;
	}
	
	bool is_keyframe = ( this->private_->codec_cxt_->coded_frame->key_frame != 0 );

	output->set_is_keyframe( is_keyframe );
	output->set_codec_type( this->get_codec_type() );
	output->set_codec_serial_id( this->get_codec_serial_id() );
	
	return true;
}

} // end namespace Core
