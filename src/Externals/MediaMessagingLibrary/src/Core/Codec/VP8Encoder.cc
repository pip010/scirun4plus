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

#include <algorithm>
#include <exception>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include <boost/cstdint.hpp>
#include <boost/lexical_cast.hpp>


#include <Core/Codec/VP8Encoder.h>
#include <Core/Utils/Buffer.h>
#include <Core/Utils/Log.h>

// FFMPEG includes
extern "C" {
	
#include <libswscale/swscale.h>

}

#define VPX_INTEGER_H
#include <vpx/vpx_encoder.h>
#include <vpx/vp8cx.h>

namespace Core
{
	
class VP8EncoderPrivate
{

	// Member variables
public:
	// Pointer back to the main class
	VP8Encoder* encoder_;
	
	bool is_initialized_;

	// The size of the image for which we initialize the codec
	unsigned int image_width_;
	unsigned int image_height_;
	
	vpx_codec_ctx_t* codec_ctx_;
	vpx_codec_enc_cfg_t* codec_cfg_;
	vpx_codec_pts_t frame_count_;
	vpx_image_t* rgb_img_;
	vpx_image_t* yuv_img_;
	BufferHandle yuv_buffer_;
	SwsContext* img_convert_ctx_;
	
	// Private functions
public:	
	VP8EncoderPrivate() :
		encoder_( 0 ),
		is_initialized_( false ),
		image_width_( 0 ),
		image_height_( 0 ),
		codec_ctx_( 0 ),
		codec_cfg_( 0 ),
		frame_count_( 0 ),
		rgb_img_( 0 ),
		yuv_img_( 0 ),
		img_convert_ctx_( 0 )
	{
	}
	
	~VP8EncoderPrivate()
	{
		clear_codec();
	}
	
	bool initialize_codec( unsigned int image_width, unsigned int image_height );
	
	void clear_codec();
	
};

bool VP8EncoderPrivate::initialize_codec( unsigned int image_width, unsigned int image_height )
{
	this->clear_codec();
	
	// Set the size for the new image series
	this->image_width_ = image_width;
	this->image_height_ = image_height;
	
	this->rgb_img_ = new vpx_image_t;
	
	vpx_codec_err_t codec_err;
	
	// Allocate the encoder configuration file.
	codec_cfg_ = new vpx_codec_enc_cfg_t;
	
	// Plugin the default VP8 configuration.
	codec_err = vpx_codec_enc_config_default( vpx_codec_vp8_cx(), 
											  codec_cfg_, 0 );
	if ( codec_err != 0 )
	{
		CORE_LOG_ERROR( "Failed to initialize VP8 with default configuration." );
		return false;
	}
	
	codec_cfg_->rc_end_usage = VPX_VBR;
	codec_cfg_->g_pass = VPX_RC_ONE_PASS;
	codec_cfg_->rc_target_bitrate = static_cast< unsigned int >( encoder_->get_bitrate() / 1000 );
	codec_cfg_->g_w = static_cast< unsigned int >( image_width );
	codec_cfg_->g_h = static_cast< unsigned int >( image_height );
	codec_cfg_->kf_max_dist = static_cast< unsigned int >( encoder_->get_keyframe_interval() );
	codec_cfg_->rc_resize_allowed = 0;
	
	// Make a new context
	try {
		codec_ctx_ = new vpx_codec_ctx_t;
	} catch ( ... ) {
		// We could not allocate the context
		CORE_LOG_ERROR( "Failed to allocate VP8 context." );
		return false;
	}
	codec_err = vpx_codec_enc_init( codec_ctx_, vpx_codec_vp8_cx(),
								    codec_cfg_, 0 );
	
	if ( codec_err != 0 )
	{
		// We could not initialize the context
		CORE_LOG_ERROR( "Failed to initialize VP8 context." );
		return false;
	}
	
	// I420 is a 12 bit per pixel format.  It consists of a full-resolution
	// sampling of the Y channel, and a half-resolution subsampling of the U
	// and V channels.
	// Because I don't want to depend on the internal representation of I420,
	// and because the VP8 library -- unlike the FFMPEG library -- does not 
	// have a function for determining the amount of space required, I allocate
	// a very conservative amount -- 8 bits per channel -- to cover my bases.
	this->yuv_buffer_ = Buffer::New( image_width * image_height * 4 );
	if ( !this->yuv_buffer_ )
	{
		CORE_LOG_ERROR( "Could not allocate YUV buffer object for VP8 codec." );
		return false;
	}
	
	yuv_img_ = vpx_img_wrap( NULL, VPX_IMG_FMT_I420, image_width, image_height, 1, 
							 reinterpret_cast<unsigned char*>( this->yuv_buffer_->get_data_ptr() ) );
	
	if ( !this->yuv_img_ )
	{
		CORE_LOG_ERROR( "Could not allocate YUV image for VP8 codec." );
		return false;
	}
		
	this->img_convert_ctx_ = sws_getContext( this->image_width_, this->image_height_, 
											 PIX_FMT_RGB24, 
											 this->image_width_, this->image_height_,
											 PIX_FMT_YUV420P, SWS_BILINEAR, 
											 NULL, NULL, NULL );
	
	this->encoder_->update_codec_serial_id();
	this->encoder_->reset_parameters_changed();
	this->frame_count_ = 0;

	this->is_initialized_ = true;
	
	return true;
}
	
void VP8EncoderPrivate::clear_codec()
{
	if ( !this->is_initialized_ ) return;
	
	vpx_codec_destroy( this->codec_ctx_ );
	delete( this->codec_ctx_ );
	this->codec_ctx_ = 0;
	
	delete( this->codec_cfg_ );
	this->codec_cfg_ = 0;
	
	vpx_img_free( this->rgb_img_ );
	this->rgb_img_ = 0;
	vpx_img_free( this->yuv_img_ );
	this->yuv_img_ = 0;
	this->yuv_buffer_.reset();
	
	sws_freeContext( this->img_convert_ctx_ ); 
	this->img_convert_ctx_ = 0;
	
	delete( this->rgb_img_ );
	
	this->is_initialized_ = false;
}

//------ VP8Encoder class

VP8Encoder::VP8Encoder() :
	Encoder( VP8_CODEC_TYPE_C ), 
	private_( new VP8EncoderPrivate )
{
	// Fill in the pointer back to this class
	this->private_->encoder_ = this;	
}

VP8Encoder::~VP8Encoder() {}
	
bool VP8Encoder::encode_image( Image3Handle input, CompressedImage3Handle& output )
{
	lock_type lock( this->get_mutex() );

	// At the outset we free the output object.
	output.reset();

	if ( !this->private_->is_initialized_ || 
		 input->get_width() != this->private_->image_width_ || 
		 input->get_height() != this->private_->image_height_ || 
		 this->get_parameters_changed() ) 
	{
	
		if (! this->private_->initialize_codec( input->get_width(), input->get_height() ) )
		{
			CORE_LOG_ERROR( "Failed to initialize codec." );
			return false;
		}
	}
	
	// wrap the input buffer's data in a vpx image.
	vpx_img_wrap( this->private_->rgb_img_, VPX_IMG_FMT_RGB24, input->get_width(), 
				  input->get_height(), 1, 
				  reinterpret_cast<unsigned char*>( input->get_buffer()->get_data_ptr() ) );

	// Convert the RGB image to YUV.
	sws_scale( this->private_->img_convert_ctx_, this->private_->rgb_img_->planes, 
			   this->private_->rgb_img_->stride, 0, 
			   input->get_height(), this->private_->yuv_img_->planes, 
			   this->private_->yuv_img_->stride );

	vpx_codec_err_t err = vpx_codec_encode( this->private_->codec_ctx_, 
										    this->private_->yuv_img_, 
										    this->private_->frame_count_, 
										    1, 0, VPX_DL_REALTIME );
	if ( err != 0 )
	{
		output.reset();
		CORE_LOG_ERROR("Failed to encode video stream.");
		return false;
	}
	
	++this->private_->frame_count_;
	vpx_codec_iter_t iter = NULL;
	const vpx_codec_cx_pkt_t *pkt;

	while ( ( pkt = vpx_codec_get_cx_data( this->private_->codec_ctx_, &iter ) ) != 0 )
	{
		switch ( pkt->kind )
		{
			case VPX_CODEC_CX_FRAME_PKT:
			{

				BufferHandle output_buffer = Buffer::New( pkt->data.frame.sz );
				
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
				
				memcpy( output_buffer->get_data_ptr(), pkt->data.frame.buf, pkt->data.frame.sz );
				
				bool is_keyframe = ( pkt->data.frame.flags & VPX_FRAME_IS_KEY );
				
				output->set_is_keyframe( is_keyframe );
				output->set_codec_type( this->get_codec_type() );
				output->set_codec_serial_id( this->get_codec_serial_id() );
				
				break;
			}
			default:
				break;
		}
	}
	
	return true;
}

} // end namespace Core
