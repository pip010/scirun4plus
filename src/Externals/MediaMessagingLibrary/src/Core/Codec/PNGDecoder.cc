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

#include <vector>

extern "C" 
{
#include <png.h>
}

#include <Core/Codec/PNGDecoder.h>
#include <Core/Codec/PNGEncoder.h>
#include <Core/Math/MathFunctions.h>
#include <Core/Utils/Log.h>

namespace Core
{

class StatefulBuffer
{
public:
	BufferHandle buffer_;
	size_t current_offset_;
};

typedef png_uint_32 size_type;

class PNGDecoderPrivate
{
public:
	PNGDecoderPrivate() :
		png_ptr_( NULL ),
		info_ptr_( NULL )
	{
	}
	
	~PNGDecoderPrivate()
	{
		this->clear_codec();
	}

	void clear_codec();
	bool initialize_codec();

	png_structp png_ptr_;
	png_infop info_ptr_;
	std::vector< png_bytep > row_ptrs_;
};

void PNGDecoderPrivate::clear_codec()
{	
	if ( this->png_ptr_ )
	{
		png_destroy_read_struct( &this->png_ptr_, &this->info_ptr_, NULL );
	}
	this->png_ptr_ = NULL;
	this->info_ptr_ = NULL;
}
	
bool PNGDecoderPrivate::initialize_codec()
{
	clear_codec();
	
	this->png_ptr_ = png_create_read_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
	if ( !this->png_ptr_ )
	{
		return false;
	}

	this->info_ptr_ = png_create_info_struct( this->png_ptr_ );
	if ( !this->info_ptr_ )
	{
		png_destroy_read_struct( &this->png_ptr_, NULL, NULL );
		return false;
	}

	return true;
}

static void ReadPNGData( png_structp png_ptr, png_bytep data, png_size_t len )
{
	StatefulBuffer* buffer = reinterpret_cast< StatefulBuffer* >( png_get_io_ptr( png_ptr ) );
	png_bytep src_data = reinterpret_cast< png_bytep >( buffer->buffer_->get_data_ptr() );
	size_t read_len = Core::Min( buffer->buffer_->get_data_size() - buffer->current_offset_, len );
	memcpy( data, src_data + buffer->current_offset_, read_len );
	buffer->current_offset_ += read_len;
}

PNGDecoder::PNGDecoder() :
	Decoder( PNG_CODEC_TYPE_C ),
	private_( new PNGDecoderPrivate )
{
}

PNGDecoder::~PNGDecoder()
{
}

bool PNGDecoder::decode_image( CompressedImage3Handle input, Image3Handle& output )
{
	lock_type lock( this->get_mutex() );

	void* data_ptr = input->get_buffer()->get_data_ptr();

	// Validate the input to make sure it's really PNG
	if ( input->get_buffer()->get_data_size() < 8 ||
		png_sig_cmp( reinterpret_cast< png_bytep >( data_ptr ), 0, 8 ) != 0 )
	{
		return false;
	}

	
	if ( !this->private_->initialize_codec() ) return false;

	if ( setjmp( png_jmpbuf( this->private_->png_ptr_ ) ) )
	{
		CORE_LOG_ERROR( "Failed to decode PNG." ); 
		return false;
	}

	StatefulBuffer read_buffer;
	read_buffer.buffer_ = input->get_buffer();
	read_buffer.current_offset_ = 0;

	png_set_read_fn( this->private_->png_ptr_, &read_buffer, ReadPNGData );

	// Read the PNG header
	png_read_info( this->private_->png_ptr_, this->private_->info_ptr_ );
	size_type input_width = png_get_image_width( this->private_->png_ptr_, this->private_->info_ptr_ );
	size_type input_height = png_get_image_height( this->private_->png_ptr_, this->private_->info_ptr_ );
	size_type color_type = png_get_color_type( this->private_->png_ptr_, this->private_->info_ptr_ );
	size_type bit_depth = png_get_bit_depth( this->private_->png_ptr_, this->private_->info_ptr_ );
	// Only deal with RGB24 format for now
	if ( color_type != PNG_COLOR_TYPE_RGB || bit_depth != 8 )
	{
		return false;
	}

	// Allocate enough row pointers
	this->private_->row_ptrs_.resize( input_height );
	// Allocate the output buffer
	BufferHandle out_buffer = Buffer::New( input_width * input_height * 3 );
	png_bytep dst_ptr = reinterpret_cast< png_bytep >( out_buffer->get_data_ptr() );
	size_type row_size = input_width * 3;
	for ( size_type i = 0; i < input_height; ++i )
	{
		this->private_->row_ptrs_[ i ] = dst_ptr + row_size * i;
	}

	png_read_image( this->private_->png_ptr_, &this->private_->row_ptrs_[ 0 ] );
	png_read_end( this->private_->png_ptr_, NULL );
	output = Image3::New( input_width, input_height, out_buffer );

	return true;
}

} // namespace Core