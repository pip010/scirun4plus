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

extern "C" 
{
#include <png.h>
}

#include <Core/Codec/PNGEncoder.h>
#include <Core/Utils/Log.h>

namespace Core
{
	
typedef png_uint_32 size_type;

class PNGEncoderPrivate
{
public:

	PNGEncoderPrivate() :
		png_ptr_( NULL ), 
		info_ptr_( NULL ),
		row_pointers_( NULL ),
		image_width_( 0 ),
		image_height_( 0 ),
		is_initialized_( false )
	{
	}
	
	~PNGEncoderPrivate()
	{
		clear_codec();
	}
	
	png_structp png_ptr_;
	png_infop info_ptr_;
	png_byte** row_pointers_;
	size_type image_width_;
	size_type image_height_;
	
	bool is_initialized_;
	
	bool initialize_codec( size_type image_width, size_type image_height );	
	void clear_codec();

	
	static void my_png_write_data( png_structp png_ptr, png_bytep data, png_size_t length );
	
}; // class PNGEncoderPrivate

void PNGEncoderPrivate::clear_codec()
{
	if (! this->is_initialized_ ) 
		return;
	
	png_destroy_write_struct( &this->png_ptr_, &this->info_ptr_ );
	delete [] this->row_pointers_;
	this->row_pointers_ = 0;
	this->image_width_	= 0;
	this->image_height_ = 0;

	this->is_initialized_ = false;
}

bool PNGEncoderPrivate::initialize_codec( size_type image_width, size_type image_height )
{
	this->clear_codec();
	
	this->image_width_  = image_width;
	this->image_height_ = image_height;
	
	this->png_ptr_ = png_create_write_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
	if ( this->png_ptr_ == NULL )
	{
		CORE_LOG_ERROR( "Failed to create write struct." );
		return false;
	}
	
	this->info_ptr_ = png_create_info_struct( this->png_ptr_ );
	if ( this->info_ptr_ == NULL )
	{
		png_destroy_write_struct( &this->png_ptr_, NULL );
		CORE_LOG_ERROR( "Failed to create info struct." );
		return false;
	}
	
	// Set the return point in case png_set_IHDR fails below.
	if ( setjmp( png_jmpbuf( this->png_ptr_ ) ) ) {
		png_destroy_write_struct( &this->png_ptr_, &this->info_ptr_ );
		CORE_LOG_ERROR( "Failed to set IHDR for PNG." ); 
		return false;
	}
	
	png_set_IHDR( this->png_ptr_, this->info_ptr_, image_width, image_height,
				  8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
				  PNG_FILTER_TYPE_DEFAULT );

	row_pointers_ = new png_byte*[image_height];
	
	// Best compression....
	this->png_ptr_->compression = 9;
	
	png_set_rows( this->png_ptr_, this->info_ptr_, this->row_pointers_ );
	
	this->is_initialized_ = true;
	
	return true;
}

// static function passed into the PNG library specifying how to write the results.
void PNGEncoderPrivate::my_png_write_data( png_structp png_ptr, png_bytep data, png_size_t length )
{
	
	BufferHandle* buffer = reinterpret_cast<BufferHandle *>( png_ptr->io_ptr );
	size_t old_size = 0;
	
	if (! *buffer )
	{
		// allocate buffer for the first time.
		*buffer = Buffer::New( length );
		if (! *buffer )
		{
			// indicates a probable out-of-memory condition
			CORE_LOG_ERROR( "Unable to create Buffer." );
			// throw the PNG library error.
			png_error( png_ptr, "Unable to create buffer for PNG." );
			
			return; // control should never get here.
		}
	}
	else
	{
		old_size = ( *buffer )->get_data_size();
		size_t new_size = old_size + length;
		
		if (! ( *buffer )->resize( new_size ) )
		{
			// if we were unable to resize the Buffer, log an error
			// and throw the PNG library error.
			CORE_LOG_ERROR( "Unable to resize Buffer." );
			png_error( png_ptr, "Unable to resize buffer to accomodate PNG." );
			
			return; // control should never get here.
		}
	}
	
	char* buf = reinterpret_cast<char*>( ( *buffer )->get_data_ptr() );
	
	// append the new data to the end of the old buffer.
	memcpy( buf + old_size, data, length );
}

// Begin PNGEncoder

PNGEncoder::PNGEncoder() : 
	Encoder( PNG_CODEC_TYPE_C ),
	private_( new PNGEncoderPrivate )
{
}

PNGEncoder::~PNGEncoder()
{
}

bool PNGEncoder::encode_image( Image3Handle input, CompressedImage3Handle& output )
{
	lock_type lock( this->get_mutex() );

	// invalidate output
	output.reset();
	
	// if our encoder is not initialized or the input parameters have changed,
//	// we must initialize.
//	if ( !this->private_->is_initialized_ || 
//		 input->get_width() != this->private_->image_width_ || 
//		 input->get_height() != this-> private_->image_height_ || 
//		 this->get_parameters_changed() ) 
//	{
		if ( !this->private_->initialize_codec( static_cast< size_type >( input->get_width() ), 
			static_cast< size_type >( input->get_height() ) ) )
		{
			CORE_LOG_ERROR( "Failed to initialize PNG Codec." );
			return false;
		}
//	}
	
	// Make the png data structure row pointers point into
	// our input data structure.  *3 because it is an image
	// with 3 components (RGB)
	unsigned int row_size = static_cast< size_type >( input->get_width() * 3 );
	void* in_buf = input->get_buffer()->get_data_ptr();
	for ( unsigned int i = 0; i < input->get_height(); ++i )
	{
		this->private_->row_pointers_[i] = 
			&reinterpret_cast<png_byte*>( in_buf )[ i * row_size ];
	}

	BufferHandle output_buffer;
	
	// Set up the error condition return point for the png_write_png function call.
	// Log an error, cleanup, exit.
	if ( setjmp( png_jmpbuf( this->private_->png_ptr_ ) ) )
	{
		CORE_LOG_ERROR( "Failed to encode PNG." ); 
		return false;
	}

	png_set_write_fn( this->private_->png_ptr_, &output_buffer, 
					  &PNGEncoderPrivate::my_png_write_data, 
					  NULL );
	
	png_write_png( this->private_->png_ptr_, this->private_->info_ptr_, PNG_TRANSFORM_IDENTITY, 
				   NULL );
	
	// Put the buffer with the encoded PNG into the output image.
	output = CompressedImage3::New( input->get_width(), input->get_height(), output_buffer );
	output->set_is_keyframe( true );
	output->set_codec_type( this->get_codec_type() );
	output->set_codec_serial_id( this->get_codec_serial_id() );
	
	return true;

}

} // end namespace Core