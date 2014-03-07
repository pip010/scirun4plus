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

#include <Core/Codec/Image3.h>
#include <string.h>
#include <vector>

namespace Core {

Image3::Image3( size_t width, size_t height, BufferHandle buf ) :
	width_( width ), height_( height ), buffer_( buf )
{}
	
Image3::~Image3()
{}

BufferHandle Image3::get_buffer() const
{
	return this->buffer_;
}

unsigned char* Image3::get_image_data() const
{
	return reinterpret_cast< unsigned char *>( this->buffer_->get_data_ptr() );
}


size_t Image3::get_width() const
{
	return this->width_;
}

size_t Image3::get_height() const
{
	return this->height_;
}

bool Image3::flip_vertical()
{
	std::vector<unsigned char> tmp_row_buffer( 3 * this->width_ );
	unsigned char *top_row, *bot_row;     
	unsigned char *tmp_row = &tmp_row_buffer[ 0 ];

	size_t top, bot;
	for ( top = this->height_ - 1, bot = 0; bot < this->height_/2; top--, bot++){
		top_row = reinterpret_cast<unsigned char *>( this->buffer_->get_data_ptr() ) + this->width_ * top * 3;
		bot_row = reinterpret_cast<unsigned char *>( this->buffer_->get_data_ptr() ) + this->width_ * bot * 3;
		memcpy( tmp_row, top_row, this->width_ * 3 );
		memcpy( top_row, bot_row, this->width_ * 3 );
		memcpy( bot_row, tmp_row, this->width_ * 3 );
	}
	
	return true;
}
	
Image3Handle Image3::New(size_t width, size_t height, BufferHandle buf)
{
	try {
		return Image3Handle( new Image3( width, height, buf ) );
	} catch ( ... ) {
		return Image3Handle();
	}
}

	
} // namespace Core
