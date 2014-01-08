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

#include <iostream>
#include <fstream>
#include <stdio.h>

#include <openssl/rsa.h>
#include <openssl/pem.h>

#include "Message.h"

namespace BioMesh3d
{

Message::Message() :
	m_message_type_( MessageTypes::INVALID_E )
{
}


Message::~Message()
{
}


void Message::reset_message()
{
	boost::mutex::scoped_lock lock( this->m_message_values_mutex_ );
	m_message_values_.clear();
	set_message_type( MessageTypes::INVALID_E );
}


void Message::set_binary_size( boost::uint32_t bin_size ) 
{
	if ( !this->buffer_ ) this->buffer_ = Core::Buffer::New( static_cast<size_t>( bin_size ) );
	this->buffer_->resize( static_cast<size_t>( bin_size ) );
}


void Message::copy_binary_data( boost::uint8_t* data, size_t data_size )
{
	if ( !this->buffer_ ) this->buffer_ = Core::Buffer::New( static_cast<size_t>( data_size ) );
	this->buffer_->resize( static_cast<size_t>( data_size ) );

	memcpy( this->buffer_->get_data_ptr() , data, data_size );
}

boost::uint8_t* Message::get_binary_buf()
{
	if ( this->buffer_ )
	{
		return static_cast<boost::uint8_t* >( this->buffer_->get_data_ptr() );
	}
	
	return 0;
}

boost::uint32_t Message::get_binary_size()
{
	if ( this->buffer_ )
	{
		return static_cast<boost::uint32_t>( this->buffer_->get_data_size() );
	}
	
	return 0;
}

MessageTypes::MessageType Message::get_message_type() 
{ 
	return this->m_message_type_; 
}

void Message::set_message_type( MessageTypes::MessageType message_type ) 
{ 
	this->m_message_type_ = message_type; 
}


Core::BufferHandle Message::get_binary_buffer()
{
	return this->buffer_;
}

void Message::set_binary_buffer( Core::BufferHandle buffer )
{
	this->buffer_ = buffer;
}


std::string Message::serialize_parameters()
{
	boost::mutex::scoped_lock lock( m_message_values_mutex_ );

	std::string serialized_string = "";
	std::map< int, std::string >::iterator iter;
	for ( iter = this->m_message_values_.begin(); iter != this->m_message_values_.end(); ++iter )
	{
		try
		{
			serialized_string += boost::lexical_cast< std::string >( iter->first ) + "=" + iter->second + "&";
		}
		catch ( std::exception& exc )
		{
			std::cout << "Lexical cast error: " << exc.what() << std::endl;
		}
	}

	std::string::size_type pos = serialized_string.find_last_of("&");
	if ( pos != std::string::npos ) // if the string has data in it then remove the last ampersand and copy it 
	{								// to m_raw_message_
		serialized_string.erase( pos );
	}
	return serialized_string;
}

}