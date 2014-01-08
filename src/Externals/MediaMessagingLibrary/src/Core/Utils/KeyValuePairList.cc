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

// STL includes
#include <fstream>
#include <exception>

// COre includes
#include <Core/Utils/KeyValuePairList.h>
#include <Core/Utils/Log.h>
#include <Core/Utils/StringUtil.h>
#include <Core/Utils/StringParser.h>

namespace Core
{


KeyValuePairList::KeyValuePairList() :
	fail_on_new_key_( false )
{
}


KeyValuePairList::KeyValuePairList( const std::string& key_value_list ) :
	fail_on_new_key_( false )
{
	if ( ! this->import_from_string( key_value_list ) )
	{
		throw std::logic_error( "Could not convert string to key/value pair list." );
	}
}


void KeyValuePairList::fail_on_new_key( bool enable )
{
	this->fail_on_new_key_ = enable;
}


std::string KeyValuePairList::export_to_string() const
{
	std::string key_value_string;
	std::map<std::string,std::string>::const_iterator it = this->key_value_.begin();
	std::map<std::string,std::string>::const_iterator it_end = this->key_value_.end();
	
	while ( it != it_end )
	{
		key_value_string += (*it).first + '=' + (*it).second + "\r\n";
		++it;
	}
	
	return key_value_string;
}


bool KeyValuePairList::import_from_string( const std::string& key_value_string )
{
	std::string error;
	return this->import_from_string( key_value_string, error );
}


bool KeyValuePairList::import_from_string( const std::string& key_value_string, std::string& error )
{
	std::string::size_type start = 0;
	std::string key;
	std::string value;
	
	while ( start < key_value_string.size() )
	{
		if (! ScanKeyValuePair( key_value_string, start, key, value, error ) )
		{
			this->key_value_.clear();
			return false;
		}
		
		if ( !key.empty() )
		{
			if ( this->fail_on_new_key_ )
			{
				std::map<std::string,std::string>::iterator it = this->key_value_.find( key );
				if ( this->key_value_.end() == it )
				{
					error = std::string( "Could not parse input as key '"  ) + key + 
						"' is not known.";
					return false;
				}
				else
				{	
					// Since we already did the lookup reuse it
					(*it).second = value;
					continue;
				}
				
			}
			this->key_value_[ key ] = value;		
		}
	}
	
	return true;
}


bool KeyValuePairList::import_from_string( const std::string& key_value_string,
	std::set<std::string>& changed_keys )
{
	std::string error;
	return this->import_from_string( key_value_string, changed_keys, error );
}


bool KeyValuePairList::import_from_string( const std::string& key_value_string, 
	std::set<std::string>& changed_keys, std::string& error )
{
	std::string::size_type start = 0;
	std::string key;
	std::string value;
	
	while ( start < key_value_string.size() )
	{
		if (! ScanKeyValuePair( key_value_string, start, key, value, error ) )
		{
			this->key_value_.clear();
			return false;
		}
		
		if ( !key.empty() )
		{
			std::map<std::string,std::string>::iterator it = this->key_value_.find( key );

			if ( this->fail_on_new_key_ )
			{
				if ( this->key_value_.end() == it )
				{
					error = std::string( "Could not parse input as key '"  ) + key + 
						"' is not known.";
					return false;
				}
			}
			
			if ( it != this->key_value_.end() )
			{
				// Only set it if it was changed
				if ( value != (*it).second )
				{
					changed_keys.insert( key );
					(*it).second = value;	
				}
			}
			else
			{
				// It's a new key
				this->key_value_[ key ] = value;		
				changed_keys.insert( key );
			}
		}
	}
	
	return true;
}


bool KeyValuePairList::export_to_file( const boost::filesystem::path& filename ) const
{
	std::string error;
	return this->export_to_file( filename, error );
}


bool KeyValuePairList::export_to_file( const boost::filesystem::path& filename, std::string& error ) const
{
	error.clear();
	try
	{
		std::ofstream file( filename.string().c_str() );
		file << this->export_to_string();
		return true;
	}
	catch ( ... )
	{
		error = std::string( "Could not open file '" ) + filename.string() + "'.";
		return false;
	}
}


bool KeyValuePairList::import_from_file( const boost::filesystem::path& filename )
{
	std::string error;
	return this->import_from_file( filename, error );
}
	
	
bool KeyValuePairList::import_from_file( const boost::filesystem::path& filename, std::string& error )
{
	error.clear();
	
	if (! boost::filesystem::exists( filename ) )
	{
		error = std::string( "File '" ) + filename.string() + "' does not exist.";
		return false;
	}
	
	try
	{
		std::ifstream iss( filename.string().c_str() );

		while ( !iss.eof() )
		{
			std::string line;
			std::getline( iss, line );
				
			std::string::size_type start = 0;
			std::string key;
			std::string value;
			std::string error;
			
			while ( start < line.size() )
			{
				if (! ScanKeyValuePair( line, start, key, value, error ) )
				{
					this->key_value_.clear();
					return false;
				}
				
				if ( !key.empty() )
				{
					if ( this->fail_on_new_key_ )
					{
						std::map<std::string,std::string>::iterator it = this->key_value_.find( key );
						if ( this->key_value_.end() == it )
						{
							error = std::string( "Could not parse input as key '"  ) + key + 
								"' is not known.";
							return false;
						}						
					}
					this->key_value_[ key ] = value;		
				}
			}
		}
	}
	catch ( ... )
	{
		error = std::string( "Could not open file '" ) + filename.string() + "'.";
		return false;
	}

	return true;
}

bool KeyValuePairList::import_from_arguments( int argc, char** argv )
{
	std::string error;
	return this->import_from_arguments( argc, argv, error );
}


bool KeyValuePairList::import_from_arguments( int argc, char** argv, std::string& error )
{
	std::string argument_string;

	for ( int j = 1; j < argc; j++ )
	{
		argument_string += std::string( argv[ j ] ) + " ";
	}
	
	return this->import_from_string( argument_string, error );
}


bool KeyValuePairList::get_value( const std::string& key, std::string& value ) const
{
	std::map<std::string,std::string>::const_iterator it = this->key_value_.find( key );
	if ( it == this->key_value_.end() )
	{
		return false;
	}
	
	value = (*it).second;
	return true;
}
	
bool KeyValuePairList::set_value( const std::string& key, const std::string& value )
{
	this->key_value_[ key ] = value;
	return true;
}


void KeyValuePairList::clear()
{
	this->key_value_.clear();
	this->fail_on_new_key_ = false;
}

} // end namespace Core
