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

// STL inlcudes
#include <map>
#include <fstream>

// Boost includes
#include <boost/filesystem.hpp>

// Core includes
#include <Core/Utils/Lockable.h>
#include <Core/Utils/Config.h>
#include <Core/Utils/StringUtil.h>

namespace Core
{

class ConfigPrivate : public RecursiveLockable
{
public:
	// Map that contains the key value pairs for the configuration file
	std::map<std::string,std::string> params_;
};
	
Config::Config() :
	private_( new ConfigPrivate )
{
}

Config::~Config()
{
}

bool Config::set_defaults( const std::string& defaults )
{
	// Lock internals of the class
	ConfigPrivate::lock_type lock( this->private_->get_mutex() );
	
	std::vector<std::string> values = SplitString( defaults, "|" );
	
	for ( size_t j = 0; j <  values.size(); j++ )
	{
		std::string::size_type idx = values[ j ].find( "=" );
		// ignore entries without equal signs
		if ( idx == std::string::npos ) continue;
		std::string key = values[ j ].substr( 0, idx );
		std::string value = values[ j ].substr( idx + 1 );
		
		StripSurroundingSpaces( key );
		StripSurroundingSpaces( value );
		if ( key.empty() ) continue;
		// ignore comments
		if ( key[ 0 ] == '%' || key[ 0 ] == '/' || key[ 0 ] == '#' ) continue;
		
		key = StringToLower( key );
		this->private_->params_[ key ] = value;
	}
	
	return true;
}
	
bool Config::load_config_file( const std::string& config_file )
{
	// Lock internals of the class
	ConfigPrivate::lock_type lock( this->private_->get_mutex() );

	boost::filesystem::path cfile( config_file );
	try
	{
		cfile = boost::filesystem::absolute( cfile );
		if ( ! boost::filesystem::exists( cfile ) ) return false;
	}
	catch ( ... )
	{
		// File does not exist
		return false;
	}

	try
	{
		std::ifstream input_file( cfile.string().c_str() );

		while ( ! input_file.eof() )
		{
			std::string line;
			std::getline( input_file, line );
			
			std::string::size_type idx = line.find( "=" );
			// ignore entries without equal signs
			if ( idx == std::string::npos ) continue;
			std::string key = line.substr( 0, idx );
			std::string value = line.substr( idx + 1 );
			
			StripSurroundingSpaces( key );
			StripSurroundingSpaces( value );
			if ( key.empty() ) continue;
			// ignore comments
			if ( key[ 0 ] == '%' || key[ 0 ] == '/' || key[ 0 ] == '#' ) continue;
			
			key = StringToLower( key );
			this->private_->params_[ key ] = value;
		}
	}
	catch ( ... )
	{
		// Could not read file
		return false;
	}
	
	return true;
}
	

std::string Config::get_param( const std::string& key ) const
{	
	// Lock internals of the class
	ConfigPrivate::lock_type lock( this->private_->get_mutex() );
	
	return this->private_->params_[ StringToLower( key ) ];
}
	
	
void Config::set_param( const std::string& key, const std::string& value )
{
	// Lock internals of the class
	ConfigPrivate::lock_type lock( this->private_->get_mutex() );

	this->private_->params_[ StringToLower( key ) ] = value;
}

} // namespace Core
