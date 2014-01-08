/*
 For more information, please see: http://software.sci.utah.edu

 The MIT License

 Copyright (c) 2009 Scientific Computing and Imaging Institute,
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

#ifndef CORE_UTILS_KEYVALUEPAIRLIST_H
#define CORE_UTILS_KEYVALUEPAIRLIST_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

// Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

// STL includes
#include <string>
#include <map>
#include <set>

// Core includes
#include <Core/Utils/StringUtil.h>

namespace Core
{

class KeyValuePairList;
typedef boost::shared_ptr<KeyValuePairList> KeyValuePairListHandle;

class KeyValuePairList
{
	// -- constructors --
public:
	KeyValuePairList();
	KeyValuePairList( const std::string& );

	// -- options --
public:
	// FAIL_ON_NEW_KEYS:
	// After setting defaults do not allow for new keys to be added.
	// NOTE: This is useful to parse command line parameters and warn the users about options
	// that the program does not know about
	void fail_on_new_key( bool enable );

	// -- serialize and unserialize --
public:
	// EXPORT_TO_STRING:
	// Serialize the key value pairs into string
	std::string export_to_string() const;

	// IMPORT_FROM_STRING:
	// Import a series from list of key value pairs
	bool import_from_string( const std::string& key_value_string );

	// IMPORT_FROM_STRING:
	// Import a series from list of key value pairs
	// NOTE: The only difference with the one above is that this one returns an error string
	bool import_from_string( const std::string& key_value_string, std::string& error );

	// IMPORT_FROM_STRING:
	// Import a series from list of key value pairs, also determine which entries were changed
	bool import_from_string( const std::string& key_value_string, std::set<std::string>& changed_keys );

	// IMPORT_FROM_STRING:
	// Import a series from list of key value pairs, also determine which entries were changed
	// NOTE: The only difference with the one above is that this one returns an error string
	bool import_from_string( const std::string& key_value_string, std::set<std::string>& changed_keys, std::string& error );
	
	// -- file IO --
public:
	// EXPORT_TO_FILE:
	// Export the contents of this key/value pair list to disk
	bool export_to_file( const boost::filesystem::path& filename ) const;

	// EXPORT_TO_FILE:
	// Export the contents of this key/value pair list to disk
	// NOTE: The only difference with the one above is that this one returns an error string
	bool export_to_file( const boost::filesystem::path& filename, std::string& error ) const;

	// IMPORT_FROM_FILE:
	// Import the contents of a file into this key/value pair list
	bool import_from_file( const boost::filesystem::path& filename );

	// IMPORT_FROM_FILE:
	// Import the contents of a file into this key/value pair list
	// NOTE: The only difference with the one above is that this one returns an error string
	bool import_from_file( const boost::filesystem::path& filename, std::string& error );


	// -- input arguments --
public:
	// IMPORT_FROM_ARGUMENTS:
	// Import from standard input arguments.
	bool import_from_arguments( int argc, char** argv );
	
	// IMPORT_FROM_ARGUMENTS:
	// Import from standard input arguments.
	// NOTE: The only difference with the one above is that this one returns an error string
	bool import_from_arguments( int argc, char** argv, std::string& error );
	
	// -- get parameters --
public:
	// GET_VALUE:
	// Get a value for a certain key.
	bool get_value( const std::string& key, std::string& value ) const;
	
	// SET_VALUE:
	// Set a value for a certain key.
	bool set_value( const std::string& key, const std::string& value );

	// GET:
	// Templated function that retrieves value from the string in which it is strored.
	template< class T >
	bool get( const std::string& key, T& value ) const
	{
		std::string value_string;
		if (! this->get_value( key, value_string ) ) return false;
		return ImportFromString( value_string, value );
	}

	// GET:
	// Templated function that retrieves value from the string in which it is strored.
	// NOTE: This version return an error if the value could not be converted
	template< class T >
	bool get( const std::string& key, T& value, std::string& error ) const
	{
		error.clear();
		
		std::string value_string;
		if (! this->get_value( key, value_string ) ) return false;
		if (! ImportFromString( value_string, value ) )
		{
			error = std::string( "Could not interpret value '" ) + value_string + "'.";
			return false;
		}
		
		return true;
	}


	bool is_key( const std::string& key ) const
	{
		return ( this->key_value_.find( key ) != this->key_value_.end() );
	}


	// GET_AS_STRING:
	std::string get_as_string( const std::string& key )
	{
		std::string value;
		if (! get( key, value ) ) return "";
		return value;
	}

	// GET_AS_BOOL:
	bool get_as_bool( const std::string& key )
	{
		bool value;
		if (! get( key, value ) ) return false;
		return value;
	}

	// GET_AS_INT:
	int get_as_int( const std::string& key )
	{
		int value;
		if (! get( key, value ) ) return 0;
		return value;
	}

	// GET_AS_DOUBLE:
	double get_as_double( const std::string& key )
	{
		double value;
		if (! get( key, value ) ) return 0.0;
		return value;
	}

	// SET:
	// Set and convert the value to a string.
	template< class T >
	bool set( const std::string & key, const T& value )
	{
		return set_value( key, ExportToString( value ) );
	}

	// -- clear --
public:
	void clear();

	// -- internals --
private:
	std::map<std::string,std::string> key_value_;
	
	bool fail_on_new_key_;
};


} // end namespace Core

#endif
