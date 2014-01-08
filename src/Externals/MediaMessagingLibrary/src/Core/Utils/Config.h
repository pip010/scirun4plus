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

#ifndef CORE_UTILS_CONFIG_H
#define CORE_UTILS_CONFIG_H

// Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

namespace Core
{

class Config;
class ConfigPrivate;
typedef boost::shared_ptr<Config> ConfigHandle;
typedef boost::shared_ptr<ConfigPrivate> ConfigPrivateHandle;

class Config : public boost::noncopyable
{

public:
	Config();
	virtual ~Config();

public:
	// SET_DEFAULTS:
	// Load the configuration file with a series of defaults.
	// The default string is key/value pairs that are separated
	// by the '|' character 
	bool set_defaults( const std::string& defaults );
	
	// LOAD_CONFIG_FILE
	// Load values from the configuration file into the config class
	bool load_config_file( const std::string& config_file );

	// Access the configuration values
public:
	// GET_PARAM:
	// Get a configuration option 
	std::string get_param( const std::string& key ) const;
	
	// SET_PARAM:
	// Set the configuration option
	void set_param( const std::string& key, const std::string& value );
	
	// -- internals --
private:
	ConfigPrivateHandle private_;
	
};

} // namespace Core

#endif
