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

// Boost includes
#include <boost/filesystem.hpp>

// STL includes
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <fstream>
#include <iomanip>

// Core includes
#include <Core/Utils/Log.h>
#include <Core/Utils/StringUtil.h>
#include <Core/Math/MathFunctions.h>

namespace Core
{

// Convert a value into a string

template< class T >
bool FromString( const std::string &str, T &value )
{
	std::string data = str + " ";
	for ( size_t j = 0; j < data.size(); j++ )
		if ( ( data[ j ] == '\t' ) || ( data[ j ] == '\r' ) || ( data[ j ] == '\n' ) || ( data[ j ]
																						 == '"' ) || ( data[ j ] == ',' ) || ( data[ j ] == '[' ) || ( data[ j ] == ']' )
			|| ( data[ j ] == '(' ) || ( data[ j ] == ')' ) ) data[ j ] = ' ';
	
	std::istringstream iss( data );
	iss.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
	try
	{
		iss >> value;
		return true;
	}
	catch ( ... )
	{
		return false;
	}
}

// Convert multiple values in a string into a vector with numbers

template< class T >
bool MultipleFromString( const std::string &str, std::vector< T > &values )
{
	values.clear();

	// Clear out any markup of the numbers that make it easier to read and
	// replace it all with spaces.
	std::string data = str;
	for ( size_t j = 0; j < data.size(); j++ )
		if ( ( data[ j ] == '\t' ) || ( data[ j ] == '\r' ) || ( data[ j ] == '\n' ) || ( data[ j ]
		    == '"' ) || ( data[ j ] == ',' ) || ( data[ j ] == '[' ) || ( data[ j ] == ']' )
		    || ( data[ j ] == '(' ) || ( data[ j ] == ')' ) ) data[ j ] = ' ';

	// Loop over the data and extract all numbers from it.
	for ( size_t p = 0; p < data.size(); )
	{
		// find where the number starts
		while ( ( p < data.size() ) && ( data[ p ] == ' ' ) )
			p++;
		// Exit if we are at the end of the series
		if ( p >= data.size() ) break;

		// strip of the next number
		std::string::size_type next_space = data.find( ' ', p );
		if ( next_space == std::string::npos ) next_space = data.size();

		// Extract the number
		T value;
		if ( FromString( data.substr( p, next_space - p ), value ) ) values.push_back( value );
		p = next_space;

		if ( p >= data.size() ) break;
	}

	// If no numbers were extracted return false
	return true;
}


// Export a value to a string
template< class T >
std::string ToString( T val )
{
	std::ostringstream oss;
	oss << val;
	return oss.str();
}

inline std::string ToString( float val )
{
	if ( IsNan( val ) ) return "NaN";
	if ( IsInfinite( val ) )
	{
		if ( val > 0 ) return "INF";
		return "-INF";
	}
	
	std::ostringstream oss;
	oss << std::showpoint << val;
	return oss.str();
}

inline std::string ToString( double val )
{
	if ( IsNan( val ) ) return "NaN";
	if ( IsInfinite( val ) )
	{
		if ( val > 0 ) return "INF";
		return "-INF";
	}

	std::ostringstream oss;
	oss << std::showpoint << val;
	return oss.str();
}

// Export a value to a string with precision control

template< class T >
std::string ToString( T val, int precision )
{
	std::ostringstream oss;
	oss.precision( precision );
	oss << std::showpoint << val;
	return oss.str();
}


inline std::string ToString( unsigned char val, int precision )
{
	std::ostringstream oss;
	
	oss.precision( precision );
	oss << std::right << std::setfill( '0' ) << std::setw( precision ) << val;
	return oss.str();
}


inline std::string ToString( unsigned short val, int precision )
{
	std::ostringstream oss;
	
	oss.precision( precision );
	oss << std::right << std::setfill( '0' ) << std::setw( precision ) << val;
	return oss.str();
}


inline std::string ToString( unsigned int val, int precision )
{
	std::ostringstream oss;
	
	oss.precision( precision );
	oss << std::right << std::setfill( '0' ) << std::setw( precision ) << val;
	return oss.str();
}


inline std::string ToString( unsigned long long val, int precision )
{
	std::ostringstream oss;
	
	oss.precision( precision );
	oss << std::right << std::setfill( '0' ) << std::setw( precision ) << val;
	return oss.str();
}



inline std::string ToString( float val, int precision )
{
	if ( IsNan( val ) ) return "NaN";
	if ( IsInfinite( val ) )
	{
		if ( val > 0 ) return "INF";
		return "-INF";
	}

	std::ostringstream oss;
	oss.precision( precision );
	oss << std::fixed << val;
	return oss.str();
}

inline std::string ToString( double val, int precision )
{
	if ( IsNan( val ) ) return "NaN";
	if ( IsInfinite( val ) )
	{
		if ( val > 0 ) return "INF";
		return "-INF";
	}

	std::ostringstream oss;
	oss.precision( precision );
	oss << std::fixed << val;
	return oss.str();
}

inline std::string ToString( double val, size_t digits )
{
	if ( IsNan( val ) ) return "NaN";
	if ( IsInfinite( val ) )
	{
		if ( val > 0 ) return "INF";
		return "-INF";
	}
	
	std::ostringstream oss;
	
	oss << std::fixed << std::setprecision( digits ) << val;
	return oss.str();
}


std::string StringToUpper( std::string str )
{
	std::string::iterator iter = str.begin();
	std::string::iterator iend = str.end();
	for ( ; iter != iend; ++iter )
		*iter = toupper( *iter );
	return str;
}

std::string StringToLower( std::string str )
{
	std::string::iterator iter = str.begin();
	std::string::iterator iend = str.end();
	for ( ; iter != iend; ++iter )
		*iter = tolower( *iter );
	return str;
}

bool FromString( const std::string &str, double &value )
{
	// Clear out any markup of the numbers that make it easier to read and
	// replace it all with spaces.
	std::string data = str + " ";
	for ( size_t j = 0; j < data.size(); j++ )
		if ( ( data[ j ] == '\t' ) || ( data[ j ] == '\r' ) || ( data[ j ] == '\n' ) || ( data[ j ]
		    == '"' ) || ( data[ j ] == ',' ) || ( data[ j ] == '[' ) || ( data[ j ] == ']' )
		    || ( data[ j ] == '(' ) || ( data[ j ] == ')' ) ) data[ j ] = ' ';

	// if empty just return
	if ( data.size() == 0 ) return ( false );

	// Handle special cases: nan, inf, and -inf

	// handle nan
	if ( data.size() > 2 && ( data[ 0 ] == 'n' || data[ 0 ] == 'N' ) && ( data[ 1 ] == 'a'
	    || data[ 1 ] == 'A' ) && ( data[ 2 ] == 'n' || data[ 2 ] == 'N' ) )
	{
		value = std::numeric_limits< double >::quiet_NaN();
		return ( true );
	}
	// handle inf
	else if ( data.size() > 2 && ( data[ 0 ] == 'i' || data[ 0 ] == 'I' ) && ( data[ 1 ] == 'n'
	    || data[ 1 ] == 'N' ) && ( data[ 2 ] == 'f' || data[ 2 ] == 'F' ) )
	{
		value = std::numeric_limits< double >::infinity();
		return ( true );
	}
	// handle +inf and -inf
	else if ( data.size() > 3 && ( data[ 0 ] == '-' || data[ 0 ] == '+' ) && ( data[ 1 ] == 'i'
	    || data[ 1 ] == 'I' ) && ( data[ 2 ] == 'n' || data[ 2 ] == 'N' ) && ( data[ 3 ] == 'f'
	    || data[ 3 ] == 'F' ) )
	{
		if ( data[ 0 ] == '-' )
		{
			value = -std::numeric_limits< double >::infinity();
		}
		else
		{
			value = std::numeric_limits< double >::infinity();
		}

		return ( true );
	}

	std::istringstream iss( data );
	iss.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
	try
	{
		iss >> value;
		return ( true );
	}
	catch ( std::istringstream::failure e )
	{
		CORE_LOG_DEBUG(e.what());
		return ( false );
	}
}

bool FromString( const std::string &str, float &value )
{
	// Clear out any markup of the numbers that make it easier to read and
	// replace it all with spaces.
	std::string data = str + " ";
	for ( size_t j = 0; j < data.size(); j++ )
		if ( ( data[ j ] == '\t' ) || ( data[ j ] == '\r' ) || ( data[ j ] == '\n' ) || ( data[ j ]
		    == '"' ) || ( data[ j ] == ',' ) || ( data[ j ] == '[' ) || ( data[ j ] == ']' )
		    || ( data[ j ] == '(' ) || ( data[ j ] == ')' ) ) data[ j ] = ' ';

	// if empty just return
	if ( data.size() == 0 ) return ( false );

	// Handle special cases: nan, inf, and -inf
	// Handle special cases: nan, inf, and -inf

	// handle nan
	if ( data.size() > 2 && ( data[ 0 ] == 'n' || data[ 0 ] == 'N' ) && ( data[ 1 ] == 'a'
	    || data[ 1 ] == 'A' ) && ( data[ 2 ] == 'n' || data[ 2 ] == 'N' ) )
	{
		value = std::numeric_limits< float >::quiet_NaN();
		return ( true );
	}
	// handle inf
	else if ( data.size() > 2 && ( data[ 0 ] == 'i' || data[ 0 ] == 'I' ) && ( data[ 1 ] == 'n'
	    || data[ 1 ] == 'N' ) && ( data[ 2 ] == 'f' || data[ 2 ] == 'F' ) )
	{
		value = std::numeric_limits< float >::infinity();
		return ( true );
	}
	// handle +inf and -inf
	else if ( data.size() > 3 && ( data[ 0 ] == '-' || data[ 0 ] == '+' ) && ( data[ 1 ] == 'i'
	    || data[ 1 ] == 'I' ) && ( data[ 2 ] == 'n' || data[ 2 ] == 'N' ) && ( data[ 3 ] == 'f'
	    || data[ 3 ] == 'F' ) )
	{
		if ( data[ 0 ] == '-' )
		{
			value = -std::numeric_limits< float >::infinity();
		}
		else
		{
			value = std::numeric_limits< float >::infinity();
		}

		return ( true );
	}

	std::istringstream iss( data );
	iss.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
	try
	{
		iss >> value;
		return ( true );
	}
	catch ( ... )
	{
		return ( false );
	}
}

// Strip out space at the start and at the end of the string
void StripSpaces( std::string& str )
{
	size_t esize = str.size();
	size_t idx = 0;

	// Strip out spaces at the start of the str
	while ( ( idx < esize ) && ( ( str[ idx ] == ' ' ) || ( str[ idx ] == '\t' ) || ( str[ idx ]
	    == '\n' ) || ( str[ idx ] == '\r' ) || ( str[ idx ] == '\f' ) || ( str[ idx ] == '\v' ) ) )
		idx++;

	// Get the substring without spaces at the start 
	str = str.substr( idx, ( str.size() - idx ) );
}

// Strip out space at the start and at the end of the string
void StripSurroundingSpaces( std::string& str )
{
	size_t esize = str.size();
	size_t idx = 0;

	// Strip out spaces at the start of the str
	while ( ( idx < esize ) && ( ( str[ idx ] == ' ' ) || ( str[ idx ] == '\t' ) || ( str[ idx ]
	    == '\n' ) || ( str[ idx ] == '\r' ) || ( str[ idx ] == '\f' ) || ( str[ idx ] == '\v' ) ) )
		idx++;

	size_t ridx = 0;
	if ( str.size() ) ridx = str.size() - 1;

	// Strip out spaces at the end of the str
	while ( ( ridx > 0 ) && ( ( str[ ridx ] == ' ' ) || ( str[ ridx ] == '\t' ) || ( str[ ridx ]
	    == '\n' ) || ( str[ ridx ] == '\r' ) || ( str[ ridx ] == '\f' ) || ( str[ ridx ] == '\v' ) ) )
		ridx--;

	// Get the substring without spaces at the start or at the end
	str = str.substr( idx, ( ridx - idx + 1 ) );
}

// Function to split a list of options delimited by a character into a vector of
// strings
std::vector<std::string> SplitString( const std::string& str, const std::string& delimiter )
{
	std::string option_list_string = str;
	std::vector<std::string> option_list;
	while ( 1 )
	{
		size_t loc = option_list_string.find( delimiter );
		if ( loc >= option_list_string.size() )
		{
			option_list.push_back( option_list_string );
			break;
		}
		option_list.push_back( option_list_string.substr( 0, loc ) );
		option_list_string = option_list_string.substr( loc + delimiter.size() );
	}

	return option_list;
}


std::string ExportToString( bool value )
{
	if ( value ) return ( std::string( "True" ) );
	else return ( std::string( "False" ) );
}

std::string ExportToString( char value )
{
	return ToString( value );
}

std::string ExportToString( unsigned char value )
{
	return ToString(value);
}

std::string ExportToString( short value )
{
	return ToString( value );
}

std::string ExportToString( unsigned short value )
{
	return ToString(value);
}

std::string ExportToString( int value )
{
	return ToString( value );
}

std::string ExportToString( unsigned int value )
{
	return ToString(value);
}

std::string ExportToString( long value )
{
	return ToString( value );
}

std::string ExportToString( unsigned long value )
{
	return ToString( value );
}

std::string ExportToString( long long value )
{
	return ToString( value );
}

std::string ExportToString( unsigned long long value )
{
	return ToString( value );
}

std::string ExportToString( float value )
{
	return ToString( value );
}

std::string ExportToString( double value )
{
	return ToString( value );
}

std::string ExportToString( unsigned char value, int precision )
{
	return ToString( value, precision );
}

std::string ExportToString( unsigned short value, int precision )
{
	return ToString( value, precision );
}

std::string ExportToString( unsigned int value, int precision )
{
	return ToString( value, precision );
}

std::string ExportToString( unsigned long long value, int precision )
{
	return ToString( value, precision );
}

std::string ExportToString( float value, int precision )
{
	return ToString( value, precision );
}

std::string ExportToString( double value, int precision )
{
	return ToString( value, precision );
}

std::string ExportToString( const double& value, size_t digits )
{
	return ToString( value, digits );
}

std::string ExportToString( const std::string& value )
{
	bool need_quotes = true;
//	for ( size_t j = 0; j < value.size(); j++)
//	{
//		if ( value[j] == ' ' || value[j] == '\t' || value[j] == '[' || value[j] == ']' ||
//			value[j] == '(' || value[j] == ')' || value[j] == ',' || value[j] == '"' ) 
//		{
//			need_quotes = true;
//		}
//	}
//	if ( value.size() == 0) need_quotes = true;
	
	if ( need_quotes ) return std::string(1,'"') + value + std::string(1,'"');
	else return value;
}

std::string ExportToString( const std::vector< char >& value )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
			result += ToString( value[ j ] ) + ", ";
		 result += ToString( value[ value.size() - 1 ] );
	}
	result += "]";
	return result;
}
	
std::string ExportToString( const std::vector< std::string >& value )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ExportToString( value[ j ] ) + ", ";
		}
		result += ExportToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString(const std::vector< unsigned char >& value)
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ] ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< short >& value )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ] ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< unsigned short >& value )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ] ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< int >& value )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ] ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< unsigned int >& value)
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ] ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< long >& value )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ] ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< unsigned long >& value )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ] ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< long long >& value )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ] ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< unsigned long long >& value )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ] ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< float >& value )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ] ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< double >& value )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ] ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ] );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< float >& value, int precision )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ], precision  ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ], precision  );	
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::vector< double >& value, int precision )
{
	std::string result( "[ " );
	if ( value.size() )
	{
		for ( size_t j = 0; j < value.size() - 1; j++ )
		{
			result += ToString( value[ j ], precision  ) + ", ";
		}
		result += ToString( value[ value.size() - 1 ], precision  );	
	}
	result += " ]";
	return result;
}


std::string ExportToString( const std::set< std::string >& value )
{
	std::string result( "[ " );
	std::set< std::string >::const_iterator it = value.begin();
	
	if ( it != value.end() )
	{
		result += ExportToString( *it );
		++it;
	}

	while ( it != value.end() )
	{
		result += ", " + ExportToString( *it );
		++it;
	}
	result += "]";
	return result;
}


std::string ExportToString( const std::set< int >& value )
{
	std::string result( "[ " );
	std::set< int >::const_iterator it = value.begin();

	if ( it != value.end() )
	{
		result += ToString( *it );
		++it;
	}

	while ( it != value.end() )
	{
		result += ", " + ToString( *it );
		++it;
	}
	result += "]";
	return result;
}


std::string ExportToString( const std::set< unsigned int >& value )
{
	std::string result( "[ " );
	std::set< unsigned int >::const_iterator it = value.begin();

	if ( it != value.end() )
	{
		result += ToString( *it );
		++it;
	}

	while ( it != value.end() )
	{
		result += ", " + ToString( *it );
		++it;
	}
	result += "]";
	return result;
}

std::string ExportToString( const std::set< long long >& value )
{
	std::string result( "[ " );
	std::set< long long >::const_iterator it = value.begin();

	if ( it != value.end() )
	{
		result += ToString( *it );
		++it;
	}

	while ( it != value.end() )
	{
		result += ", " + ToString( *it );
		++it;
	}
	result += "]";
	return result;
}


std::string ExportToString( const std::set< unsigned long long >& value )
{
	std::string result( "[ " );
	std::set< unsigned long long >::const_iterator it = value.begin();

	if ( it != value.end() )
	{
		result += ToString( *it );
		++it;
	}

	while ( it != value.end() )
	{
		result += ", " + ToString( *it );
		++it;
	}
	result += "]";
	return result;
}


bool ImportFromString( const std::string& str, bool& value )
{
	std::string tmpstr( str );
	StripSurroundingSpaces( tmpstr );
	tmpstr = StringToLower( tmpstr );
	
	if ( ( tmpstr == "0" ) || ( tmpstr == "false" ) || ( tmpstr == "off" ) || ( tmpstr == "no" ) )
	{
		value = false;
		return ( true );
	}
	else if ( ( tmpstr == "1" ) || ( tmpstr == "true" ) || ( tmpstr == "on" ) || ( tmpstr == "yes" ) )
	{
		value = true;
		return ( true );
	}
	return ( false );
}

bool ImportFromString( const std::string& str, char& value )
{
	return ( FromString( str, value ) );
}

bool ImportFromString( const std::string& str, unsigned char& value )
{
	return ( FromString( str, value ) );
}

bool ImportFromString( const std::string& str, short& value )
{
	return ( FromString( str, value ) );
}

bool ImportFromString(const std::string& str, unsigned short& value)
{
	return ( FromString( str, value ) );
}

bool ImportFromString( const std::string& str, int& value )
{
	return ( FromString( str, value ) );
}

bool ImportFromString( const std::string& str, unsigned int& value )
{
	return ( FromString( str, value ) );
}

bool ImportFromString( const std::string& str, long& value )
{
	return ( FromString( str, value ) );
}

bool ImportFromString( const std::string& str, unsigned long& value )
{
	return ( FromString( str, value ) );
}

bool ImportFromString( const std::string& str, long long& value )
{
	return ( FromString( str, value ) );
}

bool ImportFromString( const std::string& str, unsigned long long& value )
{
	return ( FromString( str, value ) );
}

bool ImportFromString( const std::string& str, float& value )
{
	return ( FromString( str, value ) );
}

bool ImportFromString( const std::string& str, double& value )
{
	return ( FromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector< char >& value )
{
	return ( MultipleFromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector< unsigned char >& value )
{
	return ( MultipleFromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector< short >& value )
{
	return ( MultipleFromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector< unsigned short >& value )
{
	return ( MultipleFromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector< int >& value )
{
	return ( MultipleFromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector<unsigned int>& value )
{
	return ( MultipleFromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector<long>& value )
{
	return ( MultipleFromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector<unsigned long>& value )
{
	return ( MultipleFromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector<long long>& value )
{
	return ( MultipleFromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector<unsigned long long>& value )
{
	return ( MultipleFromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector< float >& value )
{
	return ( MultipleFromString( str, value ) );
}

bool ImportFromString( const std::string& str, std::vector< double >& value )
{
	return ( MultipleFromString( str, value ) );
}
	
bool ImportFromString( const std::string& str, std::vector< std::string >& value )
{
	std::string data = str;
	
	StripSurroundingSpaces( data );

	if ( data.size() > 1 )
	{
		if ( ( data[ 0 ] == '[' ) && ( data[ data.size() - 1 ] == ']' ) )
		{
			data = data.substr( 1, data.size() - 2 );
		}
	}
	
	value.clear();
	
	size_t j = 0;
	while ( j < data.size() )
	{
		while ( j < data.size() && ( ( data[ j ] == ' ') || ( data[ j ] == '\t' ) || 
			( data[ j ] == '\r' ) || ( data[ j ] == '\n' ) || ( data[ j ] == ',') ) ) j++;

		if ( j == data.size() ) return true;
		if ( data[ j ] == '[' )
		{
			j++;
			size_t start = j;
			size_t paren_count = 0;

			while ( j < data.size() && ( data[ j ] != ']' || paren_count > 0 ) )
			{
				if ( data[ j ] == '[' ) paren_count++;
				if ( data[ j ] == ']' ) paren_count--;
				j++;
			}

			// if there is no end quotation mark
			if ( j == data.size() ) return false;
			value.push_back( data.substr( start, j-start ) );
			j++;
		}
		else if ( data[ j ] == '"' )
		{
			j++;
			size_t start = j;

			while ( j < data.size() && data[ j ] != '"' )
			{
				j++;
			}

			// if there is no end quotation mark
			if ( j == data.size() ) return false;
			value.push_back( data.substr( start, j-start ) );
			j++;
		}		
		else
		{
			size_t start = j;
			while ( j < data.size() && ( ( data[ j ] != ' ' )  && ( data[ j ] != '\t' ) &&
				( data[ j ] != '\r' ) && ( data[ j ] != '\n' )  && ( data[ j ] != ',' ) ) ) j++;
		
			value.push_back( data.substr( start, j-start ) );
		}
	}

	return true;
}
	

bool ImportFromString( const std::string& str, std::set< std::string >& value )
{
	std::string data = str;

	StripSurroundingSpaces( data );

	if ( data.size() > 1 )
	{
		if ( ( data[ 0 ] == '[' ) && ( data[ data.size() - 1 ] == ']' ) )
		{
			data = data.substr( 1, data.size() - 2 );
		}
	}
	
	value.clear();
	
	size_t j = 0;
	while ( j < data.size() )
	{
		while ( j < data.size() && ( ( data[ j ] == ' ') || ( data[ j ] == '\t' ) || 
			( data[ j ] == '\r' ) || ( data[ j ] == '\n' )|| ( data[ j ] == ',') ) ) j++;

		if ( j == data.size() ) return true;
		if ( data[ j ] == '[' )
		{
			j++;
			size_t start = j;
			size_t paren_count = 0;

			while ( j < data.size() && ( data[ j ] != ']' || paren_count > 0 ) )
			{
				if ( data[ j ] == '[' ) paren_count++;
				if ( data[ j ] == ']' ) paren_count--;
				j++;
			}

			// if there is no end quotation mark
			if ( j == data.size() ) return false;
			value.insert( data.substr( start, j-start ) );
			j++;
		}
		else if ( data[ j ] == '"' )
		{
			j++;
			size_t start = j;

			while ( j < data.size() && data[ j ] != '"' )
			{
				j++;
			}

			// if there is no end quotation mark
			if ( j == data.size() ) return false;
			value.insert( data.substr( start, j-start ) );
			j++;
		}
		else
		{
			size_t start = j;
			while ( j < data.size() && ( ( data[ j ] != ' ' )  && ( data[ j ] != '\t' ) &&
				( data[ j ] != '\r' ) && ( data[ j ] != '\n' ) && ( data[ j ] != ',' ) ) ) j++;
		
			value.insert( data.substr( start, j-start ) );
		}
	}

	return true;
}
	

bool ImportFromString( const std::string& str, std::string& value )
{
	value = str;
	StripSurroundingSpaces( value );
	
	// Remove quotes if needed
	if ( value.size() >= 2 )
	{
		if ( ( value[0] == '"' ) && ( value[ value.size() - 1 ] == '"' ) )
		{
			value = value.substr( 1, value.size() - 2 );
		}
	}
	
	return true;
}


bool ImportFromString( const std::string& str, std::set< int >& value )
{
	value.clear();
	std::vector< int > tmp;
	if ( MultipleFromString( str, tmp ) )
	{
		value.insert( tmp.begin(), tmp.end() );
		return true;
	}
	return false;
}

bool ImportFromString( const std::string& str, std::set< unsigned int >& value )
{
	value.clear();
	std::vector< unsigned int > tmp;
	if ( MultipleFromString( str, tmp ) )
	{
		value.insert( tmp.begin(), tmp.end() );
		return true;
	}
	return false;
}

bool ImportFromString( const std::string& str, std::set< long long >& value )
{
	value.clear();
	std::vector< long long > tmp;
	if ( MultipleFromString( str, tmp ) )
	{
		value.insert( tmp.begin(), tmp.end() );
		return true;
	}
	return false;
}

bool ImportFromString( const std::string& str, std::set< unsigned long long >& value )
{
	value.clear();
	std::vector< unsigned long long > tmp;
	if ( MultipleFromString( str, tmp ) )
	{
		value.insert( tmp.begin(), tmp.end() );
		return true;
	}
	return false;
}

std::string EscapeString( const std::string& str )
{
	std::string escaped_str;
	escaped_str.reserve( str.size() );
	
	for ( size_t j = 0; j < str.size(); j++ )
	{
		char c = str[ j ];
		switch ( c )
		{
			case '\n':
				escaped_str += "\\n";
				break;
			case '\t':
				escaped_str += "\\t";
				break;
			case '\b':
				escaped_str += "\\b";
				break;
			case '\v':
				escaped_str += "\\v";
				break;
			case '\r':
				escaped_str += "\\r";
				break;
			case '\0':
				escaped_str += "\\0";
				break;
			case '\\':
				escaped_str += "\\\\";
				break;
			case '\"':
				escaped_str += "\\\"";
				break;
			default:
				escaped_str += c;		
		}
	}

	return escaped_str;
}

static char* ascii_table = "................................"
					" !\"#$%&'()*+,-./0123456789:;<=>?"
					"@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_"
					"`abcdefghijklmnopqrstuvwxyz{|}~."
					"................................"
					"................................"
					"................................"
					"..................................";

std::string BinaryToString( const unsigned char* data, size_t size )
{
	std::string binary_string;
	
	binary_string.reserve( size );
	for( size_t j = 0; j < size; j ++ )
	{
		binary_string += ascii_table[ data[ j ] ];
	}
	return binary_string;
}


bool ImportFromTextFile( const std::string& filename, std::string& contents )
{
	contents.clear();
	try
	{
		std::ifstream infile( filename.c_str() );
		if ( !infile )
		{
			return false;
		}
		infile.exceptions( std::ifstream::badbit );
		std::vector<char> buffer( 1024 );
		while ( !infile.eof() )
		{
			infile.read( &buffer[ 0 ], buffer.size() );
			contents += std::string( buffer.begin(), buffer.begin() + infile.gcount() );
		}
	}
	catch ( ... )
	{
		return false;
	}
	
	return true;
}

bool ExportToTextFile( const std::string& filename, const std::string& contents )
{
	try
	{
		std::ofstream outfile( filename.c_str() );
		outfile << contents;
	}
	catch ( ... )
	{
		return false;
	}
	
	return true;	
}

bool ExportToTextFileSecurely( const std::string& filename, const std::string& contents )
{
	boost::system::error_code ec;
	boost::filesystem::path file_path( filename );
	
	boost::filesystem::path file_path_new = file_path.parent_path() / 
		( file_path.filename().string() + ".new" );
	
	if (! ExportToTextFile( file_path_new.string(), contents ) )
	{
		return false;
	}
	
	// Remove the old file
	boost::filesystem::remove( file_path, ec );
	if ( ec )
	{
		return false;
	}

	// Replace the old file with the new file
	boost::filesystem::rename( file_path_new, file_path, ec );
	if ( ec )
	{
		return false;
	}
	
	return true;
}


bool ImportFromTextFileSecurely( const std::string& filename, std::string& contents )
{
	boost::system::error_code ec;
	boost::filesystem::path file_path( filename );
	
	// If file does not exist, the .new copy may, if so rename that one and put it
	// on top of the older version
	if (! boost::filesystem::exists( file_path, ec ) )
	{
		boost::filesystem::path file_path_new = file_path.parent_path() / 
			( file_path.filename().string() + ".new" );
		
		if ( boost::filesystem::exists( file_path_new, ec ) )
		{
			boost::filesystem::rename( file_path_new, file_path, ec );
			if ( ec )
			{
				return false;
			}
		}
		else
		{	
			// Neither old file or new file does exist
			return false;
		}
	}

	return ImportFromTextFile( filename, contents );
}

} // End namespace Core

