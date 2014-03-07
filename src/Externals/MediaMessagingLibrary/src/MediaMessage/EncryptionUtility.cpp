#include <iostream>

#include "EncryptionUtility.h"

namespace BioMesh3d
{
	EncryptionUtility* EncryptionUtility::s_encryption_utility_handle_ = new EncryptionUtility();

	EncryptionUtility::EncryptionUtility() : m_private_rsa_( NULL ), m_public_rsa_( NULL )
	{
		generate_rsa();
	}

	EncryptionUtility::~EncryptionUtility()
	{
		RSA_free( this->m_private_rsa_ );
		RSA_free( this->m_public_rsa_ );
	}

	EncryptionUtility* BioMesh3d::EncryptionUtility::instance()
	{
		return s_encryption_utility_handle_;
	}

	void EncryptionUtility::generate_rsa()
	{
		boost::mutex::scoped_lock lock( this->m_private_rsa_mutex_ );
		RAND_status();
		this->m_private_rsa_ = RSA_generate_key( 1024, RSA_F4, NULL, NULL );
		if ( RSA_check_key( this->m_private_rsa_ ) != 1 )
		{
			std::cout << "There was an error while generating the keys";
		}
	}

	std::vector< char > EncryptionUtility::encrypt_string_rsa( const std::string& public_key_string,
		const std::string& plain_text )
	{
		set_public_key_rsa( public_key_string );
		boost::mutex::scoped_lock lock( this->m_public_rsa_mutex_ );
		std::vector< char > cypher_text;
		cypher_text.resize( 10000 );
		int success = RSA_public_encrypt( plain_text.size(), ( unsigned char* )plain_text.c_str(), 
			( unsigned char* )&cypher_text[0], this->m_public_rsa_, RSA_PKCS1_PADDING );
		if ( success == -1 )
		{
			std::cout << "There was a problem while encrypting the string" << std::endl;
		}
		cypher_text.resize( success );
		return cypher_text;
	}

	std::string EncryptionUtility::decrypt_string_rsa( const std::vector< char >& cypher_text )
	{
		boost::mutex::scoped_lock lock( this->m_private_rsa_mutex_ );
		
		std::vector< unsigned char > translated_text;
		translated_text.resize( 10000 );
		int success = RSA_private_decrypt( cypher_text.size(), ( unsigned char* )&cypher_text[0], 
			( unsigned char* )&translated_text[0], this->m_private_rsa_, RSA_PKCS1_PADDING );
		if ( success == -1 )
		{
			std::cout << "There was a problem while decrypting the string" << std::endl;
		}
		std::string decrypted_string( ( char* )&translated_text[0], success );
		return decrypted_string;
	}

	void EncryptionUtility::set_public_key_rsa( const std::string& public_key_string )
	{
		boost::mutex::scoped_lock lock( this->m_public_rsa_mutex_ );
		if (  this->m_public_key_pem_.compare( public_key_string ) == 0 ) { return; } // we don't need to generate a new public key
		BIO* mem_ptr = NULL;
		EVP_PKEY* pkey = NULL;
		RSA_free( this->m_public_rsa_ );

		mem_ptr = BIO_new( BIO_s_mem() );
		BIO_puts( mem_ptr, public_key_string.c_str() );
		pkey = PEM_read_bio_PUBKEY( mem_ptr, NULL, NULL, NULL );
		this->m_public_rsa_ = EVP_PKEY_get1_RSA( pkey );
		BIO_free( mem_ptr );
	}

	std::string EncryptionUtility::get_pem_public_key_rsa()
	{
		boost::mutex::scoped_lock lock( this->m_private_rsa_mutex_ );
		BIO* mem_ptr = NULL;
		mem_ptr = BIO_new( BIO_s_mem() );
//		BIO_set_close( mem_ptr, BIO_NOCLOSE );
		PEM_write_bio_RSA_PUBKEY( mem_ptr, this->m_private_rsa_ );
		char* public_pem;
		int pem_size = BIO_get_mem_data( mem_ptr, &public_pem );
		std::string public_pem_string( public_pem, pem_size );
		BIO_free( mem_ptr );
		return public_pem_string;
	}

	std::vector< unsigned char > EncryptionUtility::SHA512_hash( const std::string& message )
	{
		SHA512_CTX context;
		std::vector< unsigned char > md;
		md.resize( SHA512_DIGEST_LENGTH );

		SHA512_Init( &context );
		SHA512_Update( &context, ( unsigned char* )message.c_str(), message.size() );
		SHA512_Final( ( unsigned char* )&md[0], &context );
		return md;
	}

}
