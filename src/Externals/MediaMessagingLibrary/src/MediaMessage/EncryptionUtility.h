#ifndef ENCRYPTION_UTILITY_H
#define ENCRYPTION_UTILITY_H

#define WIN32_LEAN_AND_MEAN
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <openssl/rand.h>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/thread/condition.hpp>

#include <vector>

namespace BioMesh3d
{
	class EncryptionUtility : public boost::noncopyable
	{
	public:
		static EncryptionUtility* instance();

		std::vector< char > encrypt_string_rsa( const std::string& public_key_string, 
			const std::string& plain_text );
		std::string decrypt_string_rsa( const std::vector< char >& cypher_text );
		std::string get_pem_public_key_rsa();

		std::vector< unsigned char > SHA512_hash( const std::string& message );
	private:
		EncryptionUtility();
		~EncryptionUtility();

		void set_public_key_rsa( const std::string& public_key_string );
		void generate_rsa();

		boost::mutex m_private_rsa_mutex_;
		RSA* m_private_rsa_;

		boost::mutex m_public_rsa_mutex_;
		RSA* m_public_rsa_;
		std::string m_public_key_pem_;

		static EncryptionUtility* s_encryption_utility_handle_;
	};
}

#endif