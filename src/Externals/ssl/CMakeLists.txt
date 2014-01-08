CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(VERSION 2.6)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(SSL)

# set defaults
SET(SCI_CRYPTO_LIBRARY "" CACHE INTERNAL "crypto library name." FORCE)
SET(SCI_CRYPTO_INCLUDE "" CACHE INTERNAL "crypto include directories." FORCE)
SET(SCI_SSL_LIBRARY "" CACHE INTERNAL "crypto library name." FORCE)
SET(SCI_SSL_INCLUDE "" CACHE INTERNAL "crypto include directories." FORCE)
SET(SCI_SSL_FOUND FALSE CACHE INTERNAL "crypto library found" FORCE)

#FIND_PACKAGE(OPENSSL)

IF( OPENSSL_FOUND )

	SET(SCI_CRYPTO_LIBRARY ${OPENSSL_LIBRARIES} CACHE INTERNAL "crypto library name." FORCE)
	SET(SCI_CRYPTO_INCLUDE ${OPENSSL_INCLUDE_DIR} CACHE INTERNAL "crypto include directories." FORCE)
	SET(SCI_SSL_LIBRARY ${OPENSSL_LIBRARIES} CACHE INTERNAL "crypto library name." FORCE)
	SET(SCI_SSL_INCLUDE ${OPENSSL_INCLUDE_DIR} CACHE INTERNAL "crypto include directories." FORCE)
	SET(SCI_SSL_FOUND TRUE CACHE INTERNAL "crypto library found" FORCE)

ELSE( OPENSSL_FOUND )
	SET(OPENSSL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/openssl-1.0.0d)
	OPTION(DOWNLOAD_OPENSSL_SOURCE "Download the openssl source code." OFF)
	IF(NOT DOWNLOAD_OPENSSL_SOURCE)
		MESSAGE(WARNING "OpenSSL package not found. Set DOWNLOAD_OPENSSL_SOURCE to ON to download.")
	ENDIF(NOT DOWNLOAD_OPENSSL_SOURCE)
	
	# Download the source tarball if it doesn't exist
	IF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/openssl-1.0.0d.tar.gz AND DOWNLOAD_OPENSSL_SOURCE)
		FILE(DOWNLOAD http://www.openssl.org/source/openssl-1.0.0d.tar.gz 
					${CMAKE_CURRENT_SOURCE_DIR}/openssl-1.0.0d.tar.gz
					STATUS DOWNLOAD_OPENSSL_STATUS
					EXPECTED_MD5 40b6ea380cc8a5bf9734c2f8bf7e701e
					SHOW_PROGRESS
		)
		LIST(GET ${DOWNLOAD_OPENSSL_STATUS} 0 DOWNLOAD_ERROR_CODE)
		IF(${DOWNLOAD_ERROR_CODE})
			MESSAGE(FATAL_ERROR "Failed to download OpenSSL source code.")
		ENDIF(${DOWNLOAD_ERROR_CODE})
		
		# unpack the tarball
		EXECUTE_PROCESS(
			COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/openssl-1.0.0d.tar.gz
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		)
	ENDIF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/openssl-1.0.0d.tar.gz AND DOWNLOAD_OPENSSL_SOURCE)

	IF( EXISTS ${OPENSSL_SOURCE_DIR}/crypto )

		INCLUDE (CMakeBackwardCompatibilityC)
		SET(CRYPTO_PATH ${OPENSSL_SOURCE_DIR}/crypto )
		SET(SSL_PATH ${OPENSSL_SOURCE_DIR}/ssl )

# IN GENERAL SYMBOLIC LINKS DO NOT ZIP VERY WELL
		
		FILE( GLOB_RECURSE HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h )

		FOREACH( HEADER_FILE ${HEADER_FILES} )
			
			GET_FILENAME_COMPONENT( HEADER_PATH ${HEADER_FILE} PATH )
			GET_FILENAME_COMPONENT( HEADER_PATH_NAME ${HEADER_PATH} NAME )

			IF (NOT ${HEADER_PATH_NAME} STREQUAL "include" )
			IF (NOT ${HEADER_PATH_NAME} STREQUAL "MacOS" )
			IF (NOT ${HEADER_PATH_NAME} STREQUAL "ms" )
				GET_FILENAME_COMPONENT( HEADER_NAME ${HEADER_FILE} NAME )
				CONFIGURE_FILE(${HEADER_FILE}
					${OPENSSL_SOURCE_DIR}/include/openssl/${HEADER_NAME}  COPYONLY )
			ENDIF (NOT ${HEADER_PATH_NAME} STREQUAL "ms" )
			ENDIF (NOT ${HEADER_PATH_NAME} STREQUAL "MacOS" )
			ENDIF (NOT ${HEADER_PATH_NAME} STREQUAL "include" )


		ENDFOREACH( HEADER_FILE ${HEADER_FILES} )

		CONFIGURE_FILE(${OPENSSL_SOURCE_DIR}/crypto/opensslconf.h
			${OPENSSL_SOURCE_DIR}/include/openssl/opensslconf.h  COPYONLY )

		CONFIGURE_FILE(${OPENSSL_SOURCE_DIR}/crypto/opensslv.h
			${OPENSSL_SOURCE_DIR}/include/openssl/opensslv.h  COPYONLY )


		SET(CRYPTO_LIB_SRCS
		  ${CRYPTO_PATH}/cryptlib.c
		  ${CRYPTO_PATH}/mem.c
		  ${CRYPTO_PATH}/mem_clr.c
		  ${CRYPTO_PATH}/mem_dbg.c
		  ${CRYPTO_PATH}/cversion.c
		  ${CRYPTO_PATH}/ex_data.c
		  ${CRYPTO_PATH}/cpt_err.c
		  ${CRYPTO_PATH}/ebcdic.c
		  ${CRYPTO_PATH}/uid.c
		  ${CRYPTO_PATH}/o_time.c
		  ${CRYPTO_PATH}/o_str.c
		  ${CRYPTO_PATH}/o_dir.c
		)

		SET(CRYPTO_LIB_AES_SRCS
		  ${CRYPTO_PATH}/aes/aes_core.c
		  ${CRYPTO_PATH}/aes/aes_misc.c
		  ${CRYPTO_PATH}/aes/aes_ecb.c
		  ${CRYPTO_PATH}/aes/aes_cbc.c
		  ${CRYPTO_PATH}/aes/aes_cfb.c
		  ${CRYPTO_PATH}/aes/aes_ofb.c
		  ${CRYPTO_PATH}/aes/aes_ctr.c
		  ${CRYPTO_PATH}/aes/aes_ige.c
		  ${CRYPTO_PATH}/aes/aes_wrap.c
		)

		SET(CRYPTO_LIB_ASN1_SRCS
		  ${CRYPTO_PATH}/asn1/a_object.c
		  ${CRYPTO_PATH}/asn1/a_bitstr.c 
		  ${CRYPTO_PATH}/asn1/a_utctm.c 
		  ${CRYPTO_PATH}/asn1/a_gentm.c 
		  ${CRYPTO_PATH}/asn1/a_time.c 
		  ${CRYPTO_PATH}/asn1/a_int.c 
		  ${CRYPTO_PATH}/asn1/a_octet.c
		  ${CRYPTO_PATH}/asn1/a_print.c
		  ${CRYPTO_PATH}/asn1/a_type.c
		  ${CRYPTO_PATH}/asn1/a_set.c
		  ${CRYPTO_PATH}/asn1/a_dup.c
		  ${CRYPTO_PATH}/asn1/a_d2i_fp.c
		  ${CRYPTO_PATH}/asn1/a_i2d_fp.c
		  ${CRYPTO_PATH}/asn1/a_enum.c
		  ${CRYPTO_PATH}/asn1/a_utf8.c
		  ${CRYPTO_PATH}/asn1/a_sign.c
		  ${CRYPTO_PATH}/asn1/a_digest.c
		  ${CRYPTO_PATH}/asn1/a_verify.c
		  ${CRYPTO_PATH}/asn1/a_mbstr.c
		  ${CRYPTO_PATH}/asn1/a_strex.c
		  ${CRYPTO_PATH}/asn1/x_algor.c
		  ${CRYPTO_PATH}/asn1/x_val.c
		  ${CRYPTO_PATH}/asn1/x_pubkey.c
		  ${CRYPTO_PATH}/asn1/x_sig.c
		  ${CRYPTO_PATH}/asn1/x_req.c
		  ${CRYPTO_PATH}/asn1/x_attrib.c
		  ${CRYPTO_PATH}/asn1/x_bignum.c
		  ${CRYPTO_PATH}/asn1/x_long.c
		  ${CRYPTO_PATH}/asn1/x_name.c
		  ${CRYPTO_PATH}/asn1/x_x509.c
		  ${CRYPTO_PATH}/asn1/x_x509a.c
		  ${CRYPTO_PATH}/asn1/x_crl.c
		  ${CRYPTO_PATH}/asn1/x_info.c
		  ${CRYPTO_PATH}/asn1/x_spki.c
		  ${CRYPTO_PATH}/asn1/nsseq.c
		  ${CRYPTO_PATH}/asn1/x_nx509.c
		  ${CRYPTO_PATH}/asn1/d2i_pu.c
		  ${CRYPTO_PATH}/asn1/d2i_pr.c
		  ${CRYPTO_PATH}/asn1/i2d_pu.c
		  ${CRYPTO_PATH}/asn1/i2d_pr.c
		  ${CRYPTO_PATH}/asn1/t_req.c
		  ${CRYPTO_PATH}/asn1/t_x509.c
		  ${CRYPTO_PATH}/asn1/t_x509a.c
		  ${CRYPTO_PATH}/asn1/t_crl.c
		  ${CRYPTO_PATH}/asn1/t_pkey.c
		  ${CRYPTO_PATH}/asn1/t_spki.c
		  ${CRYPTO_PATH}/asn1/t_bitst.c
		  ${CRYPTO_PATH}/asn1/tasn_new.c
		  ${CRYPTO_PATH}/asn1/tasn_fre.c
		  ${CRYPTO_PATH}/asn1/tasn_enc.c
		  ${CRYPTO_PATH}/asn1/tasn_dec.c
		  ${CRYPTO_PATH}/asn1/tasn_utl.c
		  ${CRYPTO_PATH}/asn1/tasn_typ.c
		  ${CRYPTO_PATH}/asn1/tasn_prn.c
		  ${CRYPTO_PATH}/asn1/ameth_lib.c
		  ${CRYPTO_PATH}/asn1/f_int.c
		  ${CRYPTO_PATH}/asn1/f_string.c
		  ${CRYPTO_PATH}/asn1/n_pkey.c
		  ${CRYPTO_PATH}/asn1/f_enum.c
		  ${CRYPTO_PATH}/asn1/x_pkey.c
		  ${CRYPTO_PATH}/asn1/a_bool.c
		  ${CRYPTO_PATH}/asn1/x_exten.c
		  ${CRYPTO_PATH}/asn1/bio_asn1.c
		  ${CRYPTO_PATH}/asn1/bio_ndef.c
		  ${CRYPTO_PATH}/asn1/asn_mime.c
		  ${CRYPTO_PATH}/asn1/asn1_gen.c
		  ${CRYPTO_PATH}/asn1/asn1_par.c
		  ${CRYPTO_PATH}/asn1/asn1_lib.c
		  ${CRYPTO_PATH}/asn1/asn1_err.c
		  ${CRYPTO_PATH}/asn1/a_bytes.c
		  ${CRYPTO_PATH}/asn1/a_strnid.c
		  ${CRYPTO_PATH}/asn1/evp_asn1.c
		  ${CRYPTO_PATH}/asn1/p5_pbe.c
		  ${CRYPTO_PATH}/asn1/p5_pbev2.c
		  ${CRYPTO_PATH}/asn1/p8_pkey.c
		  ${CRYPTO_PATH}/asn1/asn_pack.c
		  ${CRYPTO_PATH}/asn1/asn_moid.c
		  ${CRYPTO_PATH}/asn1/asn_mime.c
		)

		SET(CRYPTO_LIB_BF_SRCS
		  ${CRYPTO_PATH}/bf/bf_skey.c 
		  ${CRYPTO_PATH}/bf/bf_ecb.c 
		  ${CRYPTO_PATH}/bf/bf_enc.c 
		  ${CRYPTO_PATH}/bf/bf_cfb64.c 
		  ${CRYPTO_PATH}/bf/bf_ofb64.c 
		)

		SET(CRYPTO_LIB_BN_SRCS
		  ${CRYPTO_PATH}/bn/bn_add.c 
		  ${CRYPTO_PATH}/bn/bn_div.c
		  ${CRYPTO_PATH}/bn/bn_exp.c
		  ${CRYPTO_PATH}/bn/bn_lib.c
		  ${CRYPTO_PATH}/bn/bn_ctx.c
		  ${CRYPTO_PATH}/bn/bn_mul.c
		  ${CRYPTO_PATH}/bn/bn_mod.c
		  ${CRYPTO_PATH}/bn/bn_print.c
		  ${CRYPTO_PATH}/bn/bn_rand.c
		  ${CRYPTO_PATH}/bn/bn_shift.c  
		  ${CRYPTO_PATH}/bn/bn_word.c
		  ${CRYPTO_PATH}/bn/bn_blind.c
		  ${CRYPTO_PATH}/bn/bn_kron.c
		  ${CRYPTO_PATH}/bn/bn_sqrt.c
		  ${CRYPTO_PATH}/bn/bn_gcd.c
		  ${CRYPTO_PATH}/bn/bn_prime.c
		  ${CRYPTO_PATH}/bn/bn_err.c
		  ${CRYPTO_PATH}/bn/bn_sqr.c
		  ${CRYPTO_PATH}/bn/bn_asm.c
		  ${CRYPTO_PATH}/bn/bn_recp.c
		  ${CRYPTO_PATH}/bn/bn_mont.c
		  ${CRYPTO_PATH}/bn/bn_mpi.c
		  ${CRYPTO_PATH}/bn/bn_exp2.c
		  ${CRYPTO_PATH}/bn/bn_gf2m.c
		  ${CRYPTO_PATH}/bn/bn_nist.c
		  ${CRYPTO_PATH}/bn/bn_depr.c
		  ${CRYPTO_PATH}/bn/bn_const.c
		)

		SET(CRYPTO_LIB_BUFFER_SRCS
		  ${CRYPTO_PATH}/buffer/buffer.c 
		  ${CRYPTO_PATH}/buffer/buf_err.c
		)

		SET(CRYPTO_LIB_CMS_SRCS
		  ${CRYPTO_PATH}/cms/cms_lib.c
		  ${CRYPTO_PATH}/cms/cms_asn1.c
		  ${CRYPTO_PATH}/cms/cms_att.c
		  ${CRYPTO_PATH}/cms/cms_io.c
		  ${CRYPTO_PATH}/cms/cms_smime.c
		  ${CRYPTO_PATH}/cms/cms_err.c
		  ${CRYPTO_PATH}/cms/cms_sd.c
		  ${CRYPTO_PATH}/cms/cms_dd.c
		  ${CRYPTO_PATH}/cms/cms_cd.c
		  ${CRYPTO_PATH}/cms/cms_env.c
		  ${CRYPTO_PATH}/cms/cms_enc.c
		  ${CRYPTO_PATH}/cms/cms_ess.c
		)

		SET(CRYPTO_LIB_CONF_SRCS
		  ${CRYPTO_PATH}/conf/conf_err.c
		  ${CRYPTO_PATH}/conf/conf_lib.c
		  ${CRYPTO_PATH}/conf/conf_api.c
		  ${CRYPTO_PATH}/conf/conf_def.c
		  ${CRYPTO_PATH}/conf/conf_mod.c
		  ${CRYPTO_PATH}/conf/conf_mall.c
		  ${CRYPTO_PATH}/conf/conf_sap.c
		)

		SET(CRYPTO_LIB_COMP_SRCS
		  ${CRYPTO_PATH}/comp/comp_lib.c 
		  ${CRYPTO_PATH}/comp/comp_err.c
		  ${CRYPTO_PATH}/comp/c_rle.c
		  ${CRYPTO_PATH}/comp/c_zlib.c
		)

		SET(CRYPTO_LIB_DES_SRCS
		  ${CRYPTO_PATH}/des/cbc_cksm.c
		  ${CRYPTO_PATH}/des/cbc_enc.c
		  ${CRYPTO_PATH}/des/cfb64enc.c
		  ${CRYPTO_PATH}/des/cfb_enc.c
		  ${CRYPTO_PATH}/des/ecb3_enc.c
		  ${CRYPTO_PATH}/des/ecb_enc.c
		  ${CRYPTO_PATH}/des/enc_read.c
		  ${CRYPTO_PATH}/des/enc_writ.c 
		  ${CRYPTO_PATH}/des/fcrypt.c
		  ${CRYPTO_PATH}/des/ofb64enc.c
		  ${CRYPTO_PATH}/des/ofb_enc.c
		  ${CRYPTO_PATH}/des/pcbc_enc.c
		  ${CRYPTO_PATH}/des/qud_cksm.c
		  ${CRYPTO_PATH}/des/rand_key.c
		  ${CRYPTO_PATH}/des/rpc_enc.c
		  ${CRYPTO_PATH}/des/set_key.c
		  ${CRYPTO_PATH}/des/des_enc.c
		  ${CRYPTO_PATH}/des/fcrypt_b.c
		  ${CRYPTO_PATH}/des/xcbc_enc.c
		  ${CRYPTO_PATH}/des/str2key.c
		  ${CRYPTO_PATH}/des/cfb64ede.c
		  ${CRYPTO_PATH}/des/ofb64ede.c
		  ${CRYPTO_PATH}/des/ede_cbcm_enc.c
		  ${CRYPTO_PATH}/des/des_old.c
		  ${CRYPTO_PATH}/des/des_old2.c
		  ${CRYPTO_PATH}/des/read2pwd.c
		)

		SET(CRYPTO_LIB_DH_SRCS
		  ${CRYPTO_PATH}/dh/dh_asn1.c
		  ${CRYPTO_PATH}/dh/dh_gen.c
		  ${CRYPTO_PATH}/dh/dh_key.c
		  ${CRYPTO_PATH}/dh/dh_lib.c
		  ${CRYPTO_PATH}/dh/dh_check.c
		  ${CRYPTO_PATH}/dh/dh_err.c
		  ${CRYPTO_PATH}/dh/dh_depr.c
		  ${CRYPTO_PATH}/dh/dh_ameth.c 
		  ${CRYPTO_PATH}/dh/dh_pmeth.c 
		  ${CRYPTO_PATH}/dh/dh_prn.c
		)

		SET(CRYPTO_LIB_DSA_SRCS
		  ${CRYPTO_PATH}/dsa/dsa_gen.c
		  ${CRYPTO_PATH}/dsa/dsa_key.c
		  ${CRYPTO_PATH}/dsa/dsa_lib.c
		  ${CRYPTO_PATH}/dsa/dsa_asn1.c
		  ${CRYPTO_PATH}/dsa/dsa_vrf.c
		  ${CRYPTO_PATH}/dsa/dsa_sign.c
		  ${CRYPTO_PATH}/dsa/dsa_err.c
		  ${CRYPTO_PATH}/dsa/dsa_ossl.c
		  ${CRYPTO_PATH}/dsa/dsa_depr.c
		  ${CRYPTO_PATH}/dsa/dsa_ameth.c
		  ${CRYPTO_PATH}/dsa/dsa_pmeth.c
		  ${CRYPTO_PATH}/dsa/dsa_prn.c
		)

		SET(CRYPTO_LIB_DSO_SRCS
		  ${CRYPTO_PATH}/dso/dso_dl.c
		  ${CRYPTO_PATH}/dso/dso_dlfcn.c
		  ${CRYPTO_PATH}/dso/dso_err.c
		  ${CRYPTO_PATH}/dso/dso_lib.c
		  ${CRYPTO_PATH}/dso/dso_null.c
		  ${CRYPTO_PATH}/dso/dso_openssl.c
		  ${CRYPTO_PATH}/dso/dso_win32.c
		  ${CRYPTO_PATH}/dso/dso_vms.c
		  ${CRYPTO_PATH}/dso/dso_beos.c
		)

		SET(CRYPTO_LIB_EC_SRCS
		  ${CRYPTO_PATH}/ec/ec_lib.c
		  ${CRYPTO_PATH}/ec/ecp_smpl.c
		  ${CRYPTO_PATH}/ec/ecp_mont.c
		  ${CRYPTO_PATH}/ec/ecp_nist.c
		  ${CRYPTO_PATH}/ec/ec_cvt.c
		  ${CRYPTO_PATH}/ec/ec_mult.c
		  ${CRYPTO_PATH}/ec/ec_err.c
		  ${CRYPTO_PATH}/ec/ec_curve.c
		  ${CRYPTO_PATH}/ec/ec_check.c
		  ${CRYPTO_PATH}/ec/ec_print.c
		  ${CRYPTO_PATH}/ec/ec_asn1.c
		  ${CRYPTO_PATH}/ec/ec_key.c
		  ${CRYPTO_PATH}/ec/ec2_smpl.c
		  ${CRYPTO_PATH}/ec/ec2_mult.c
		  ${CRYPTO_PATH}/ec/ec_ameth.c
		  ${CRYPTO_PATH}/ec/ec_pmeth.c
		  ${CRYPTO_PATH}/ec/eck_prn.c
		)

		SET(CRYPTO_LIB_ECDH_SRCS
		  ${CRYPTO_PATH}/ecdh/ech_lib.c 
		  ${CRYPTO_PATH}/ecdh/ech_ossl.c
		  ${CRYPTO_PATH}/ecdh/ech_key.c
		  ${CRYPTO_PATH}/ecdh/ech_err.c
		)

		SET(CRYPTO_LIB_ECDSA_SRCS
		  ${CRYPTO_PATH}/ecdsa/ecs_lib.c
		  ${CRYPTO_PATH}/ecdsa/ecs_asn1.c
		  ${CRYPTO_PATH}/ecdsa/ecs_ossl.c
		  ${CRYPTO_PATH}/ecdsa/ecs_sign.c
		  ${CRYPTO_PATH}/ecdsa/ecs_vrf.c
		  ${CRYPTO_PATH}/ecdsa/ecs_err.c
		)

		SET(CRYPTO_LIB_RSA_SRCS
		  ${CRYPTO_PATH}/rsa/rsa_eay.c
		  ${CRYPTO_PATH}/rsa/rsa_gen.c
		  ${CRYPTO_PATH}/rsa/rsa_lib.c
		  ${CRYPTO_PATH}/rsa/rsa_sign.c
		  ${CRYPTO_PATH}/rsa/rsa_saos.c
		  ${CRYPTO_PATH}/rsa/rsa_err.c
		  ${CRYPTO_PATH}/rsa/rsa_pss.c
		  ${CRYPTO_PATH}/rsa/rsa_ssl.c
		  ${CRYPTO_PATH}/rsa/rsa_none.c
		  ${CRYPTO_PATH}/rsa/rsa_oaep.c
		  ${CRYPTO_PATH}/rsa/rsa_chk.c
		  ${CRYPTO_PATH}/rsa/rsa_null.c
		  ${CRYPTO_PATH}/rsa/rsa_pss.c
		  ${CRYPTO_PATH}/rsa/rsa_x931.c
		  ${CRYPTO_PATH}/rsa/rsa_asn1.c
		  ${CRYPTO_PATH}/rsa/rsa_depr.c
		  ${CRYPTO_PATH}/rsa/rsa_ameth.c
		  ${CRYPTO_PATH}/rsa/rsa_prn.c
		  ${CRYPTO_PATH}/rsa/rsa_pmeth.c
		  ${CRYPTO_PATH}/rsa/rsa_pk1.c
		)

		SET(CRYPTO_LIB_SHA_SRCS
		  ${CRYPTO_PATH}/sha/sha_dgst.c 
		  ${CRYPTO_PATH}/sha/sha1dgst.c 
		  ${CRYPTO_PATH}/sha/sha_one.c 
		  ${CRYPTO_PATH}/sha/sha1_one.c 
		  ${CRYPTO_PATH}/sha/sha256.c 
		  ${CRYPTO_PATH}/sha/sha512.c
		)

		SET(CRYPTO_LIB_SEED_SRCS
		  ${CRYPTO_PATH}/seed/seed.c
		  ${CRYPTO_PATH}/seed/seed_ecb.c
		  ${CRYPTO_PATH}/seed/seed_cbc.c
		  ${CRYPTO_PATH}/seed/seed_cfb.c
		  ${CRYPTO_PATH}/seed/seed_ofb.c
		)

		SET(CRYPTO_LIB_BIO_SRCS
		  ${CRYPTO_PATH}/bio/bio_lib.c 
		  ${CRYPTO_PATH}/bio/bio_cb.c 
		  ${CRYPTO_PATH}/bio/bio_err.c
		  ${CRYPTO_PATH}/bio/bss_mem.c 
		  ${CRYPTO_PATH}/bio/bss_null.c 
		  ${CRYPTO_PATH}/bio/bss_fd.c
		  ${CRYPTO_PATH}/bio/bss_file.c
		  ${CRYPTO_PATH}/bio/bss_sock.c
		  ${CRYPTO_PATH}/bio/bss_conn.c
		  ${CRYPTO_PATH}/bio/bf_null.c 
		  ${CRYPTO_PATH}/bio/bf_buff.c 
		  ${CRYPTO_PATH}/bio/b_print.c 
		  ${CRYPTO_PATH}/bio/b_dump.c
		  ${CRYPTO_PATH}/bio/b_sock.c 
		  ${CRYPTO_PATH}/bio/bss_acpt.c 
		  ${CRYPTO_PATH}/bio/bf_nbio.c 
		  ${CRYPTO_PATH}/bio/bss_log.c 
		  ${CRYPTO_PATH}/bio/bss_bio.c
		  ${CRYPTO_PATH}/bio/bss_dgram.c
		)

		SET(CRYPTO_LIB_OBJECTS_SRCS
		  ${CRYPTO_PATH}/objects/o_names.c 
		  ${CRYPTO_PATH}/objects/obj_dat.c 
		  ${CRYPTO_PATH}/objects/obj_lib.c 
		  ${CRYPTO_PATH}/objects/obj_err.c 
		  ${CRYPTO_PATH}/objects/obj_xref.c
		)

		SET(CRYPTO_LIB_ENGINE_SRCS
		  ${CRYPTO_PATH}/engine/eng_err.c
		  ${CRYPTO_PATH}/engine/eng_lib.c
		  ${CRYPTO_PATH}/engine/eng_list.c
		  ${CRYPTO_PATH}/engine/eng_init.c
		  ${CRYPTO_PATH}/engine/eng_ctrl.c
		  ${CRYPTO_PATH}/engine/eng_table.c
		  ${CRYPTO_PATH}/engine/eng_pkey.c
		  ${CRYPTO_PATH}/engine/eng_fat.c
		  ${CRYPTO_PATH}/engine/eng_all.c
		  ${CRYPTO_PATH}/engine/tb_rsa.c
		  ${CRYPTO_PATH}/engine/tb_dsa.c
		  ${CRYPTO_PATH}/engine/tb_ecdsa.c
		  ${CRYPTO_PATH}/engine/tb_dh.c
		  ${CRYPTO_PATH}/engine/tb_ecdh.c
		  ${CRYPTO_PATH}/engine/tb_rand.c
		  ${CRYPTO_PATH}/engine/tb_store.c
		  ${CRYPTO_PATH}/engine/tb_cipher.c
		  ${CRYPTO_PATH}/engine/tb_digest.c
		  ${CRYPTO_PATH}/engine/tb_pkmeth.c
		  ${CRYPTO_PATH}/engine/tb_asnmth.c
		)

		SET(CRYPTO_LIB_ERR_SRCS
		  ${CRYPTO_PATH}/err/err.c 
		  ${CRYPTO_PATH}/err/err_all.c 
		  ${CRYPTO_PATH}/err/err_prn.c
		)

		SET(CRYPTO_LIB_EVP_SRCS
		  ${CRYPTO_PATH}/evp/encode.c 
		  ${CRYPTO_PATH}/evp/digest.c 
		  ${CRYPTO_PATH}/evp/evp_enc.c 
		  ${CRYPTO_PATH}/evp/evp_key.c 
		  ${CRYPTO_PATH}/evp/evp_acnf.c
		  ${CRYPTO_PATH}/evp/e_des.c
		  ${CRYPTO_PATH}/evp/e_bf.c
		  ${CRYPTO_PATH}/evp/e_idea.c
		  ${CRYPTO_PATH}/evp/e_des3.c
		  ${CRYPTO_PATH}/evp/e_camellia.c
		  ${CRYPTO_PATH}/evp/e_rc4.c
		  ${CRYPTO_PATH}/evp/e_aes.c
		  ${CRYPTO_PATH}/evp/names.c
		  ${CRYPTO_PATH}/evp/e_seed.c
		  ${CRYPTO_PATH}/evp/e_xcbc_d.c
		  ${CRYPTO_PATH}/evp/e_rc2.c
		  ${CRYPTO_PATH}/evp/e_cast.c
		  ${CRYPTO_PATH}/evp/e_rc5.c
		  ${CRYPTO_PATH}/evp/m_null.c
		  ${CRYPTO_PATH}/evp/m_md2.c
		  ${CRYPTO_PATH}/evp/m_md4.c
		  ${CRYPTO_PATH}/evp/m_md5.c
		  ${CRYPTO_PATH}/evp/m_sha.c
		  ${CRYPTO_PATH}/evp/m_sha1.c
		  ${CRYPTO_PATH}/evp/m_wp.c
		  ${CRYPTO_PATH}/evp/m_dss.c
		  ${CRYPTO_PATH}/evp/m_dss1.c
		  ${CRYPTO_PATH}/evp/m_mdc2.c
		  ${CRYPTO_PATH}/evp/m_ripemd.c
		  ${CRYPTO_PATH}/evp/m_ecdsa.c
		  ${CRYPTO_PATH}/evp/p_open.c
		  ${CRYPTO_PATH}/evp/p_seal.c
		  ${CRYPTO_PATH}/evp/p_sign.c
		  ${CRYPTO_PATH}/evp/p_verify.c
		  ${CRYPTO_PATH}/evp/p_lib.c
		  ${CRYPTO_PATH}/evp/p_enc.c
		  ${CRYPTO_PATH}/evp/p_dec.c
		  ${CRYPTO_PATH}/evp/bio_md.c
		  ${CRYPTO_PATH}/evp/bio_b64.c
		  ${CRYPTO_PATH}/evp/bio_enc.c
		  ${CRYPTO_PATH}/evp/evp_err.c
		  ${CRYPTO_PATH}/evp/e_null.c
		  ${CRYPTO_PATH}/evp/c_all.c
		  ${CRYPTO_PATH}/evp/c_allc.c
		  ${CRYPTO_PATH}/evp/c_alld.c
		  ${CRYPTO_PATH}/evp/evp_lib.c
		  ${CRYPTO_PATH}/evp/bio_ok.c
		  ${CRYPTO_PATH}/evp/evp_pkey.c
		  ${CRYPTO_PATH}/evp/evp_pbe.c
		  ${CRYPTO_PATH}/evp/p5_crpt.c
		  ${CRYPTO_PATH}/evp/p5_crpt2.c
		  ${CRYPTO_PATH}/evp/e_old.c
		  ${CRYPTO_PATH}/evp/pmeth_lib.c
		  ${CRYPTO_PATH}/evp/pmeth_fn.c
		  ${CRYPTO_PATH}/evp/pmeth_gn.c
		  ${CRYPTO_PATH}/evp/m_sigver.c
		)

		SET(CRYPTO_LIB_MD5_SRCS
		  ${CRYPTO_PATH}/md5/md5_dgst.c 
		  ${CRYPTO_PATH}/md5/md5_one.c
		)

		SET(CRYPTO_LIB_MODES_SRCS
		  ${CRYPTO_PATH}/modes/cbc128.c 
		  ${CRYPTO_PATH}/modes/ctr128.c 
		  ${CRYPTO_PATH}/modes/cts128.c 
		  ${CRYPTO_PATH}/modes/cfb128.c 
		  ${CRYPTO_PATH}/modes/ofb128.c
		)

		SET(CRYPTO_LIB_HMAC_SRCS
		  ${CRYPTO_PATH}/hmac/hmac.c
		  ${CRYPTO_PATH}/hmac/hm_ameth.c
		  ${CRYPTO_PATH}/hmac/hm_pmeth.c
		)

		SET(CRYPTO_LIB_LHASH_SRCS
		  ${CRYPTO_PATH}/lhash/lhash.c 
		  ${CRYPTO_PATH}/lhash/lh_stats.c
		)

		SET(CRYPTO_LIB_OCSP_SRCS
		  ${CRYPTO_PATH}/ocsp/ocsp_asn.c
		  ${CRYPTO_PATH}/ocsp/ocsp_ext.c
		  ${CRYPTO_PATH}/ocsp/ocsp_ht.c
		  ${CRYPTO_PATH}/ocsp/ocsp_lib.c
		  ${CRYPTO_PATH}/ocsp/ocsp_cl.c
		  ${CRYPTO_PATH}/ocsp/ocsp_srv.c
		  ${CRYPTO_PATH}/ocsp/ocsp_prn.c
		  ${CRYPTO_PATH}/ocsp/ocsp_vfy.c
		  ${CRYPTO_PATH}/ocsp/ocsp_err.c
		)

		SET(CRYPTO_LIB_PEM_SRCS
		  ${CRYPTO_PATH}/pem/pem_sign.c
		  ${CRYPTO_PATH}/pem/pem_seal.c
		  ${CRYPTO_PATH}/pem/pem_info.c
		  ${CRYPTO_PATH}/pem/pem_lib.c
		  ${CRYPTO_PATH}/pem/pem_all.c
		  ${CRYPTO_PATH}/pem/pem_err.c
		  ${CRYPTO_PATH}/pem/pem_x509.c
		  ${CRYPTO_PATH}/pem/pem_xaux.c
		  ${CRYPTO_PATH}/pem/pem_oth.c
		  ${CRYPTO_PATH}/pem/pem_pk8.c
		  ${CRYPTO_PATH}/pem/pem_pkey.c
		  ${CRYPTO_PATH}/pem/pvkfmt.c
		)

		SET(CRYPTO_LIB_PKCS7_SRCS
		  ${CRYPTO_PATH}/pkcs7/pk7_asn1.c
		  ${CRYPTO_PATH}/pkcs7/pk7_lib.c
		  ${CRYPTO_PATH}/pkcs7/pkcs7err.c
		  ${CRYPTO_PATH}/pkcs7/pk7_doit.c
		  ${CRYPTO_PATH}/pkcs7/pk7_smime.c
		  ${CRYPTO_PATH}/pkcs7/pk7_attr.c
		  ${CRYPTO_PATH}/pkcs7/pk7_mime.c
		  ${CRYPTO_PATH}/pkcs7/bio_pk7.c
		)

		SET(CRYPTO_LIB_PKCS12_SRCS
		  ${CRYPTO_PATH}/pkcs12/p12_add.c
		  ${CRYPTO_PATH}/pkcs12/p12_asn.c
		  ${CRYPTO_PATH}/pkcs12/p12_attr.c
		  ${CRYPTO_PATH}/pkcs12/p12_crpt.c
		  ${CRYPTO_PATH}/pkcs12/p12_crt.c
		  ${CRYPTO_PATH}/pkcs12/p12_decr.c
		  ${CRYPTO_PATH}/pkcs12/p12_init.c
		  ${CRYPTO_PATH}/pkcs12/p12_key.c
		  ${CRYPTO_PATH}/pkcs12/p12_kiss.c
		  ${CRYPTO_PATH}/pkcs12/p12_mutl.c
		  ${CRYPTO_PATH}/pkcs12/p12_utl.c
		  ${CRYPTO_PATH}/pkcs12/p12_npas.c
		  ${CRYPTO_PATH}/pkcs12/pk12err.c
		  ${CRYPTO_PATH}/pkcs12/p12_p8d.c
		  ${CRYPTO_PATH}/pkcs12/p12_p8e.c
		)

		SET(CRYPTO_LIB_RAND_SRCS
		  ${CRYPTO_PATH}/rand/md_rand.c
		  ${CRYPTO_PATH}/rand/randfile.c
		  ${CRYPTO_PATH}/rand/rand_lib.c
		  ${CRYPTO_PATH}/rand/rand_err.c
		  ${CRYPTO_PATH}/rand/rand_egd.c
		  ${CRYPTO_PATH}/rand/rand_win.c
		  ${CRYPTO_PATH}/rand/rand_unix.c
		  ${CRYPTO_PATH}/rand/rand_os2.c
		  ${CRYPTO_PATH}/rand/rand_nw.c
		)

		SET(CRYPTO_LIB_STACK_SRCS
		  ${CRYPTO_PATH}/stack/stack.c
		)

		SET(CRYPTO_LIB_TS_SRCS
		  ${CRYPTO_PATH}/ts/ts_err.c
		  ${CRYPTO_PATH}/ts/ts_req_utils.c
		  ${CRYPTO_PATH}/ts/ts_req_print.c
		  ${CRYPTO_PATH}/ts/ts_rsp_utils.c
		  ${CRYPTO_PATH}/ts/ts_rsp_print.c
		  ${CRYPTO_PATH}/ts/ts_rsp_sign.c 
		  ${CRYPTO_PATH}/ts/ts_rsp_verify.c
		  ${CRYPTO_PATH}/ts/ts_verify_ctx.c 
		  ${CRYPTO_PATH}/ts/ts_lib.c 
		  ${CRYPTO_PATH}/ts/ts_conf.c
		  ${CRYPTO_PATH}/ts/ts_asn1.c
		)

		SET(CRYPTO_LIB_UI_SRCS
		  ${CRYPTO_PATH}/ui/ui_err.c
		  ${CRYPTO_PATH}/ui/ui_lib.c
		  ${CRYPTO_PATH}/ui/ui_openssl.c
		  ${CRYPTO_PATH}/ui/ui_util.c
		  ${CRYPTO_PATH}/ui/ui_compat.c
		)

		SET(CRYPTO_LIB_X509_SRCS
		  ${CRYPTO_PATH}/x509/x509_def.c
		  ${CRYPTO_PATH}/x509/x509_d2.c
		  ${CRYPTO_PATH}/x509/x509_r2x.c
		  ${CRYPTO_PATH}/x509/x509_cmp.c
		  ${CRYPTO_PATH}/x509/x509_obj.c
		  ${CRYPTO_PATH}/x509/x509_req.c
		  ${CRYPTO_PATH}/x509/x509spki.c
		  ${CRYPTO_PATH}/x509/x509_vfy.c
		  ${CRYPTO_PATH}/x509/x509_set.c
		  ${CRYPTO_PATH}/x509/x509cset.c
		  ${CRYPTO_PATH}/x509/x509rset.c
		  ${CRYPTO_PATH}/x509/x509_err.c
		  ${CRYPTO_PATH}/x509/x509name.c
		  ${CRYPTO_PATH}/x509/x509_v3.c
		  ${CRYPTO_PATH}/x509/x509_ext.c
		  ${CRYPTO_PATH}/x509/x509_att.c
		  ${CRYPTO_PATH}/x509/x509type.c
		  ${CRYPTO_PATH}/x509/x509_lu.c
		  ${CRYPTO_PATH}/x509/x_all.c
		  ${CRYPTO_PATH}/x509/x509_txt.c
		  ${CRYPTO_PATH}/x509/x509_trs.c
		  ${CRYPTO_PATH}/x509/by_file.c
		  ${CRYPTO_PATH}/x509/by_dir.c
		  ${CRYPTO_PATH}/x509/x509_vpm.c
		)

		SET(CRYPTO_LIB_X509V3_SRCS
		  ${CRYPTO_PATH}/x509v3/v3_bcons.c
		  ${CRYPTO_PATH}/x509v3/v3_bitst.c
		  ${CRYPTO_PATH}/x509v3/v3_conf.c
		  ${CRYPTO_PATH}/x509v3/v3_extku.c
		  ${CRYPTO_PATH}/x509v3/v3_ia5.c
		  ${CRYPTO_PATH}/x509v3/v3_lib.c
		  ${CRYPTO_PATH}/x509v3/v3_prn.c
		  ${CRYPTO_PATH}/x509v3/v3_utl.c
		  ${CRYPTO_PATH}/x509v3/v3err.c
		  ${CRYPTO_PATH}/x509v3/v3_genn.c
		  ${CRYPTO_PATH}/x509v3/v3_alt.c
		  ${CRYPTO_PATH}/x509v3/v3_skey.c
		  ${CRYPTO_PATH}/x509v3/v3_akey.c
		  ${CRYPTO_PATH}/x509v3/v3_pku.c
		  ${CRYPTO_PATH}/x509v3/v3_int.c
		  ${CRYPTO_PATH}/x509v3/v3_enum.c
		  ${CRYPTO_PATH}/x509v3/v3_sxnet.c
		  ${CRYPTO_PATH}/x509v3/v3_cpols.c
		  ${CRYPTO_PATH}/x509v3/v3_crld.c
		  ${CRYPTO_PATH}/x509v3/v3_purp.c
		  ${CRYPTO_PATH}/x509v3/v3_info.c
		  ${CRYPTO_PATH}/x509v3/v3_ocsp.c
		  ${CRYPTO_PATH}/x509v3/v3_akeya.c
		  ${CRYPTO_PATH}/x509v3/v3_pmaps.c
		  ${CRYPTO_PATH}/x509v3/v3_pcons.c
		  ${CRYPTO_PATH}/x509v3/v3_ncons.c
		  ${CRYPTO_PATH}/x509v3/v3_pcia.c
		  ${CRYPTO_PATH}/x509v3/v3_pci.c
		  ${CRYPTO_PATH}/x509v3/pcy_cache.c
		  ${CRYPTO_PATH}/x509v3/pcy_node.c
		  ${CRYPTO_PATH}/x509v3/pcy_data.c
		  ${CRYPTO_PATH}/x509v3/pcy_map.c
		  ${CRYPTO_PATH}/x509v3/pcy_tree.c
		  ${CRYPTO_PATH}/x509v3/pcy_lib.c
		  ${CRYPTO_PATH}/x509v3/v3_asid.c 
		  ${CRYPTO_PATH}/x509v3/v3_addr.c
		)


		INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
		INCLUDE_DIRECTORIES(${OPENSSL_SOURCE_DIR})
		INCLUDE_DIRECTORIES(${OPENSSL_SOURCE_DIR}/ssl)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH})
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/eas1)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/asn1)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/bf)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/bio)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/bn)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/buffer)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/camelia)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/cast)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/cms)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/comp)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/conf)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/des)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/dh)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/dsa)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/dso)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/ec)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/ecdh)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/ecdsa)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/engine)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/err)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/evp)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/hmac)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/idea)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/jpake)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/kbr5)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/lhash)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/md2)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/md4)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/md5)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/mdc2)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/modes)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/objects)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/ocsp)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/pem)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/pkc7)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/pkc12)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/pqueue)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/rand)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/rc2)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/rc4)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/rc5)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/ripemd)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/rsa)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/sha)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/seed)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/stack)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/store)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/threads)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/ts)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/txt_db)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/ui)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/whrlpool)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/x509)
		INCLUDE_DIRECTORIES(${CRYPTO_PATH}/x509v3)
		INCLUDE_DIRECTORIES(${OPENSSL_SOURCE_DIR}/include)
		INCLUDE_DIRECTORIES(${SCI_ZLIB_INCLUDE})

		SET(SCI_CRYPTO_LIBRARY crypto CACHE INTERNAL "crypto library name." FORCE)

		ADD_DEFINITIONS(-DOPENSSLDIR="${CRYPTO_PATH}/..")
		ADD_DEFINITIONS(-DENGINESDIR="${CRYPTO_PATH}/engine")
		ADD_DEFINITIONS(-DZLIB)
		IF(APPLE)
		   ADD_DEFINITIONS(-DOPENSSL_NO_INLINE_ASM)
		ENDIF(APPLE)

		IF(DEFINED SCI_ZLIB_MANGLED)
		  ADD_DEFINITIONS(-DZ_PREFIX)
		ENDIF(DEFINED SCI_ZLIB_MANGLED)

		TEST_BIG_ENDIAN(BIG_ENDIAN)

		IF(BIG_ENDIAN)
		  ADD_DEFINITIONS(-DB_ENDIAN)
		ELSE(BIG_ENDIAN)
		  ADD_DEFINITIONS(-DL_ENDIAN)
		ENDIF(BIG_ENDIAN)

		IF(WIN32)
		  ADD_DEFINITIONS(-DWINDOWS)
		  ADD_DEFINITIONS(-DOPENSSL_IMPLEMENTS_strncasecmp)
		  ADD_DEFINITIONS(-DOPENSSL_SYSNAME_WIN32)
		  ADD_DEFINITIONS(-DWIN32_LEAN_AND_MEAN)
		  ADD_DEFINITIONS(-D_BIND_TO_CURRENT_VCLIBS_VERSION=1)
		ENDIF(WIN32)

		CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/buildinf.h.in ${CMAKE_CURRENT_BINARY_DIR}/buildinf.h)

		ADD_LIBRARY(${SCI_CRYPTO_LIBRARY} STATIC
		  ${CRYPTO_LIB_SRCS} 
		  ${CRYPTO_LIB_AES_SRCS} 
		  ${CRYPTO_LIB_ASN1_SRCS} 
		  ${CRYPTO_LIB_BIO_SRCS}
		  ${CRYPTO_LIB_BF_SRCS}
		  ${CRYPTO_LIB_BUFFER_SRCS}
		  ${CRYPTO_LIB_BN_SRCS}
		  ${CRYPTO_LIB_CMS_SRCS}
		  ${CRYPTO_LIB_CONF_SRCS}
		  ${CRYPTO_LIB_COMP_SRCS}
		  ${CRYPTO_LIB_DES_SRCS}
		  ${CRYPTO_LIB_DH_SRCS}
		  ${CRYPTO_LIB_DSA_SRCS}
		  ${CRYPTO_LIB_DSO_SRCS}
		  ${CRYPTO_LIB_OBJECTS_SRCS}
		  ${CRYPTO_LIB_ERR_SRCS}
		  ${CRYPTO_LIB_ENGINE_SRCS}
		  ${CRYPTO_LIB_EVP_SRCS}
		  ${CRYPTO_LIB_EC_SRCS}
		  ${CRYPTO_LIB_ECDH_SRCS}
		  ${CRYPTO_LIB_ECDSA_SRCS}
		  ${CRYPTO_LIB_HMAC_SRCS}
		  ${CRYPTO_LIB_MD5_SRCS}
		  ${CRYPTO_LIB_MODES_SRCS}
		  ${CRYPTO_LIB_LHASH_SRCS}  
		  ${CRYPTO_LIB_OCSP_SRCS}  
		  ${CRYPTO_LIB_PEM_SRCS}  
		  ${CRYPTO_LIB_PKCS7_SRCS}  
		  ${CRYPTO_LIB_PKCS12_SRCS}  
		  ${CRYPTO_LIB_RAND_SRCS}  
		  ${CRYPTO_LIB_RSA_SRCS}
		  ${CRYPTO_LIB_SHA_SRCS}
		  ${CRYPTO_LIB_SEED_SRCS}
		  ${CRYPTO_LIB_STACK_SRCS}
		  ${CRYPTO_LIB_TS_SRCS}
		  ${CRYPTO_LIB_UI_SRCS}
		  ${CRYPTO_LIB_X509_SRCS}
		  ${CRYPTO_LIB_X509V3_SRCS}
		)

		TARGET_LINK_LIBRARIES(${SCI_CRYPTO_LIBRARY} ${SCI_ZLIB_LIBRARY})  

		IF(WIN32)
		  SET_TARGET_PROPERTIES(${SCI_CRYPTO_LIBRARY}
								PROPERTIES
								COMPILE_FLAGS "-W3 -Gs0 -GF -Gy -nologo -wd4267 -wd4244 -wd4334")
		ENDIF(WIN32)

		IF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
		  SET(ENABLE_PIC ON)
		ENDIF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)

		IF(CMAKE_COMPILER_IS_GNUCC)
		IF(ENABLE_PIC)
		  SET_TARGET_PROPERTIES(${SCI_CRYPTO_LIBRARY}
								PROPERTIES
								  COMPILE_FLAGS "-w -fPIC"
							   )
		ELSE(ENABLE_PIC)
		  SET_TARGET_PROPERTIES(${SCI_CRYPTO_LIBRARY}
								PROPERTIES
								COMPILE_FLAGS "-w")
		ENDIF(ENABLE_PIC)
		ENDIF(CMAKE_COMPILER_IS_GNUCC)

		SET(SCI_CRYPTO_INCLUDE ${CRYPTO_PATH}/../include;${CRYPTO_PATH} CACHE INTERNAL "crypto include directories." FORCE)


		SET(SSL_LIB_SRCS
			${SSL_PATH}/s2_meth.c
			${SSL_PATH}/s2_srvr.c
			${SSL_PATH}/s2_clnt.c
			${SSL_PATH}/s2_lib.c
			${SSL_PATH}/s2_enc.c
			${SSL_PATH}/s2_pkt.c
			${SSL_PATH}/s3_meth.c
			${SSL_PATH}/s3_srvr.c
			${SSL_PATH}/s3_clnt.c 
			${SSL_PATH}/s3_lib.c
			${SSL_PATH}/s3_enc.c
			${SSL_PATH}/s3_pkt.c
			${SSL_PATH}/s3_both.c
			${SSL_PATH}/s23_meth.c 
			${SSL_PATH}/s23_srvr.c
			${SSL_PATH}/s23_clnt.c
			${SSL_PATH}/s23_lib.c
			${SSL_PATH}/s23_pkt.c
			${SSL_PATH}/t1_meth.c
			${SSL_PATH}/t1_srvr.c
			${SSL_PATH}/t1_clnt.c
			${SSL_PATH}/t1_lib.c
			${SSL_PATH}/t1_enc.c
			${SSL_PATH}/d1_meth.c
			${SSL_PATH}/d1_srvr.c
			${SSL_PATH}/d1_clnt.c
			${SSL_PATH}/d1_lib.c
			${SSL_PATH}/d1_pkt.c
			${SSL_PATH}/d1_both.c
			${SSL_PATH}/d1_enc.c
			${SSL_PATH}/ssl_lib.c
			${SSL_PATH}/ssl_err2.c
			${SSL_PATH}/ssl_cert.c
			${SSL_PATH}/ssl_sess.c
			${SSL_PATH}/ssl_ciph.c
			${SSL_PATH}/ssl_stat.c
			${SSL_PATH}/ssl_rsa.c
			${SSL_PATH}/ssl_asn1.c
			${SSL_PATH}/ssl_txt.c
			${SSL_PATH}/ssl_algs.c
			${SSL_PATH}/bio_ssl.c
			${SSL_PATH}/ssl_err.c
			${SSL_PATH}/kssl.c
			${SSL_PATH}/t1_reneg.c)



		SET(SCI_SSL_LIBRARY "ssl;${SCI_CRYPTO_LIBRARY}" CACHE INTERNAL "ssl library name." FORCE)

		ADD_DEFINITIONS(-DOPENSSLDIR="${SSL_PATH}/..")
		ADD_DEFINITIONS(-DENGINESDIR="${SSL_PATH}/engine")
		ADD_DEFINITIONS(-DZLIB)
		ADD_DEFINITIONS(-DZ_PREFIX)

		TEST_BIG_ENDIAN(BIG_ENDIAN)

		IF(BIG_ENDIAN)
		  ADD_DEFINITIONS(-DB_ENDIAN)
		ELSE(BIG_ENDIAN)
		  ADD_DEFINITIONS(-DL_ENDIAN)
		ENDIF(BIG_ENDIAN)

		IF(WIN32)
		  ADD_DEFINITIONS(-DOPENSSL_IMPLEMENTS_strncasecmp)
		  ADD_DEFINITIONS(-DOPENSSL_SYSNAME_WIN32)
		  ADD_DEFINITIONS(-DWIN32_LEAN_AND_MEAN)
		  ADD_DEFINITIONS(-D_BIND_TO_CURRENT_VCLIBS_VERSION=1)
		ENDIF(WIN32)

		ADD_LIBRARY(ssl STATIC
		  ${SSL_LIB_SRCS} 
		)

		IF(WIN32)
		  SET_TARGET_PROPERTIES(ssl
								PROPERTIES
								COMPILE_FLAGS "-W3 -Gs0 -GF -Gy -nologo -wd4267 -wd4244 -wd4334")
		ENDIF(WIN32)

		SET(SCI_SSL_INCLUDE ${OPENSSL_SOURCE_DIR}/include CACHE INTERNAL "crypto include directories." FORCE)
		SET(SCI_SSL_FOUND TRUE CACHE INTERNAL "crypto library found" FORCE)
	ELSE( EXISTS ${OPENSSL_SOURCE_DIR}/crypto )
		MESSAGE( WARNING "Warning: Could not find openssl, disabling SSL options.")
		MESSAGE( WARNING "Warning: To enable SSL please download openssl 1.0 source code and extract the code into 'src/Externals/ssl/openssl_code/'. Once the code is extracted there, please reconfigure using cmake.")
	ENDIF( EXISTS ${OPENSSL_SOURCE_DIR}/crypto )

ENDIF( OPENSSL_FOUND )

