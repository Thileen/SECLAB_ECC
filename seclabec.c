/*



*/
#include <stdio.h>
#include <stdlib.h>

#include "ec_domain_parameter.h"
//#include "seclabec.h"

void shiftArrayLeft_fast_reduction(bit_32 *a, bit_16 len, bit_16 *carry);
void shiftArrayRightBySteps(bit_32 *a, bit_16 len, bit_16 steps);
void shiftArrayLeftBySteps(bit_32 *a, bit_16 len, bit_16 steps);
void gmp_add(bit_32 *a, bit_32 *b, bit_32 *c, bit_16 *carry, bit_16 len);
void gmp_sub(bit_32 *a, bit_32 *b, bit_32 *c, bit_16 *carry, bit_16 len);
void gmp_add_mod(bit_32 *a, bit_32 *b, bit_32 *c);
void gmp_sub_mod(bit_32 *a, bit_32 *b, bit_32 *c);
void gmp_mul(bit_32 *a, bit_32 *b, bit_32 *c, bit_16 len);
void gmp_mul_mod(bit_32 *a, bit_32 *b, bit_32 *c);
void omp_div(bit_32 *a, bit_16 lenA, bit_32 *b, bit_16 lenB, bit_32 *c, bit_32 *r);
void br_mod(bit_32 *a, bit_32 *c);
void fast_reduction_p256(bit_32 *a,bit_32 *c);
bit_16 compareAA(bit_32 *a, bit_32 *b, bit_16 len);
bit_16 compareAI(bit_32 *a, bit_32 b, bit_16 len);

static const DOMAIN_PARAMETER *dp;

// var
bit_16 i;

// memory
char keySize;
char keySizeM1;
char keySizeP1;
char keySize2;
char keySize4;
bit_16 pMSBpos;

// add_mod
bit_32 *add_tmp1, *add_tmp2, *add_tmp3;
// sub_mod
bit_32 *sub_tmp;
// mul_mod
bit_32 *mul_tmp1, *mul_tmp2;	
// fast reduction
bit_32 *s1, *s2, *s3, *s4, *s5, *s6, *s7, *s8, *s9;
// barrett reduction
// br_mod
bit_32 *br_q, *br_tmp1111, *br_tmp11, *br_tmp22, *br_tmp33, *br_tmp44;
bit_32 br_k;
// point affine and edward
bit_32 *p_s, *p_tmp1, *p_tmp2, *p_tmp3, *p_tmp4, *p_tmp5;
// shift left mod
bit_32 *sl_tmp1, *sl_tmp2;
// montgomery inversion 
bit_32 *mi_u, *mi_v, *mi_x1, *mi_x2, *mi_tmp;
// point projective
bit_32 *p_u, *p_uu, *p_v, *p_R, *p_A; 
// edward curve point in projective plane
bit_32 *p_B, *p_C, *p_D, *p_E, *p_F, *p_G;


// konstanten

bit_32 *c1;
bit_32 *c3;

bit_32 *cp2; // primzahl in einem Array der größte keysize*2
	
	
// random permutation 
bit_16 rp[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

// array of points for sliding window
// index 0 is empty!!!
point *oArrPoint; // original array of points
point *arrPoint; // array of points
	
// barrett reduktion variablen
bit_32 *bkp1;
bit_32 *bkp1M1;	// bkp1-1
bit_32 *mikro;
 
static const DOMAIN_PARAMETER p_192 = {
	KEY_SIZE_192,																		// key_size
	(bit_32[]){0xffffffff, 0xffffffff, 0xffffffff, 0xfffffffe, 0xffffffff, 0xffffffff},	// p
	(bit_32[]){0xffffffff, 0xffffffff, 0xffffffff, 0xfffffffe, 0xffffffff, 0xfffffffc},	// a
	(bit_32[]){0x64210519, 0xe59c80e7, 0x0fa7e9ab, 0x72243049, 0xfeb8deec, 0xc146b9b1},	// b
	(bit_32[]){0x188da80e, 0xb03090f6, 0x7cbf20eb, 0x43a18800, 0xf4ff0afd, 0x82ff1012},	// x
	(bit_32[]){0x07192b95, 0xffc8da78, 0x631011ed, 0x6b24cdd5, 0x73f977a1, 0x1e794811},	// y
	(bit_32[]){0xffffffff, 0xffffffff, 0xffffffff, 0x99def836, 0x146bc9b1, 0xb4d22831},	// n 
	0x1};

static const DOMAIN_PARAMETER p_224 = {
	KEY_SIZE_224,																		// key_size
	(bit_32[]){0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 0x00000000, 0x00000001},	// p
	(bit_32[]){0xffffffff, 0xffffffff, 0xffffffff, 0xfffffffe, 0xffffffff, 0xffffffff, 0xfffffffe},	// a
	(bit_32[]){0xb4050a85, 0x0c04b3ab, 0xf5413256, 0x5044b0b7, 0xd7bfd8ba, 0x270b3943, 0x2355ffb4},	// b
	(bit_32[]){0xb70e0cbd, 0x6bb4bf7f, 0x321390b9, 0x4a03c1d3, 0x56c21122, 0x343280d6, 0x115c1d21},	// x
	(bit_32[]){0xbd376388, 0xb5f723fb, 0x4c22dfe6, 0xcd4375a0, 0x5a074764, 0x44d58199, 0x85007e34},	// y
	(bit_32[]){0xffffffff, 0xffffffff, 0xffffffff, 0xffff16a2, 0xe0b8f03e, 0x13dd2945, 0x5c5c2a3d},	// n 
	0x1};	
//http://csrc.nist.gov/publications/fips/fips186-3/fips_186-3.pdf
static const DOMAIN_PARAMETER p_256 = {
	KEY_SIZE_256,																		// key_size
	(bit_32[]){0xffffffff, 0x00000001, 0x00000000, 0x00000000, 0x00000000, 0xffffffff, 0xffffffff, 0xffffffff},	// p
	(bit_32[]){0xffffffff, 0x00000001, 0x00000000, 0x00000000, 0x00000000, 0xffffffff, 0xffffffff, 0xfffffffc},	// a
	(bit_32[]){0x5ac635d8, 0xaa3a93e7, 0xb3ebbd55, 0x769886bc, 0x651d06b0, 0xcc53b0f6, 0x3bce3c3e, 0x27d2604b},	// b
	(bit_32[]){0x6b17d1f2, 0xe12c4247, 0xf8bce6e5, 0x63a440f2, 0x77037d81, 0x2deb33a0, 0xf4a13945, 0xd898c296},	// x
	(bit_32[]){0x4fe342e2, 0xfe1a7f9b, 0x8ee7eb4a, 0x7c0f9e16, 0x2bce3357, 0x6b315ece, 0xcbb64068, 0x37bf51f5},	// y
	(bit_32[]){0xffffffff, 0x00000000, 0xffffffff, 0xffffffff, 0xbce6faad, 0xa7179e84, 0xf3b9cac2, 0xfc632551},	// n 
	0x1};	

static const DOMAIN_PARAMETER p_384 = {
	KEY_SIZE_384,																		// key_size
	(bit_32[]){0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xfffffffe, 0xffffffff, 0x00000000, 0x00000000, 0xffffffff},	// p
	(bit_32[]){0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xfffffffe, 0xffffffff, 0x00000000, 0x00000000, 0xfffffffc},	// a
	(bit_32[]){0xb3312fa7, 0xe23ee7e4, 0x988e056b, 0xe3f82d19, 0x181d9c6e, 0xfe814112, 0x0314088f, 0x5013875a, 0xc656398d, 0x8a2ed19d, 0x2a85c8ed, 0xd3ec2aef},	// b
	(bit_32[]){0xaa87ca22, 0xbe8b0537, 0x8eb1c71e, 0xf320ad74, 0x6e1d3b62, 0x8ba79b98, 0x59f741e0, 0x82542a38, 0x5502f25d, 0xbf55296c, 0x3a545e38, 0x72760aB7},	// x
	(bit_32[]){0x3617de4a, 0x96262c6f, 0x5d9e98bf, 0x9292dc29, 0xf8f41dbd, 0x289a147c, 0xe9da3113, 0xb5f0b8c0, 0x0a60b1ce, 0x1d7e819d, 0x7a431d7c, 0x90ea0e5F},	// y
	(bit_32[]){0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xc7634d81, 0xf4372ddf, 0x581a0db2, 0x48b0a77a, 0xecec196a, 0xccc52973},	// n 
	0x1};	

// 
static const DOMAIN_PARAMETER p_521 = {
	KEY_SIZE_521,
	(bit_32[]){0x000001ff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff},	// p
	(bit_32[]){0x000001ff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xfffffffc},	// a
	(bit_32[]){0x00000051, 0x953eb961, 0x8e1c9a1f, 0x929a21a0, 0xb68540ee, 0xa2da725b, 0x99b315f3, 0xb8b48991, 0x8ef109e1, 0x56193951, 0xec7e937b, 0x1652c0bd, 0x3bb1bf07, 0x3573df88, 0x3d2c34f1, 0xef451fd4, 0x6b503f00},	// b
	(bit_32[]){0x000000c6, 0x858e06b7, 0x0404e9cd, 0x9e3ecb66, 0x2395b442, 0x9c648139, 0x053fb521, 0xf828af60, 0x6b4d3dba, 0xa14b5e77, 0xefe75928, 0xfe1dc127, 0xa2ffa8de, 0x3348b3c1, 0x856a429b, 0xf97e7e31, 0xc2e5bd66},	// x
	(bit_32[]){0x00000118, 0x39296a78, 0x9a3bc004, 0x5c8a5fb4, 0x2c7d1bd9, 0x98f54449, 0x579b4468, 0x17afbd17, 0x273e662c, 0x97ee7299, 0x5ef42640, 0xc550b901, 0x3fad0761, 0x353c7086, 0xa272c240, 0x88be9476, 0x9fd16650},	// y	
	(bit_32[]){0x000001ff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xfffffffa, 0x51868783, 0xbf2f966b, 0x7fcc0148, 0xf709a5d0, 0x3bb5c9b8, 0x899c47ae, 0xbb6fb71e, 0x91386409},	// n
	0x01};
	
// https://tools.ietf.org/html/rfc5639
static const DOMAIN_PARAMETER brainpoolP256r1 = {
	KEY_SIZE_256,
	(bit_32[]){0xA9FB57DB, 0xA1EEA9BC, 0x3E660A90, 0x9D838D72, 0x6E3BF623, 0xD5262028, 0x2013481D, 0x1F6E5377}, // p
	(bit_32[]){0x7D5A0975, 0xFC2C3057, 0xEEF67530, 0x417AFFE7, 0xFB8055C1, 0x26DC5C6C, 0xE94A4B44, 0xF330B5D9}, // a
	(bit_32[]){0x26DC5C6C, 0xE94A4B44, 0xF330B5D9, 0xBBD77CBF, 0x95841629, 0x5CF7E1CE, 0x6BCCDC18, 0xFF8C07B6}, // b
	(bit_32[]){0x8BD2AEB9, 0xCB7E57CB, 0x2C4B482F, 0xFC81B7AF, 0xB9DE27E1, 0xE3BD23C2, 0x3A4453BD, 0x9ACE3262}, // x
	(bit_32[]){0x547EF835, 0xC3DAC4FD, 0x97F8461A, 0x14611DC9, 0xC2774513, 0x2DED8E54, 0x5C1D54C7, 0x2F046997}, // y
	(bit_32[]){0xa9fb57db, 0xa1eea9bc, 0x3e660a90, 0x9d838d71, 0x8c397aa3, 0xb561a6f7, 0x901e0e82, 0x974856a7}, // n
	0x1}; 		// h
	
static const DOMAIN_PARAMETER brainpoolP320r1 = {
	KEY_SIZE_320,
	(bit_32[]){0xD35E4720, 0x36BC4FB7, 0xE13C785E, 0xD201E065, 0xF98FCFA6, 0xF6F40DEF, 0x4F92B9EC, 0x7893EC28, 0xFCD412B1, 0xF1B32E27}, // p
	(bit_32[]){0x3EE30B56, 0x8FBAB0F8, 0x83CCEBD4, 0x6D3F3BB8, 0xA2A73513, 0xF5EB79DA, 0x66190EB0, 0x85FFA9F4, 0x92F375A9, 0x7D860EB4}, // a
	(bit_32[]){0x52088394, 0x9DFDBC42, 0xD3AD1986, 0x40688A6F, 0xE13F4134, 0x9554B49A, 0xCC31DCCD, 0x88453981, 0x6F5EB4AC, 0x8FB1F1A6}, // b
	(bit_32[]){0x43BD7E9A, 0xFB53D8B8, 0x5289BCC4, 0x8EE5BFE6, 0xF20137D1, 0x0A087EB6, 0xE7871E2A, 0x10A599C7, 0x10AF8D0D, 0x39E20611}, // x
	(bit_32[]){0x14FDD055, 0x45EC1CC8, 0xAB409324, 0x7F77275E, 0x0743FFED, 0x117182EA, 0xA9C77877, 0xAAAC6AC7, 0xD35245D1, 0x692E8EE1}, // y
	(bit_32[]){0xD35E4720, 0x36BC4FB7, 0xE13C785E, 0xD201E065, 0xF98FCFA5, 0xB68F12A3, 0x2D482EC7, 0xEE8658E9, 0x8691555B, 0x44C59311}, // n
	0x1}; // h

static const DOMAIN_PARAMETER brainpoolP384r1 = {
	KEY_SIZE_384,
	(bit_32[]){0x8CB91E82, 0xA3386D28, 0x0F5D6F7E, 0x50E641DF, 0x152F7109, 0xED5456B4, 0x12B1DA19, 0x7FB71123, 0xACD3A729, 0x901D1A71, 0x87470013, 0x3107EC53}, // p
	(bit_32[]){0x7BC382C6, 0x3D8C150C, 0x3C72080A, 0xCE05AFA0, 0xC2BEA28E, 0x4FB22787, 0x139165EF, 0xBA91F90F, 0x8AA5814A, 0x503AD4EB, 0x04A8C7DD, 0x22CE2826}, // a
	(bit_32[]){0x04A8C7DD, 0x22CE2826, 0x8B39B554, 0x16F0447C, 0x2FB77DE1, 0x07DCD2A6, 0x2E880EA5, 0x3EEB62D5, 0x7CB43902, 0x95DBC994, 0x3AB78696, 0xFA504C11}, // b
	(bit_32[]){0x1D1C64F0, 0x68CF45FF, 0xA2A63A81, 0xB7C13F6B, 0x8847A3E7, 0x7EF14FE3, 0xDB7FCAFE, 0x0CBD10E8, 0xE826E034, 0x36D646AA, 0xEF87B2E2, 0x47D4AF1E}, // x
	(bit_32[]){0x8ABE1D75, 0x20F9C2A4, 0x5CB1EB8E, 0x95CFD552, 0x62B70B29, 0xFEEC5864, 0xE19C054F, 0xF9912928, 0x0E464621, 0x77918111, 0x42820341, 0x263C5315}, // y
	(bit_32[]){0x8CB91E82, 0xA3386D28, 0x0F5D6F7E, 0x50E641DF, 0x152F7109, 0xED5456B3, 0x1F166E6C, 0xAC0425A7, 0xCF3AB6AF, 0x6B7FC310, 0x3B883202, 0xE9046565}, // n
	0x1}; // h
	
static const DOMAIN_PARAMETER brainpoolP512r1 = {
	KEY_SIZE_512,
	(bit_32[]){0xAADD9DB8, 0xDBE9C48B, 0x3FD4E6AE, 0x33C9FC07, 0xCB308DB3, 0xB3C9D20E, 0xD6639CCA, 0x70330871, 0x7D4D9B00, 0x9BC66842, 0xAECDA12A, 0xE6A380E6, 0x2881FF2F, 0x2D82C685, 0x28AA6056, 0x583A48F3}, // p
	(bit_32[]){0x7830A331, 0x8B603B89, 0xE2327145, 0xAC234CC5, 0x94CBDD8D, 0x3DF91610, 0xA83441CA, 0xEA9863BC, 0x2DED5D5A, 0xA8253AA1, 0x0A2EF1C9, 0x8B9AC8B5, 0x7F1117A7, 0x2BF2C7B9, 0xE7C1AC4D, 0x77FC94CA}, // a
	(bit_32[]){0x3DF91610, 0xA83441CA, 0xEA9863BC, 0x2DED5D5A, 0xA8253AA1, 0x0A2EF1C9, 0x8B9AC8B5, 0x7F1117A7, 0x2BF2C7B9, 0xE7C1AC4D, 0x77FC94CA, 0xDC083E67, 0x984050B7, 0x5EBAE5DD, 0x2809BD63, 0x8016F723}, // b
	(bit_32[]){0x81AEE4BD, 0xD82ED964, 0x5A21322E, 0x9C4C6A93, 0x85ED9F70, 0xB5D916C1, 0xB43B62EE, 0xF4D0098E, 0xFF3B1F78, 0xE2D0D48D, 0x50D1687B, 0x93B97D5F, 0x7C6D5047, 0x406A5E68, 0x8B352209, 0xBCB9F822}, // x
	(bit_32[]){0x7DDE385D, 0x566332EC, 0xC0EABFA9, 0xCF7822FD, 0xF209F700, 0x24A57B1A, 0xA000C55B, 0x881F8111, 0xB2DCDE49, 0x4A5F485E, 0x5BCA4BD8, 0x8A2763AE, 0xD1CA2B2F, 0xA8F05406, 0x78CD1E0F, 0x3AD80892}, // y
	(bit_32[]){0xAADD9DB8, 0xDBE9C48B, 0x3FD4E6AE, 0x33C9FC07, 0xCB308DB3, 0xB3C9D20E, 0xD6639CCA, 0x70330870, 0x553E5C41, 0x4CA92619, 0x41866119, 0x7FAC1047, 0x1DB1D381, 0x085DDADD, 0xB5879682, 0x9CA90069}, // n
	0x1}; // h
	
// ##################################################################
// http://www.secg.org/sec2-v2.pdf

static const DOMAIN_PARAMETER secp128r1 = {
	KEY_SIZE_128,
	(bit_32[]){0xFFFFFFFD, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},	// p
	(bit_32[]){0xFFFFFFFD, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFC},	// a
	(bit_32[]){0xE87579C1, 0x1079F43D, 0xD824993C, 0x2CEE5ED3},	// b
	(bit_32[]){0x161FF752, 0x8B899B2D, 0x0C28607C, 0xA52C5B86},	// x
	(bit_32[]){0x27B6916A, 0x894D3AEE, 0x7106FE80, 0x5FC34B44},	// y	
	(bit_32[]){0x3FFFFFFF, 0x7FFFFFFF, 0xBE002472, 0x0613B5A3},	// n
	0x04};

static const DOMAIN_PARAMETER secp192k1 = {
	KEY_SIZE_192,																		// key_size
	(bit_32[]){0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFE, 0xFFFFEE37},	// p
	(bit_32[]){0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000},	// a
	(bit_32[]){0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000003},	// b
	(bit_32[]){0xDB4FF109, 0xC057E9DF, 0x26B07CF3, 0x80b7f3d1, 0x1DA5D1B3, 0xEAE06C7D},	// x
	(bit_32[]){0x9B2F2F6D, 0x9C5628A7, 0x844163D0, 0x15BE8634, 0x4082AA88, 0xD95E2F9D},	// y
	(bit_32[]){0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFE, 0x26F2FC17, 0x0F69466A, 0x74DEFD8D},	// n 
	0x1};																						// h
	
static const DOMAIN_PARAMETER secp256k1 = {
	KEY_SIZE_256,
	(bit_32[]){0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFE, 0xFFFFFC2F}, // p
	(bit_32[]){0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000}, // a
	(bit_32[]){0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000007}, // b
	(bit_32[]){0x79BE667E, 0xF9DCBBAC, 0x55A06295, 0xCE870B07, 0x029BFCDB, 0x2DCE28D9, 0x59F2815B, 0x16F81798}, // x
	(bit_32[]){0x483ADA77, 0x26A3C465, 0x5DA4FBFC, 0x0E1108A8, 0xFD17B448, 0xA6855419, 0x9C47D08F, 0xFB10D4B8}, // y
	(bit_32[]){0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFE, 0xBAAEDCE6, 0xAF48A03B, 0xBFD25E8C, 0xD0364141}, // n
	0x1}; // h
	
static const DOMAIN_PARAMETER secp256r1 = {
	KEY_SIZE_256,
	(bit_32[]){0xFFFFFFFF, 0x00000001, 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF}, // p
	(bit_32[]){0xFFFFFFFF, 0x00000001, 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFC}, // a
	(bit_32[]){0x5AC635D8, 0xAA3A93E7, 0xB3EBBD55, 0x769886BC, 0x651D06B0, 0xCC53B0F6, 0x3BCE3C3E, 0x27D2604B}, // b
	(bit_32[]){0x6B17D1F2, 0xE12C4247, 0xF8BCE6E5, 0x63A440F2, 0x77037D81, 0x2DEB33A0, 0xF4A13945, 0xD898C296}, // x
	(bit_32[]){0x4FE342E2, 0xFE1A7F9B, 0x8EE7EB4A, 0x7C0F9E16, 0x2BCE3357, 0x6B315ECE, 0xCBB64068, 0x37BF51F5}, // y
	(bit_32[]){0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xBCE6FAAD, 0xA7179E84, 0xF3B9CAC2, 0xFC632551}, // n
	0x1}; // h

// #####################################################################
// Edward curve domain parameter
// http://cr.yp.to/elligator/elligator-20130527.pdf
// https://eprint.iacr.org/2013/647.pdf
static const DOMAIN_PARAMETER curve1174 = {
	KEY_SIZE_256,
	(bit_32[]){0x7FFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFF7}, // p
	(bit_32[]){0, 0, 0, 0, 0, 0, 0, 0}, // a
	(bit_32[]){0, 0, 0, 0, 0, 0, 0, 0}, // b
	(bit_32[]){0, 0, 0, 0, 0, 0, 0, 0}, // x
	(bit_32[]){0, 0, 0, 0, 0, 0, 0, 0x00000004}, // y
	(bit_32[]){0x1ffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xf77965c4, 0xdfd30734, 0x8944d45f, 0xd166c971}, // n
	0,	// h
	(bit_32[]){0x7FFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFb61}}; // d
	
// https://tools.ietf.org/html/draft-josefsson-eddsa-ed25519-02	
static const DOMAIN_PARAMETER ed25519 = {
	KEY_SIZE_256,
	(bit_32[]){0x7FFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFED}, // p
	(bit_32[]){0, 0, 0, 0, 0, 0, 0, 0}, // a
	(bit_32[]){0, 0, 0, 0, 0, 0, 0, 0}, // b
	(bit_32[]){0x216936D3, 0xCD6E53FE, 0xC0A4E231, 0xFDD6DC5C, 0x692CC760, 0x9525A7B2, 0xC9562D60, 0x8F25D51A}, // x
	(bit_32[]){0x66666666, 0x66666666, 0x66666666, 0x66666666, 0x66666666, 0x66666666, 0x66666666, 0x66666658}, // y
	(bit_32[]){0x10000000, 0x00000000, 0x00000000, 0x00000000, 0x14DEF9DE, 0xA2F79CD6, 0x5812631A, 0x5CF5D3ED}, // n
	0,	// h
	(bit_32[]){0x52036CEE, 0x2B6FFE73, 0x8CC74079, 0x7779E898, 0x00700A4D, 0x4141D8AB, 0x75EB4DCA, 0x135978A3}}; // d

	
// compare array with array	
// check if a >= p
// return: a is 0=lesser, 1=equal, 2=grather
bit_16 compareAA(bit_32 *a, bit_32 *b, bit_16 len){
	bit_16 i;
	for (i = 0; i < len; i++){
		if (a[i] > b[i]){
			return 2;
		}
		if (a[i] < b[i]){
			return 0;
		}
	}
	return 1;	
}

// compare array wiht integer
// check if a >= p
// return:  0=lesser, 1=equal, 2=grather
bit_16 compareAI(bit_32 *a, bit_32 b, bit_16 len){

	bit_16 i;	
	if (a[len-1]>b){
		return 2;
	}
	else {
		for (i = len-2; i >=0; i--){
			if (a[i] > 0){
				return 2;
			}
		}
		if (a[len-1]==b){
			return 1;
		}
		else{
			return 0;
		}
	}	
}

// compare two points in affine plane
// p < q
// return 0=not equal, 1=equal 
bit_16 compareP(point *p, point *q){
	short i;
	for (i=0;i<keySize;i++){
		if ((p->x[i] != q->x[i]) || (p->y[i] != q->y[i])){
			return 0;
		}
	}
	return 1;
}
	
// copy array tmp to array c
void cpArray(bit_32 *tmp, bit_32 *c, bit_16 len){
	/*bit_16 i;
	for (i = 0; i<len;i++){
		c[i] = tmp[i]; 
	}
	*/
	while (len--) * c++ = * tmp++;
}

// copy a point p to q
void cpPoint(point *p, point *q){

	bit_16 i;
	for (i = 0; i< keySize;i++){
		q->x[i] = p->x[i];
		q->y[i] = p->y[i];
	}
}

// copy a point p to q
void cpPointZ(point *p, point *q){

	bit_16 i;
	for (i = 0; i< keySize;i++){
		q->x[i] = p->x[i];
		q->y[i] = p->y[i];
		q->z[i] = p->z[i];
	}
}

// get the position from MSB
bit_16 getMSBPos(bit_32 *a, bit_16 len){

	bit_16 i,j;
	bit_16 it = BIT_SIZE*len;
	bit_16 bit;
	for(i=0;i<len;i++){		
		for(j=BIT_SIZE-1;j>=0;j--){
			it--;
			if ((a[i]>>j)&0x1){
				return it;
			}
		}
	}
	return it;
}	

// get the ith bit position
bit_16 getIthBit(bit_32 *a,bit_16 pos){

	bit_16 idx1, idx2;
	idx1 = pos/BIT_SIZE;
	idx2 = pos%BIT_SIZE;
	return (a[keySizeM1-idx1]>>idx2)&0x1;
}

// get the ith bit position for random scalar
bit_16 getIthBit_rs(bit_32 *a, bit_16 pos){

	bit_16 idx1, idx2;
	idx1 = pos/BIT_SIZE;
	idx2 = pos%BIT_SIZE;
	return (a[keySize2-1-idx1]>>idx2)&0x1;
}

// get the w -ith bits position
bit_16 getIthBits(bit_32 *a,bit_16 lsbPos,bit_16 w){

	bit_16 idx1, idx2, ret, tmp, var;
	
	idx1 = lsbPos/BIT_SIZE;
	idx2 = lsbPos%BIT_SIZE;
	var = keySizeM1-idx1;
	if ((tmp=idx2+w)>BIT_SIZE){
		tmp=tmp-BIT_SIZE;
		ret = (a[var-1]&0x3)<<2 |(a[var]>>idx2);
	}
	else {
		ret = (a[var]>>idx2)&~(0xffffffff<<w);
	}
	
	return ret;
}


// Multiplication with 2
// shift array right MSB...LSB
void shiftArrayLeft_mod(bit_32 *a, bit_16 len){

	bit_16 carry1 = 0;
	bit_16 carry2 = 0;
	bit_16 i;
		
	for (i=len-1;i>=0;i--){
		carry1 = (bit_16)(a[i]>>31)&0x1;
		a[i] = (a[i] << 1) | (carry2 & 0x1);
		carry2 = carry1;
	}
	if (carry2){
		sl_tmp1[keySizeM1] = carry2;
		cpArray(a,&sl_tmp1[keySize],keySize);
		omp_div(sl_tmp1,keySize2,dp->p,keySize,sl_tmp2,a);
	}
	else if(compareAA(a,dp->p,keySize)==2){
		cpArray(a,sl_tmp1,keySize);
		gmp_sub_mod(sl_tmp1,dp->p,a);
	}

}

// Multiplication with 2 for fast reduction 
// shift array right MSB...LSB
void shiftArrayLeft_fast_reduction(bit_32 *a, bit_16 len, bit_16 *carry){

	bit_16 carry1 = 0;
	bit_16 carry2 = 0;
	bit_16 i;
		
	for (i=len-1;i>=0;i--){
		carry1 = (bit_16)(a[i]>>31)&0x1;
		a[i] = (a[i] << 1) | (carry2 & 0x1);
		carry2 = carry1;
	}
	*carry=carry2;

}

// Multiplication with 2 with barrett_reduction
// shift array right MSB...LSB
void shiftArrayLeft_br_mod(bit_32 *a, bit_16 len){

	bit_32 carry1 = 0;
	bit_32 carry2 = 0;
	bit_16 i;
	
	for (i=len-1;i>=0;i--){
		carry1 = (bit_16)(a[i]>>31)&0x1;
		a[i] = (a[i] << 1) | (carry2 & 0x1);
		carry2 = carry1;
		sl_tmp1[i] = 0;
	}
	if (carry2){
		sl_tmp1[keySizeM1] = carry2;
		cpArray(a,&sl_tmp1[keySize],keySize);
		br_mod(sl_tmp1,a);
		
	}
	else if(compareAA(a,dp->p,keySize)==2){
		cpArray(a,sl_tmp1,keySize);
		gmp_sub_mod(sl_tmp1,dp->p,a);
	}
}

// Multiplication with 2 with fast reduction
// shift array right MSB...LSB
void shiftArrayLeft_fr_mod(bit_32 *a, bit_16 len){

	bit_32 carry1 = 0;
	bit_32 carry2 = 0;
	bit_16 i;
	
	for (i=len-1;i>=0;i--){
		carry1 = (bit_16)(a[i]>>31)&0x1;
		a[i] = (a[i] << 1) | (carry2 & 0x1);
		carry2 = carry1;
		sl_tmp1[i] = 0;
	}
	if (carry2){
		sl_tmp1[keySizeM1] = carry2;
		cpArray(a,&sl_tmp1[keySize],keySize);
		fast_reduction_p256(sl_tmp1,a);
		
	}
	else if(compareAA(a,dp->p,keySize)==2){
		cpArray(a,sl_tmp1,keySize);
		gmp_sub_mod(sl_tmp1,dp->p,a);
	}
}

// Division with two
// shift array right MSB...LSB
void shiftArrayRight(bit_32 *a, bit_16 len){

	bit_16 carry1 = 0;
	bit_16 carry2 = 0;
	bit_16 i;
	
	for (i=0;i<len;i++){
		carry1 = (bit_16)a[i]&0x1;
		a[i] = (a[i]>>1) | carry2<<31;
		carry2 = carry1;
	}
}

// shift array by steps to right MSB...LSB
void shiftArrayRightBySteps(bit_32 *a, bit_16 len, bit_16 steps){

	bit_16 i, idx1,idx2;
	idx1 = steps/BIT_SIZE;
	idx2 = steps%BIT_SIZE;
	
	for (i=len-1;i>=idx1;i--){
		if ((i-idx1) >=0){
			a[i] = (a[i-idx1] >> (idx2));
		}
		if (((i-idx1-1) >=0) && (idx2 != 0)){
			a[i] = a[i] | (a[i-idx1-1]<<(BIT_SIZE-idx2));
		}
	}
	for (i=idx1-1;i>=0;i--){
		a[i] = 0;
	}
}

// shift array by steps to left MSB...LSB
void shiftArrayLeftBySteps(bit_32 *a, bit_16 len, bit_16 steps){

	bit_16 i, idx1,idx2;
	idx1 = steps/BIT_SIZE;
	idx2 = steps%BIT_SIZE;
	
	for (i=0;i<len-idx1;i++){
		if ((i+idx1) < len){
			a[i] = (a[i+idx1] << (idx2));
		}
		if (((i+idx1+1) < len) && (idx2 != 0)){
			a[i] = a[i] | (a[i+idx1+1]>>(BIT_SIZE-idx2));
		}
	}
	for (i=len-idx1;i<len;i++){
		a[i] = 0;
	}
}

// AND operation on an array
void andArray(bit_32 *a, bit_32 *b, bit_32 *c, bit_16 len){
	bit_16 i;
	for (i = 0;i<len;i++){
		c[i] = a[i] & b[i];
	}
}

// random permutation of memory for the Points
void rndPermPoint(bit_32 *rndNum){
	bit_16 tmp;
	
	// random permutation by knuth
	tmp = rp[0];
	rp[0] = rp[rndNum[0]&0xf];
	rp[rndNum[0]&0xf] = tmp;
	tmp = rp[1];
	rp[1] = rp[(rndNum[0]>>4)&0xf];
	rp[(rndNum[0]>>4)&0xf] = tmp;
	tmp = rp[2];
	rp[2] = rp[(rndNum[0]>>8)&0xf];
	rp[(rndNum[0]>>8)&0xf] = tmp;
	tmp = rp[3];
	rp[3] = rp[(rndNum[0]>>12)&0xf];
	rp[(rndNum[0]>>12)&0xf] = tmp;
	tmp = rp[4];
	rp[4] = rp[(rndNum[0]>>16)&0xf];
	rp[(rndNum[0]>>16)&0xf] = tmp;
	tmp = rp[5];
	rp[5] = rp[(rndNum[0]>>20)&0xf];
	rp[(rndNum[0]>>20)&0xf] = tmp;
	tmp = rp[6];
	rp[6] = rp[(rndNum[0]>>24)&0xf];
	rp[(rndNum[0]>>24)&0xf] = tmp;
	tmp = rp[7];
	rp[7] = rp[(rndNum[0]>>28)&0xf];
	rp[(rndNum[0]>>28)&0xf] = tmp;
	
	tmp = rp[8];
	rp[8] = rp[rndNum[1]&0xf];
	rp[rndNum[1]&0xf] = tmp;
	tmp = rp[9];
	rp[9] = rp[(rndNum[1]>>4)&0xf];
	rp[(rndNum[1]>>4)&0xf] = tmp;
	tmp = rp[10];
	rp[10] = rp[(rndNum[1]>>8)&0xf];
	rp[(rndNum[1]>>8)&0xf] = tmp;
	tmp = rp[11];
	rp[11] = rp[(rndNum[1]>>12)&0xf];
	rp[(rndNum[1]>>12)&0xf] = tmp;
	tmp = rp[12];
	rp[12] = rp[(rndNum[1]>>16)&0xf];
	rp[(rndNum[1]>>16)&0xf] = tmp;
	tmp = rp[13];
	rp[13] = rp[(rndNum[1]>>20)&0xf];
	rp[(rndNum[1]>>20)&0xf] = tmp;
	tmp = rp[14];
	rp[14] = rp[(rndNum[1]>>24)&0xf];
	rp[(rndNum[1]>>24)&0xf] = tmp;
	tmp = rp[15];
	rp[15] = rp[(rndNum[1]>>28)&0xf];
	rp[(rndNum[1]>>28)&0xf] = tmp;

	for (tmp = 0; tmp<16;tmp++){
		cpPointZ(&oArrPoint[tmp],&arrPoint[rp[tmp]]);
	}
}

void free_memory(){
	// add_mod
	free(add_tmp1);
	free(add_tmp2);
	free(add_tmp3);
	// sub_mod
	free(sub_tmp);
	// mul_mod
	free(mul_tmp1);
	free(mul_tmp2);	
	// fast reduction
	free(s1);
	free(s2);
	free(s3);
	free(s4);
	free(s5);
	free(s6);
	free(s7);
	free(s8);
	free(s9);
	// barrett reduction
	// br_mod
	free(br_q);
	free(br_tmp1111);
	free(br_tmp11);
	free(br_tmp22);
	free(br_tmp33);
	free(br_tmp44);
	// point_add
	free(p_s);
	free(p_tmp1);
	free(p_tmp2);
	free(p_tmp3);
	free(p_tmp4);
	free(p_tmp5);
	// shift left mod
	free(sl_tmp1);
	free(sl_tmp2);
	// montgomery inversion
	free(mi_u);
	free(mi_v);
	free(mi_x1);
	free(mi_x2);
	free(mi_tmp);
	// point projective 
	free(p_u);
	free(p_uu);
	free(p_v);
	free(p_R);
	free(p_A); 
	
	// konstanten
	free(c1);
	free(c3);
	free(cp2);
	
	// sliding window
	for (i = 0; i < 16; ++i){
		free(oArrPoint[i].x);
		free(oArrPoint[i].y);
		free(oArrPoint[i].z);
		free(arrPoint[i].x);
		free(arrPoint[i].y);
		free(arrPoint[i].z);
	}
	free(oArrPoint);
	free(arrPoint);
	
	// edward curve in projective plan
	free(p_B);
	free(p_C);
	free(p_D);
	free(p_E);
	free(p_F);
	free(p_G);
}

void init_memory(){
	bit_16 i;

	keySize   = dp->key_size;
	keySizeM1 = keySize-1;
	keySizeP1 = keySize+1;
	keySize2  = keySize*2;
	keySize4  = keySize*4;
	
	pMSBpos = getMSBPos(dp->p,keySize)+1;

	// add_mod
	add_tmp1 = calloc(keySize2,sizeof(bit_32));
	add_tmp2 = calloc(keySize2,sizeof(bit_32));
	add_tmp3 = calloc(keySize2,sizeof(bit_32));

	// sub_mod
	sub_tmp = calloc(keySize,sizeof(bit_32));
	
	// mul_mod
	mul_tmp1 = calloc(keySize2,sizeof(bit_32));
	mul_tmp2 = calloc(keySize2,sizeof(bit_32));
	
	// fast reduktion
	s1 = calloc(keySize,sizeof(bit_32));
	s2 = calloc(keySize,sizeof(bit_32));
	s3 = calloc(keySize,sizeof(bit_32));
	s4 = calloc(keySize,sizeof(bit_32));
	s5 = calloc(keySize,sizeof(bit_32));
	s6 = calloc(keySize,sizeof(bit_32));
	s7 = calloc(keySize,sizeof(bit_32));
	s8 = calloc(keySize,sizeof(bit_32));
	s9 = calloc(keySize,sizeof(bit_32));
	
	// barrett reduction
	// br_mod
	br_q = calloc(keySize2,sizeof(bit_32));
	br_tmp1111 = calloc(keySize4,sizeof(bit_32));
	br_tmp11 = calloc(keySize2,sizeof(bit_32));
	br_tmp22 = calloc(keySize2,sizeof(bit_32));
	br_tmp33 = calloc(keySize2,sizeof(bit_32));
	br_tmp44 = calloc(keySize2,sizeof(bit_32));
	
	// point 
	p_s = calloc(keySize,sizeof(bit_32));
	p_tmp1 = calloc(keySize,sizeof(bit_32));
	p_tmp2 = calloc(keySize,sizeof(bit_32));
	p_tmp3 = calloc(keySize,sizeof(bit_32));
	p_tmp4 = calloc(keySize,sizeof(bit_32));
	p_tmp5 = calloc(keySize,sizeof(bit_32));
	
	// shift left mod
	sl_tmp1 = calloc(keySize2,sizeof(bit_32));
	sl_tmp2 = calloc(keySize2,sizeof(bit_32));
	
	// montgomery inversion
	mi_u = calloc(dp->key_size,sizeof(bit_32));
	mi_v = calloc(dp->key_size,sizeof(bit_32));
	mi_x1 = calloc(dp->key_size,sizeof(bit_32));
	mi_x2 = calloc(dp->key_size,sizeof(bit_32));
	mi_tmp = calloc(dp->key_size*2,sizeof(bit_32));
	
	// point projective
	p_u = calloc(dp->key_size,sizeof(bit_32));
	p_uu = calloc(dp->key_size,sizeof(bit_32));
	p_v = calloc(dp->key_size,sizeof(bit_32));
	p_R = calloc(dp->key_size,sizeof(bit_32));
	p_A = calloc(dp->key_size,sizeof(bit_32));
	
	// konstanten
	c1 = calloc(keySize,sizeof(bit_32));
	c1[keySizeM1] = 0x1;
	c3 = calloc(keySize,sizeof(bit_32));
	c3[keySizeM1] = 0x3;
	cp2 = calloc(keySize2,sizeof(bit_32));
	cpArray(dp->p,&cp2[keySize],keySize);

	// array with points for sliding window
	oArrPoint = (point *) calloc(16,sizeof(point));
	arrPoint = (point *) calloc(16,sizeof(point));
	for (i = 0; i < 16; ++i){
		oArrPoint[i].x = calloc(dp->key_size,sizeof(bit_32));
		oArrPoint[i].y = calloc(dp->key_size,sizeof(bit_32));
		oArrPoint[i].z = calloc(dp->key_size,sizeof(bit_32));
		arrPoint[i].x = calloc(dp->key_size,sizeof(bit_32));
		arrPoint[i].y = calloc(dp->key_size,sizeof(bit_32));
		arrPoint[i].z = calloc(dp->key_size,sizeof(bit_32));
	}
	
	// Edward curve in projective plane
	
	p_B = calloc(keySize,sizeof(bit_32));
	p_C = calloc(keySize,sizeof(bit_32));
	p_D = calloc(keySize,sizeof(bit_32));
	p_E = calloc(keySize,sizeof(bit_32));
	p_F = calloc(keySize,sizeof(bit_32));
	p_G = calloc(keySize,sizeof(bit_32));
}
	
// initial elliptic curve domain parameter
bit_16 init_ecc(bit_32 name){
	switch (name){
		case SECP256K1 :
			dp = &secp256k1;
			break;
		case SECP256R1 :
			dp = &secp256r1;
			break;
		case SECP192K1 : 
			dp = &secp192k1;
			break;
			
		case P_192 :
			dp = &p_192;
			break;
		case P_224 :
			dp = &p_224;
			break;
		case P_256 :
			dp = &p_256;
			break;
		case P_384 :
			dp = &p_384;
			break;
		case P_521 :
			dp = &p_521;
			break;
			
		case BRAINPOOLP256R1 :
			dp = &brainpoolP256r1;
			break;
		case BRAINPOOLP320R1 :
			dp = &brainpoolP320r1;
			break;
		case BRAINPOOLP384R1 :
			dp = &brainpoolP384r1;
			break;
		case BRAINPOOLP512R1 :
			dp = &brainpoolP512r1;
			break;
			
		case CURVE1174 :
			dp = &curve1174;
			break;
		case ED25519 :
			dp = &ed25519;
			break;
	}
//	free_memory();
//	conclude_barrett_reduction();
	init_memory();
	init_barrett_reduction();
	return keySize;
}

// Guide to Elliptic Curve Cryptography
// Algorithm 2.5 Multiprecision addition
// carry,c = a+b
// carry,a = a+b
void gmp_add(bit_32 *a, bit_32 *b, bit_32 *c, bit_16 *carry,bit_16 len){
    
	bit_16 i;
	bit_32 tmp;
	*carry = 0;
	for (i = len-1; i >= 0 ;i--){
		
		tmp = a[i]+b[i]+*carry;
		if (tmp < a[i] || ((b[i]==0xffffffff) && (*carry==1))){
			*carry = 1;
		}
		else {
			*carry = 0;
		}
		c[i] = tmp;
	}
}

// Guide to Elliptic Curve Cryptography
// Algorithm 2.6 Multiprecision subtraction
// carry,c = a-b 
// carry,a = a-b
void gmp_sub(bit_32 *a, bit_32 *b, bit_32 *c, bit_16 *carry, bit_16 len){

    bit_16 i;
	bit_32 tmp;
	*carry = 0;
	for (i = len-1; i >= 0 ;i--){
		
		tmp = a[i]-b[i]+*carry;
		if (tmp > a[i] || ((b[i]==0xffffffff) && (*carry==-1))){
			*carry = -1;
		}
		else {
			*carry = 0;
		}
		c[i] = tmp;
	}
}

// Guide to Elliptic Curve Cryptography
// Algorithm 2.7 Addition in Fp
// c = a+b mod p
void gmp_add_mod(bit_32 *a, bit_32 *b, bit_32 *c){

	bit_16 carry;
	gmp_add(a,b,c,&carry, keySize);
	
	if (carry == 1){
		gmp_sub(c,dp->p,c,&carry, keySize);
	}
	else if (compareAA(c,dp->p,keySize)==2){
		gmp_sub(c,dp->p,c,&carry, keySize);
	}

}

// Guide to Elliptic Curve Cryptography
// Algorithm 2.8 Subtraction in Fp
// c = a-b mod p
void gmp_sub_mod(bit_32 *a, bit_32 *b, bit_32 *c){

	bit_16 carry;
/*	bit_32 *tmp; sub_tmp
	tmp = calloc(dp->key_size,sizeof(bit_32));*/
	gmp_sub(a,b,c,&carry, keySize);

	if (carry==-1){
		gmp_add(c,dp->p,c,&carry, keySize);
	}
	else if (compareAA(c,dp->p,keySize)==2){
		gmp_add(c,dp->p,c,&carry, keySize);
	}
//	else{
//		cpArray(sub_tmp,c,keySize);
//	}
//	free(tmp);
}


// Guide to Elliptic Curve Cryptography
// Algorithm 2.9 Integer multiplication (operand scanning form)
// carry,c = a*b 
void gmp_mul(bit_32 *a, bit_32 *b, bit_32 *c, bit_16 len){

	// cortex m3 multiply-accumulate instructions for 32-bit or 64-bit results

    bit_16 i,j;
	bit_64 uv = 0;
	bit_64 ta=0, tb=0, tc=0;
	bit_16 t = (len-1);
	bit_16 l = (len*2);
	
	for (i=0;i<l;i++){
		c[i] = 0x0;
	}
	
	for (i = t; i >= 0 ;i--){
		uv = uv & 0x00000000ffffffff;
		for (j = t; j >=0;j--){
			ta = (bit_64)a[i];
			tb = (bit_64)b[j];
			tc = (bit_64)c[i+j+1];
			uv = tc+ta*tb+(uv>>32);
			c[i+j+1]=(bit_32)(uv&0x00000000ffffffff);			
		}
		c[i] =  (bit_32)(uv>>32);
	}
	
}


// Guide to Elliptic Curve Cryptography
// Algorithm 2.13 Integer squaring
// carry,c = a^2 
void gmp_sqrt(bit_32 *a, bit_32 *c, bit_16 len){

	// cortex m3 multiply-accumulate instructions for 32-bit or 64-bit results
	
    bit_16 i,j,k;
	bit_16 carry=0;
	bit_32 r0=0,r1=0,r2=0,rt=0;
	bit_64 uv = 0;
	bit_64 ai=0, aj=0;
	bit_32 t = (len-1);
	bit_32 t2 = (len*2)-2;
	
	
	for (k=t2; k>=0;k--){
		for (j=t; j>=0; j--){
		
			i = k-j;
			if (i<=j && i >= 0){
			
				ai = (bit_64)a[i];
				aj = (bit_64)a[j];
				uv = ai*aj;
				if (i<j){
					// carry bit 
					if (uv>>63){
						r2 = r2 +1;
					}
					uv = uv<<1;
				}
				rt = r0;
				r0 = r0+(bit_32)(uv&0x00000000ffffffff);
				if (rt>r0){
					carry = 1;
				}
				else {
					carry = 0;
				}
				
				rt = r1;
				r1 = r1+(bit_32)(uv>>32)+carry;
				if (rt>r1){
					carry = 1;
				}
				else {
					carry = 0;
				}
				r2 = r2 + carry;
			}
		
		}
		c[k+1] = r0;
		r0 = r1;
		r1 = r2;
		r2 = 0;
		
	}
	c[0] = r0;	
}


// mp division by own
// determine q and r
// a = q * b + r
void omp_div(bit_32 *a, bit_16 lenA, bit_32 *b, bit_16 lenB, bit_32 *c, bit_32 *r){
	bit_16 t,n,tn;
	bit_16 tv1, tv2;
	bit_16 i;
	bit_16 carry;
	bit_32 *base;
	bit_32 *tmpc;
	bit_32 *tmp1;
	bit_32 *tmp11;
	bit_32 *tmp2222;
	bit_32 *x;
	bit_32 *y;
	tmpc = calloc(lenA,sizeof(bit_32));
	base = calloc(lenA,sizeof(bit_32));
	tmp1 = calloc(lenB,sizeof(bit_32));
	tmp11 = calloc(lenA,sizeof(bit_32));
	tmp2222 = calloc(lenA*2,sizeof(bit_32));
	x = calloc(lenA*2,sizeof(bit_32));
	y = calloc(lenA*2,sizeof(bit_32));

	if (compareAI(b,0x1,lenB) == 1){
		cpArray(a,c,lenA);
		cpArray(tmpc,r,lenB);
	}
	else if(compareAI(b,0x2,lenB) == 1){
		carry = (bit_16)a[lenA-1] & 0x1;
		cpArray(a,c,lenA);
		shiftArrayRight(c, lenA);
		r[dp->key_size-1] = carry;
	}
	else{
	
		cpArray(b,&tmp11[lenA-lenB],lenB);
		cpArray(a,&x[lenA],lenA);
			
		t=getMSBPos(a,lenA);
		n=getMSBPos(b,lenB);
		tn=t-n;
		// base^tn
		tv1 = tn/BIT_SIZE;
		tv2 = tn%BIT_SIZE;
		base[lenA-1-tv1] = 0x1<<tv2;
		gmp_mul(tmp11,base,y,lenA);
		for (i=0; i<tn+1;i++){
			if(compareAA(x,y,lenB*4)>=1 ){
				gmp_sub(x,y,tmp2222,&carry, lenA*2);
				cpArray(tmp2222,x,lenA*2);
				gmp_add(tmpc,base,tmp11,&carry,lenA);
				cpArray(tmp11,tmpc,lenA);				
			}
			shiftArrayRight(y, lenA*2);
			shiftArrayRight(base, lenA);
		}

		cpArray(tmpc,c,lenA);
		cpArray(&x[lenA*2-lenB],r,lenB);	
	}	
	free(tmp1);
	free(tmp11);
	free(tmp2222);
	free(x);
	free(y);
	free(base);
	free(tmpc);
}

// c = (a * b) % p
void gmp_mul_mod(bit_32 *a, bit_32 *b, bit_32 *c){

	bit_16 comp;
	
	gmp_mul(a,b,mul_tmp1,keySize);
	comp = compareAA(mul_tmp1,cp2,keySize2);
	// tmp1 > tmp2 (tmp2=p)
	if(comp==2){
		br_mod(mul_tmp1,c);
	}
	// tmp1 == tmp2
	else if (comp==1){
		c[keySize]=0;
	}
	else{
		cpArray(&mul_tmp1[keySize],c,keySize);
	}
}

// multiplication modulo fast reduction p256
// c = (a * b) % p
void gmp_mul_mod_fr(bit_32 *a, bit_32 *b, bit_32 *c){

	bit_16 comp;
	
	gmp_mul(a,b,mul_tmp1,keySize);
	comp = compareAA(mul_tmp1,cp2,keySize2);
	// tmp1 > tmp2 (tmp2=p)
	if(comp==2){
		//br_mod(mul_tmp1,c);
		fast_reduction_p256(mul_tmp1,c);
	}
	// tmp1 == tmp2
	else if (comp==1){
		c[keySize]=0;
	}
	else{
		cpArray(&mul_tmp1[keySize],c,keySize);
	}
}


// Guide to Elliptic Curve Cryptography
// Algorithm 2.20 Inversion in Fp using the extended Euclidean algorithm
// a^-1 mop p 
void gmp_inv(bit_32 *a, bit_32 *c){

	bit_32 *u;
	bit_32 *v;
	bit_32 *x;
	bit_32 *x1;
	bit_32 *x2;
	bit_32 *q;
	bit_32 *r;
	bit_32 *tmp1;
	u = calloc(dp->key_size,sizeof(bit_32));
	v = calloc(dp->key_size,sizeof(bit_32));
	x = calloc(dp->key_size,sizeof(bit_32));
	x1 = calloc(dp->key_size,sizeof(bit_32));
	x2 = calloc(dp->key_size,sizeof(bit_32));
	q = calloc(dp->key_size,sizeof(bit_32));
	r = calloc(dp->key_size,sizeof(bit_32));
	tmp1 = calloc(dp->key_size,sizeof(bit_32));
	
	// step 1
	cpArray(a,u,dp->key_size);
	cpArray(dp->p,v,dp->key_size);

	// step 2
	x1[dp->key_size-1] = 0x1;
	
	// step 3
	while (compareAI(u,0x1,dp->key_size) != 1){
				
		// step 3.1
		omp_div(v,dp->key_size,u,dp->key_size,q,r);
		gmp_mul_mod(q,x1,tmp1);
		gmp_sub_mod(x2,tmp1,x);

		// step 3.2
		// evt adresse übergeben
		cpArray(u,v,dp->key_size);
		cpArray(r,u,dp->key_size);
		cpArray(x1,x2,dp->key_size);
		cpArray(x,x1,dp->key_size);
	}
			
	cpArray(x1,c,dp->key_size);
	
	free(u);
	free(v);
	free(x);
	free(x1);
	free(x2);
	free(q);
	free(r);
	free(tmp1);

}


// Handbook of Elliptic and Hyperelliptic Curve Cryptography
// Algorithm 11.9 Prime field inversion
// 
void hmp_inv(bit_32 *a, bit_32 *c){
	short carry;
	short i;
	bit_32 *q;
	bit_32 *u;
	bit_32 *z;
	bit_32 *tmp1;
	bit_32 *tmp11;
	q = calloc(dp->key_size,sizeof(bit_32));
	u = calloc(dp->key_size,sizeof(bit_32));
	z = calloc(dp->key_size,sizeof(bit_32));
	tmp1 = calloc(dp->key_size,sizeof(bit_32));
	tmp11 = calloc(dp->key_size*2,sizeof(bit_32));
	cpArray(a,z,dp->key_size);

	u[dp->key_size-1] = 0x1;

	while (compareAI(z,0x1,dp->key_size) != 1){
			omp_div(dp->p, dp->key_size, z, dp->key_size, q, tmp11);
			gmp_mul(q, z, tmp11, dp->key_size);
			gmp_sub(dp->p, &tmp11[dp->key_size], z, &carry, dp->key_size);
			gmp_sub(dp->p,q,tmp1,&carry,dp->key_size);
			gmp_mul_mod(tmp1, u, u);
	}
	cpArray(u,c,dp->key_size);
	
	free(q);
	free(u);
	free(z);
	free(tmp1);
	free(tmp11);
}


/// Guide to Elliptic Curve Cryptography
// Partial Montgomery inversion in Fp
// a^-1 mop p 
void gmp_mont_inv(bit_32 *a, bit_32 *c){
	short carry,i;
	bit_32 k = 0;
	////////////////////////////////////////////
	// phase I
	
	for (i=0;i<keySize;i++){
		mi_u[i]  = 0;
		mi_v[i]  = 0;
		mi_x1[i] = 0;
		mi_x2[i] = 0;
	}
	
	// step 1
	cpArray(a,mi_u,keySize);
	cpArray(dp->p,mi_v,keySize);

	mi_x1[keySizeM1] = 1;
	
	// step 2
	while (compareAI(mi_v,0x0,keySize) == 2){
		
		// step 2.1
		if ((mi_v[keySizeM1]&0x1)==0){
			shiftArrayRight(mi_v, keySize);
			shiftArrayLeft_mod(mi_x1, keySize);
		}
		else if((mi_u[keySizeM1]&0x1)==0){
			shiftArrayRight(mi_u, keySize);
			shiftArrayLeft_mod(mi_x2, keySize);
		}

		else if(compareAA(mi_v, mi_u, keySize)>=1){
			gmp_sub(mi_v, mi_u, mi_tmp, &carry, keySize);
			shiftArrayRight(mi_tmp, keySize);
			cpArray(mi_tmp,mi_v,keySize);
			gmp_add(mi_x2, mi_x1, mi_tmp, &carry, keySize);
			cpArray(mi_tmp,mi_x2,keySize);
			shiftArrayLeft_mod(mi_x1, keySize);
		}
		else{
			gmp_sub(mi_u, mi_v, mi_tmp, &carry, keySize);
			shiftArrayRight(mi_tmp, keySize);
			cpArray(mi_tmp,mi_u,keySize);
			gmp_add(mi_x2, mi_x1, mi_tmp, &carry, keySize);
			cpArray(mi_tmp,mi_x1,keySize);
			shiftArrayLeft_mod(mi_x2, keySize);
		}
		
		// step 2.2
		 k=k+1;
	}	
	
	// step 3
	if (compareAI(mi_u,0x1,keySize) != 1){
		//return not invertible
		printf("not invertible\n");
	}
	// step 4
	if (compareAA(mi_x1, dp->p, keySize)==2){
		gmp_sub(mi_x1, dp->p, mi_tmp, carry, keySize);
		cpArray(mi_tmp,mi_x1,keySize);
	}

	
	////////////////////////////////////////////
	// phase II
	// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.75.8377&rep=rep1&type=pdf
	mi_tmp[keySizeM1]=0;
	while (k!=0){
		k--;

		if ((mi_x1[keySizeM1]&0x1)==0){
			shiftArrayRight(mi_x1, keySize);
		}
		else {

			gmp_add(mi_x1, dp->p, &mi_tmp[keySize], &carry, keySize);
			if (carry==1){
				mi_tmp[keySizeM1]=0x1;
			}
			shiftArrayRight(mi_tmp, keySize2);
			cpArray(&mi_tmp[keySize],mi_x1,keySize);
		}
		
	}
	cpArray(mi_x1,c,keySize);
	
}



void conclude_barrett_reduction(){
	free(bkp1);
	free(bkp1M1);
	free(mikro);
}

// Software Implementation of the NIST Elliptic Curves Over Prime Fields
// Algorithm 2.29 Fast reduction modulo
void fast_reduction_p256(bit_32 *a,bit_32 *c){

	bit_16 carry, carry1;
	
	s1[0] = a[8];
	s1[1] = a[9];
	s1[2] = a[10];
	s1[3] = a[11];
	s1[4] = a[12];
	s1[5] = a[13];
	s1[6] = a[14];
	s1[7] = a[15];
	
	s2[0] = a[0];
	s2[1] = a[1];
	s2[2] = a[2];
	s2[3] = a[3];
	s2[4] = a[4];
	s2[5] = 0;
	s2[6] = 0;
	s2[7] = 0;
	
	s3[0] = 0;
	s3[1] = a[0];
	s3[2] = a[1];
	s3[3] = a[2];
	s3[4] = a[3];
	s3[5] = 0;
	s3[6] = 0;
	s3[7] = 0;

	s4[0] = a[0];
	s4[1] = a[1];
	s4[2] = 0;
	s4[3] = 0;
	s4[4] = 0;
	s4[5] = a[5];
	s4[6] = a[6];
	s4[7] = a[7];
	
	s5[0] = a[7];
	s5[1] = a[2];
	s5[2] = a[0];
	s5[3] = a[1];
	s5[4] = a[2];
	s5[5] = a[4];
	s5[6] = a[5];
	s5[7] = a[6];
	
	s6[0] = a[5];
	s6[1] = a[7];
	s6[2] = 0;
	s6[3] = 0;
	s6[4] = 0;
	s6[5] = a[2];
	s6[6] = a[3];
	s6[7] = a[4];
	
	s7[0] = a[4];
	s7[1] = a[6];
	s7[2] = 0;
	s7[3] = 0;
	s7[4] = a[0];
	s7[5] = a[1];
	s7[6] = a[2];
	s7[7] = a[3];
	
	s8[0] = a[3];
	s8[1] = 0;
	s8[2] = a[5];
	s8[3] = a[6];
	s8[4] = a[7];
	s8[5] = a[0];
	s8[6] = a[1];
	s8[7] = a[2];
	
	s9[0] = a[2];
	s9[1] = 0;
	s9[2] = a[4];
	s9[3] = a[5];
	s9[4] = a[6];
	s9[5] = 0;
	s9[6] = a[0];
	s9[7] = a[1];
	
	// s1+2s2+2s3+s4+s5-s6-s7-s8-s9 mod p256
	
	shiftArrayLeft_fast_reduction(s2,keySize, &carry1); // 2*s2
	shiftArrayLeft_fast_reduction(s3,keySize, &carry); // 2*s3
	carry1 +=carry;
	gmp_add(s1,s2,c,&carry,keySize);
	carry1 +=carry;
	gmp_add(c,s3,c,&carry,keySize);
	carry1 +=carry;
	gmp_add(c,s4,c,&carry,keySize);
	carry1 +=carry;
	gmp_add(c,s5,c,&carry,keySize);
	carry1 +=carry;
	gmp_sub(c,s6,c,&carry,keySize);
	carry1 +=carry;
	gmp_sub(c,s7,c,&carry,keySize);
	carry1 +=carry;
	gmp_sub(c,s8,c,&carry,keySize);
	carry1 +=carry;
	gmp_sub(c,s9,c,&carry,keySize);
	carry1 +=carry;
	
	for (i=carry1;i<0;i++){
		gmp_add(c,dp->p,c,&carry,keySize);
	}
	for (i=carry1;i>0;i--){
		gmp_sub(c,dp->p,c,&carry,keySize);
	}
	
}

void init_barrett_reduction(){
	//bit_16 b = 4;
	short carry,i;
	bit_16 k;
	bit_32 *tb2k;
	bit_32 *tmp11;
	bit_32 *tmp22;
	bit_32 *tmp1111;
	tb2k  = calloc(keySize4,sizeof(bit_32));
	tmp11 = calloc(keySize2,sizeof(bit_32));
	tmp22 = calloc(keySize2,sizeof(bit_32));
	tmp1111 = calloc(keySize4,sizeof(bit_32));
	
	tb2k[keySize4-1]=1;

	bkp1  = calloc(keySize2,sizeof(bit_32));
	bkp1M1= calloc(keySize2,sizeof(bit_32));
	mikro = calloc(keySize2,sizeof(bit_32));	
	
	k=getMSBPos(dp->p,keySize)/2+1;
	br_k = k*2;
	// µ=[b^2k/p]
	shiftArrayLeftBySteps(tb2k,keySize4,(k+1)*2);
	cpArray(&tb2k[keySize2],bkp1,keySize2);	
	for (i=0;i<keySize2;i++){
		br_tmp22[i]=0;
	}
	br_tmp22[keySize2-1] = 1;
	gmp_sub(bkp1, br_tmp22, bkp1M1, &carry, keySize2);
	
	shiftArrayLeftBySteps(tb2k,keySize4,(k-1)*2);
	cpArray(dp->p,&tmp11[keySize],keySize);
	omp_div(tb2k, keySize4, tmp11, keySize2, tmp1111, tmp22);
	cpArray(&tmp1111[keySize2],mikro,keySize2);	
	
	free(tb2k);
	free(tmp11);
	free(tmp22);
	free(tmp1111);
	
}

void br_mod(bit_32 *a, bit_32 *c){
	
	bit_16 carry;	
	// step 1
	cpArray(a,br_tmp11,keySize2);
	shiftArrayRightBySteps(br_tmp11, keySize2, br_k-2);
	gmp_mul(br_tmp11, mikro, br_tmp1111, keySize2);
	shiftArrayRightBySteps(br_tmp1111, keySize4, br_k+2);
	cpArray(&br_tmp1111[keySize2],br_q,keySize2);
	
	//step 2
	andArray(a, bkp1M1, br_tmp22, keySize2); // mod b^k+1
	gmp_mul(&br_q[keySize], dp->p, br_tmp33, keySize);
	andArray(br_tmp33, bkp1M1, br_q, keySize2);
	gmp_sub(br_tmp22, br_q, br_tmp33, &carry, keySize2);
	
	// step 3
	if (carry == -1){
		gmp_add(br_tmp33,bkp1,br_tmp33,&carry,keySize2);
	}	
	
	// step 4
	while (compareAA(br_tmp33, cp2, keySize2)>=1){
		gmp_sub(br_tmp33,cp2,br_tmp33,&carry,keySize2);
	}
	
	cpArray(&br_tmp33[keySize],c,keySize);
}



// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.140.7944&rep=rep1&type=pdf
// modular division 
void mod_div(bit_32 *a, bit_32 *b, bit_32 *c){
	short i,carry;
	bit_16 com;
	
	cpArray(a,mi_u,keySize);		// u =U
	cpArray(b,mi_x1,keySize); 		// x1=A
	cpArray(dp->p,mi_x2,keySize); 	// x2=B
	for (i=0;i<keySize;i++){		// c =V
		c[i] = 0;					
	}
	while ((com=compareAA(mi_x1,mi_x2,keySize))!=1){
	if((mi_x1[keySizeM1]&0x1) == 0){
			shiftArrayRight(mi_x1, keySize);
			if ((mi_u[keySizeM1]&0x1) == 0){
				shiftArrayRight(mi_u, keySize);
			}
			else {
				gmp_add(mi_u,dp->p,mi_u,&carry,keySize);
				shiftArrayRight(mi_u, keySize);
				mi_u[0] = mi_u[0] | (carry<<31);
			}
		}
		else if((mi_x2[keySizeM1]&0x1) == 0){
			shiftArrayRight(mi_x2, keySize);
			if ((c[keySizeM1]&0x1) == 0){
				shiftArrayRight(c, keySize);
			}
			else {
				gmp_add(c,dp->p,c,&carry,keySize);
				shiftArrayRight(c, keySize);
				c[0] = c[0] | (carry<<31);
			}
		}
		else if (com == 2){
			gmp_sub(mi_x1,mi_x2,mi_x1,&carry,keySize);
			shiftArrayRight(mi_x1, keySize);
			
			gmp_sub_mod(mi_u,c,mi_u);
			
			if ((mi_u[keySizeM1]&0x1) == 0){
				shiftArrayRight(mi_u, keySize);
			}
			else {
				gmp_add(mi_u,dp->p,mi_u,&carry,keySize);
				shiftArrayRight(mi_u, keySize);
				mi_u[0] = mi_u[0] | (carry<<31);
			}
		}
		else {
			gmp_sub(mi_x2,mi_x1,mi_x2,&carry,keySize);
			shiftArrayRight(mi_x2, keySize);
			gmp_sub_mod(c,mi_u,c);
	
			if ((c[keySizeM1]&0x1) == 0){
				shiftArrayRight(c, keySize);
			}
			else {
				gmp_add(c,dp->p,c,&carry,keySize);
				shiftArrayRight(c, keySize);
				c[0] = c[0] | (carry<<31);
			}
		}
	}
	
}


// #############################################################################################
// #############################################################################################
// Weierstrass affine plane binary method

void p_add_aff_w(point *p, point *q){
	// s = Yp -Yq / Xp - Xq
	// Xr = s^2 - Xp - Xq
	// Yr = -Yp + s(Xp - Xr)
	// p=q
	if (compareP(p,q)==1){
		p_dbl_aff_w(p);
		return;
	}
	// p = -q
	else if((compareAA(p->x,q->x,keySize)==1) && (compareAA(p->y,q->y,keySize)!=1)){
		// P addiert mit inverse = Neutralesellement
		for (i=0;i<keySize;i++){
			p->x[i] = 0;
			p->y[i] = 0;
		}
		return;
	}
	// p + 0 = p
	if ((compareAI(p->x,0,keySize)==1) && (compareAI(p->y,0,keySize)==1)){
		cpPoint(q,p);
		return;
	}
	if ((compareAI(q->x,0,keySize)==1) && (compareAI(q->y,0,keySize)==1)){
		return;
	}
	// s = Yp -Yq / Xp - Xq
	gmp_sub_mod(p->y,q->y,p_tmp1);
	gmp_sub_mod(p->x,q->x,p_tmp2);
	mod_div(p_tmp1,p_tmp2,p_s);
	// Xr = s^2 - Xp - Xq
	gmp_mul_mod(p_s,p_s,p_tmp1);
	gmp_sub_mod(p_tmp1,p->x,p_tmp2);
	gmp_sub_mod(p_tmp2,q->x,p_tmp3);
	// Yr = -Yp + s(Xp - Xr)
	gmp_sub_mod(p->x,p_tmp3,p_tmp1);
	gmp_mul_mod(p_s,p_tmp1,p_tmp2);
	gmp_sub_mod(p_tmp2,p->y,p->y);
    cpArray(p_tmp3, p->x, keySize);
	
}

void p_dbl_aff_w(point *p){
	// s = (3*px^2 + a) / (2*py)
	// rx = s^2 - 2px
	// ry = -py + s(px - rx)
	if (compareAI(p->y, 0, keySize)==1){
		for (i=0;i<keySize;i++){
			p->x[i] = 0;
		}
		return;
	}
	// s = (3*px^2 + a) / (2*py)
	gmp_mul_mod(p->x,p->x,p_tmp2);
	gmp_mul_mod(c3,p_tmp2,p_tmp3);
	gmp_add_mod(p_tmp3,dp->a,p_tmp1);
	cpArray(p->y,p_tmp2,keySize);
	shiftArrayLeft_br_mod(p_tmp2,keySize);
	mod_div(p_tmp1,p_tmp2,p_s);
	
	// rx = s^2 - 2px
	gmp_mul_mod(p_s,p_s,p_tmp1);
	cpArray(p->x,p_tmp2,keySize);
	shiftArrayLeft_br_mod(p_tmp2,keySize);
	gmp_sub_mod(p_tmp1,p_tmp2,p_tmp3);
	
	// Yr = -Yp + s(Xp - Xr)
	gmp_sub_mod(p->x,p_tmp3,p_tmp1);
	gmp_mul_mod(p_s,p_tmp1,p_tmp2);
	gmp_sub_mod(p_tmp2,p->y,p_tmp1);
	
	cpArray(p_tmp3, p->x, keySize);
	cpArray(p_tmp1, p->y, keySize);
}

char p_mul_aff_w_bin(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 msbPos;
	point *baseP;
	
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	cpPoint(baseP,public_key);
	
	for (i=msbPos;i>=0;i--){
		
		p_dbl_aff_w(public_key);
		if (getIthBit(private_key,i)){
			p_add_aff_w(public_key,baseP);
		}
	}
	free(baseP);
	return 0;
}


// Weierstrass affine plane binary method
// ###############################################################################################
// ###############################################################################################

// ==============================================================================================================================

// ###############################################################################################
// #############################################################################################
// Weierstrass projective plane Brainpool256

// http://hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-mmadd-1998-cmo
// 
// add-2007-bl
void p_add_proj_w(point *p, point *q, point *r){
	
	// addition with infinity
	if ((compareAI(p->x,0,keySize)==1) && (compareAI(p->y,1,keySize)==1) && (compareAI(p->z,0,keySize)==1)){
		cpPointZ(q,r);
		return;
	}
	if ((compareAI(q->x,0,keySize)==1) && (compareAI(q->y,1,keySize)==1) && (compareAI(q->z,0,keySize)==1)){
		cpPointZ(p,r);
		return;
	}
	gmp_mul_mod(p->x,q->z,p_A);				// A = X1*Z2	U1
	gmp_mul_mod(q->x,p->z,p_B);				// B = X2*Z1	U2
	gmp_mul_mod(p->y,q->z,p_C);				// C = Y1*Z2	S1
	gmp_mul_mod(q->y,p->z,p_D);				// D = Y2*Z1	S2
	if (compareAA(p_A, p_B, keySize)==1){
		if (compareAA(p_C, p_D, keySize)!=1){
			for (i=0;i<keySize;i++){
				r->x[i] = 0;
				r->y[i] = 0;
				r->z[i] = 0;
			}
			r->y[keySizeM1] = 0x1;
			return;
		}
		else{
			cpPointZ(p,r);
			p_dbl_proj_w(r);
			return;
		}
	}
	gmp_mul_mod(p->z,q->z,p_E);				// E = Z1*Z2	ZZ
	gmp_add_mod(p_A,p_B,p_F);				// T=U1+U2			
	gmp_mul_mod(p_F,p_F,p_G);				// TT = T^2 (F^2)	
	gmp_add_mod(p_C,p_D,p_s);				// M = C+D (S1+S2) 
	gmp_mul_mod(p_A,p_B,p_tmp1);			// U1*U2 
	gmp_mul_mod(p_E,p_E,p_tmp2);			// E^2 (ZZ^2)
	gmp_mul_mod(dp->a,p_tmp2,p_tmp3);		// a*ZZ^2
	gmp_sub_mod(p_G,p_tmp1,p_tmp1);			//TT-U1*U2
	gmp_add_mod(p_tmp1,p_tmp3,p_tmp1);		// R = G-(A*B)+a*E^2 (TT-U1*U2+a*ZZ^2)
	gmp_mul_mod(p_E,p_s,p_tmp2);			// F = ZZ*M
	gmp_mul_mod(p_tmp2,p_s,p_D);			// L = F*M
	gmp_mul_mod(p_D,p_D,p_tmp3);			// LL = L^2
	gmp_add_mod(p_F,p_D,p_C);				// (T+L)
	gmp_mul_mod(p_C,p_C,p_D);				// (T+L)^2
	gmp_sub_mod(p_D,p_G,p_C);				// (T+L)^2 - TT
	gmp_sub_mod(p_C,p_tmp3,p_D);			// G = (T+L)^2 - TT - LL
	gmp_mul_mod(p_tmp1,p_tmp1,p_tmp4);		// R^2
	shiftArrayLeft_br_mod(p_tmp4,keySize);	// 2*R^2
	gmp_sub_mod(p_tmp4,p_D,p_tmp4);			// W = 2*R^2-G
	gmp_mul_mod(p_tmp2,p_tmp4,r->x);		// F*W
	shiftArrayLeft_br_mod(r->x,keySize);	// X = 2*F*W
	
	shiftArrayLeft_br_mod(p_tmp4,keySize);	// 2*W
	gmp_sub_mod(p_D,p_tmp4,p_tmp4);			// G-2*W
	gmp_mul_mod(p_tmp1,p_tmp4,p_s);			// R*(G-2*W)
	shiftArrayLeft_br_mod(p_tmp3,keySize);	// 2*LL
	gmp_sub_mod(p_s,p_tmp3,r->y);			// Y = R*(G-2*W) - 2*LL
	
	gmp_mul_mod(p_tmp2,p_tmp2,p_tmp3);		// F^2
	gmp_mul_mod(p_tmp3,p_tmp2,r->z);		// F^2*F
	shiftArrayLeft_br_mod(r->z,keySize);	// 2*F^2*F
	shiftArrayLeft_br_mod(r->z,keySize);	// 4*F^2*F
	
	
}

// Source: 2007 Bernstein–Lange.
// dbl-2007-bl
void p_dbl_proj_w(point *p){

	if (compareAI(p->y,0,keySize)==1){
		for (i=0;i<keySize;i++){
			p->x[i] = 0;
			p->y[i] = 0;
			p->z[i] = 0;
		}
		p->y[keySizeM1] = 0x1;
		return;
	}
	if ((compareAI(p->x,0,keySize)==1) && (compareAI(p->y,1,keySize)==1) && (compareAI(p->z,0,keySize)==1)){
		return;
	}
	gmp_mul_mod(p->x,p->x,p_A);				// XX = X1^2
	gmp_mul_mod(p->z,p->z,p_B);				// ZZ = Z1^2
	gmp_mul_mod(dp->a,p_B,p_tmp1);			// a*ZZ
	gmp_mul_mod(c3,p_A,p_tmp2);				// 3*XX
	gmp_add_mod(p_tmp1,p_tmp2,p_C);			// w = a*ZZ + 3*XX
	gmp_mul_mod(p->y,p->z,p_D);				// Y1*Z1
	shiftArrayLeft_br_mod(p_D,keySize);		// s = 2*Y1*Z1
	gmp_mul_mod(p_D,p_D,p_E);				// ss = s^2
	gmp_mul_mod(p_D,p_E,p->z);				// z = sss = s^2 * s
	gmp_mul_mod(p->y,p_D,p_G);				// R = Y1 * s
	gmp_mul_mod(p_G,p_G,p_tmp1);			// RR = R^2
	gmp_add_mod(p->x,p_G,p_tmp2);			// X1+R
	gmp_mul_mod(p_tmp2,p_tmp2,p_s);			// (X1+R)^2
	gmp_sub_mod(p_s,p_A,p_tmp2); 			// (X1+R)^2 - XX
	gmp_sub_mod(p_tmp2,p_tmp1,p_s); 		// B = (X1+R)^2 - XX - RR
	gmp_mul_mod(p_C,p_C,p_tmp3);			// w^2
	cpArray(p_s,p_tmp4,keySize);
	shiftArrayLeft_br_mod(p_tmp4,keySize);	// 2*B
	gmp_sub_mod(p_tmp3,p_tmp4,p_tmp3);		// h = w^2-2*B
	gmp_mul_mod(p_tmp3,p_D,p->x);			// x = h*s
	gmp_sub_mod(p_s,p_tmp3,p_tmp4);			// B-h
	gmp_mul_mod(p_C,p_tmp4,p_tmp2);			// w*(B-h)
	shiftArrayLeft_br_mod(p_tmp1,keySize);	// 2*RR
	gmp_sub_mod(p_tmp2,p_tmp1,p->y);		// w*(B-h) - 2*RR
	
	
}

void transforProjToAff(point *p){
	
	// if Z > 1
	if (compareAI(p->z,1,keySize)==2){
		mod_div(p->x,p->z,p->x);
		mod_div(p->y,p->z,p->y);	
		mod_div(p->z,p->z,p->z);
	}
}

char p_mul_proj_w_bin(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 msbPos;
	point *baseP;
	
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	
	cpPointZ(baseP,public_key);
	for (i=msbPos;i>=0;i--){
		p_dbl_proj_w(public_key);
		if (getIthBit(private_key,i)){
			p_add_proj_w(public_key,baseP,public_key);
		}
	}
	transforProjToAff(public_key);
	free(baseP);
	return 0;
}

// pre compute Points for sliding window
void preComputSW_w(){
	bit_16 i;
	point *baseP;
	
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	cpPointZ(baseP,&oArrPoint[1]);
	cpPointZ(baseP,&oArrPoint[15]);
	
	for (i=2;i<15;i++){
		p_add_proj_w(&oArrPoint[15],&oArrPoint[1],&oArrPoint[15]);
		cpPointZ(&oArrPoint[15],&oArrPoint[i]);
	}
	// point at infinity
	oArrPoint[0].y[keySizeM1] = 0x1;
	p_add_proj_w(&oArrPoint[15],&oArrPoint[1],&oArrPoint[15]);
}

// handbook of elliptic and hyperelliptic curve cryptography
// angepasste sliding window mehtode, im buch ist der rahmen nicht fest hier schon immer 4 oder 2 bits
char p_mul_proj_w_rndSW(bit_32 *private_key, point *public_key, bit_32 *rnd1, bit_32 *rnd2){
	
	bit_16 w, i, h; 
	signed int msbPos,msbPosRnd1;
	rndPermPoint(rnd1);
	cpPointZ(&arrPoint[rp[0]],public_key);
	msbPos = getMSBPos(private_key,keySize);
	msbPosRnd1 = getMSBPos(rnd2,keySize)-1;
	if (getIthBit(rnd2,msbPosRnd1) && (msbPos>=4)){
		w=4;
	}
	else {
		w=2;
	}
	msbPosRnd1--;
	msbPos = msbPos - (msbPos%w);
		
	while (msbPos >= 0){
		for (i=0;i<w;i++){
			p_dbl_proj_w(public_key);	
		}
		h = getIthBits(private_key,msbPos,w);
		p_add_proj_w(public_key,&arrPoint[rp[h]],public_key);
		if ((getIthBit(rnd2,msbPosRnd1)) && (msbPos>=4)){
			w=4;
		}
		else {
			w=2;
		}
		msbPosRnd1--;
		msbPos=msbPos-w;
	}
	transforProjToAff(public_key);
	return 0;
}

void randomize_scalar(bit_32 *private_key, bit_32 rnd, bit_32 *rs_private_key){

	// k' = k + q*n , k=secret key, q=random number (32bit), n=group order
	bit_32 *rnd_tmp;
	bit_16 carry;
	rnd_tmp = calloc(keySize2,sizeof(bit_32));
	rnd_tmp[keySize2-1]=rnd;
	gmp_mul(dp->n,&rnd_tmp[keySize],rs_private_key,keySize);
	cpArray(private_key, &rnd_tmp[keySize],keySize);
	gmp_add(rnd_tmp,rs_private_key,rs_private_key,&carry,keySize2);	
	free(rnd_tmp);
}

// random scalar
char p_mul_proj_w_rndSca(bit_32 *private_key, point *public_key, bit_32 rnd){
	
	bit_16 i;
	bit_16 msbPos;
	bit_32 *pk;
	pk = calloc(keySize2,sizeof(bit_32));
	point *baseP;
	
	randomize_scalar(private_key, rnd, pk);
	msbPos = getMSBPos(pk,keySize2)-1;
	
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;	
	cpPointZ(baseP,public_key);
	for (i=msbPos;i>=0;i--){
		p_dbl_proj_w(public_key);
		if (getIthBit_rs(pk,i)){
			p_add_proj_w(public_key,baseP,public_key);
		}	

	}
	transforProjToAff(public_key);
	free(baseP);
	free(pk);
	return 0;
}

// double and add always
char p_mul_proj_w_daa(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 it1 = 0, it2;
	bit_16 msbPos;
	point *baseP, *poiP;
	
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	
	poiP = calloc(2,sizeof(point));
	poiP[0].x = calloc(keySize,sizeof(bit_32));
	poiP[0].y = calloc(keySize,sizeof(bit_32));
	poiP[0].z = calloc(keySize,sizeof(bit_32));
			
	poiP[1].x = calloc(keySize,sizeof(bit_32));
	poiP[1].y = calloc(keySize,sizeof(bit_32));
	poiP[1].z = calloc(keySize,sizeof(bit_32));
	
	cpPointZ(baseP,&poiP[0]);
	cpPointZ(baseP,&poiP[1]);
	p_dbl_proj_w(&poiP[1]);
	for (i=msbPos;i>=0;i--){
		it2 = (1-getIthBit(private_key,i));
		p_dbl_proj_w(&poiP[it1]);
		p_add_proj_w(&poiP[it1],baseP,&poiP[it2]);
	}
	transforProjToAff(&poiP[0]);
	cpPointZ(&poiP[0],public_key);
	free(baseP);
	free(poiP);
	return 0;
}

// randomise coordiantes
// Resistance Against Differential Power Analysis For Elliptic Curve Cryptosystems
char p_mul_proj_w_rndCoor(bit_32 *private_key, point *public_key, bit_32 rnd){
	
	bit_16 i;
	bit_16 msbPos;
	point *baseP;
	bit_32 *rnd1;

	rnd1 = calloc(keySize,sizeof(bit_32));
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	rnd1[keySizeM1] = rnd;
	
	cpPointZ(baseP,public_key);
	for (i=msbPos;i>=0;i--){
		// randomise coordinate
		// fibonacci-LFSR x^32+x^23+x^7+1
		rnd1[keySizeM1] = ((rnd1[keySizeM1]>>31) ^ (rnd1[keySizeM1]>>22 & 1) ^ (rnd1[keySizeM1]>>6 & 1)) | (rnd1[keySizeM1]<<1 );
		gmp_mul_mod(public_key->x,rnd1,public_key->x);
		gmp_mul_mod(public_key->y,rnd1,public_key->y);
		gmp_mul_mod(public_key->z,rnd1,public_key->z);
		
		p_dbl_proj_w(public_key);
		if (getIthBit(private_key,i)){
			p_add_proj_w(public_key,baseP,public_key);
		}
	}
	transforProjToAff(public_key);
	free(baseP);
	free(rnd1);
	return 0;
}

// point multiplication montgomery ladder
char p_mul_proj_w_ml(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 it1 = 0, it2;
	bit_16 msbPos;
	point *baseP, *poiP;
	
	msbPos = getMSBPos(private_key,keySize);
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	
	poiP = calloc(2,sizeof(point));
	poiP[0].x = calloc(keySize,sizeof(bit_32));
	poiP[0].y = calloc(keySize,sizeof(bit_32));
	poiP[0].z = calloc(keySize,sizeof(bit_32));
			
	poiP[1].x = calloc(keySize,sizeof(bit_32));
	poiP[1].y = calloc(keySize,sizeof(bit_32));
	poiP[1].z = calloc(keySize,sizeof(bit_32));
	
	poiP[0].y[keySizeM1] = 0x1;
	cpPointZ(baseP,&poiP[1]);
	
	for (i=msbPos;i>=0;i--){
		it1 = getIthBit(private_key,i);
		it2 = 1-it1;
		p_add_proj_w(&poiP[0],&poiP[1],&poiP[it2]);
		p_dbl_proj_w(&poiP[it1]);
		
	}
	transforProjToAff(&poiP[0]);
	cpPointZ(&poiP[0],public_key);
	free(baseP);
	free(poiP);
	return 0;
}

// point multiplication randomised skalar splitting
char p_mul_proj_w_rss(bit_32 *private_key, point *public_key, bit_32 *rnd){
	
	bit_16 i;
	bit_16 msbPos1, msbPos2, carry;
	bit_32 *rnd2;
	rnd2 = calloc(keySize,sizeof(bit_32));
	point *baseP, *poiP;
	gmp_sub(private_key,rnd,rnd2,&carry,keySize);
	
	msbPos1 = getMSBPos(rnd,keySize)-1;
	msbPos2 = getMSBPos(rnd2,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	
	poiP = calloc(2,sizeof(point));
	poiP[0].x = calloc(keySize,sizeof(bit_32));
	poiP[0].y = calloc(keySize,sizeof(bit_32));
	poiP[0].z = calloc(keySize,sizeof(bit_32));
			
	poiP[1].x = calloc(keySize,sizeof(bit_32));
	poiP[1].y = calloc(keySize,sizeof(bit_32));
	poiP[1].z = calloc(keySize,sizeof(bit_32));
	
	cpPointZ(baseP,&poiP[0]);
	cpPointZ(baseP,&poiP[1]);
	
	for (i=msbPos1;i>=0;i--){
		p_dbl_proj_w(&poiP[0]);
		if (getIthBit(rnd,i)){
			p_add_proj_w(&poiP[0],baseP,&poiP[0]);
		}
	}
	
	for (i=msbPos2;i>=0;i--){
		p_dbl_proj_w(&poiP[1]);
		if (getIthBit(rnd2,i)){
			p_add_proj_w(&poiP[1],baseP,&poiP[1]);
		}
	}

	p_add_proj_w(&poiP[0],&poiP[1],public_key);
	transforProjToAff(public_key);
	free(baseP);
	free(poiP);
	return 0;
}


// Weierstrass projective plane brainpool256
// ###############################################################################################
// ###############################################################################################

// ==============================================================================================================================

// ###############################################################################################
// ############################################################################################
// projective plane in weierstrass with fast reduction p256


void p_add_aff_w_fr(point *p, point *q){
	// s = Yp -Yq / Xp - Xq
	// Xr = s^2 - Xp - Xq
	// Yr = -Yp + s(Xp - Xr)
	// p=q
	if (compareP(p,q)==1){
		p_dbl_aff_w_fr(p);
		return;
	}
	// p = -q
	else if((compareAA(p->x,q->x,keySize)==1) && (compareAA(p->y,q->y,keySize)!=1)){
		// P addiert mit inverse = Neutralesellement
		for (i=0;i<keySize;i++){
			p->x[i] = 0;
			p->y[i] = 0;
		}
		return;
	}
	// p + 0 = p
	if ((compareAI(p->x,0,keySize)==1) && (compareAI(p->y,0,keySize)==1)){
		cpPoint(q,p);
		return;
	}
	if ((compareAI(q->x,0,keySize)==1) && (compareAI(q->y,0,keySize)==1)){
		return;
	}
	// s = Yp -Yq / Xp - Xq
	gmp_sub_mod(p->y,q->y,p_tmp1);
	gmp_sub_mod(p->x,q->x,p_tmp2);
	mod_div(p_tmp1,p_tmp2,p_s);
	// Xr = s^2 - Xp - Xq
	gmp_mul_mod_fr(p_s,p_s,p_tmp1);
	gmp_sub_mod(p_tmp1,p->x,p_tmp2);
	gmp_sub_mod(p_tmp2,q->x,p_tmp3);
	// Yr = -Yp + s(Xp - Xr)
	gmp_sub_mod(p->x,p_tmp3,p_tmp1);
	gmp_mul_mod_fr(p_s,p_tmp1,p_tmp2);
	gmp_sub_mod(p_tmp2,p->y,p->y);
    cpArray(p_tmp3, p->x, keySize);
	
}

void p_dbl_aff_w_fr(point *p){
	// s = (3*px^2 + a) / (2*py)
	// rx = s^2 - 2px
	// ry = -py + s(px - rx)
	if (compareAI(p->y, 0, keySize)==1){
		for (i=0;i<keySize;i++){
			p->x[i] = 0;
		}
		return;
	}
	// s = (3*px^2 + a) / (2*py)
	gmp_mul_mod_fr(p->x,p->x,p_tmp2);
	gmp_mul_mod_fr(c3,p_tmp2,p_tmp3);
	gmp_add_mod(p_tmp3,dp->a,p_tmp1);
	cpArray(p->y,p_tmp2,keySize);
	shiftArrayLeft_fr_mod(p_tmp2,keySize);
	mod_div(p_tmp1,p_tmp2,p_s);
	
	// rx = s^2 - 2px
	gmp_mul_mod_fr(p_s,p_s,p_tmp1);
	cpArray(p->x,p_tmp2,keySize);
	shiftArrayLeft_fr_mod(p_tmp2,keySize);
	gmp_sub_mod(p_tmp1,p_tmp2,p_tmp3);
	
	// Yr = -Yp + s(Xp - Xr)
	gmp_sub_mod(p->x,p_tmp3,p_tmp1);
	gmp_mul_mod_fr(p_s,p_tmp1,p_tmp2);
	gmp_sub_mod(p_tmp2,p->y,p_tmp1);
	
	cpArray(p_tmp3, p->x, keySize);
	cpArray(p_tmp1, p->y, keySize);
}


// http://hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-mmadd-1998-cmo
// 
// add-2007-bl
void p_add_proj_w_fr(point *p, point *q, point *r){
	
	// addition with infinity
	if ((compareAI(p->x,0,keySize)==1) && (compareAI(p->y,1,keySize)==1) && (compareAI(p->z,0,keySize)==1)){
		cpPointZ(q,r);
		return;
	}
	if ((compareAI(q->x,0,keySize)==1) && (compareAI(q->y,1,keySize)==1) && (compareAI(q->z,0,keySize)==1)){
		cpPointZ(p,r);
		return;
	}
	gmp_mul_mod_fr(p->x,q->z,p_A);				// A = X1*Z2	U1
	gmp_mul_mod_fr(q->x,p->z,p_B);				// B = X2*Z1	U2
	gmp_mul_mod_fr(p->y,q->z,p_C);				// C = Y1*Z2	S1
	gmp_mul_mod_fr(q->y,p->z,p_D);				// D = Y2*Z1	S2
	if (compareAA(p_A, p_B, keySize)==1){
		if (compareAA(p_C, p_D, keySize)!=1){
			for (i=0;i<keySize;i++){
				r->x[i] = 0;
				r->y[i] = 0;
				r->z[i] = 0;
			}
			r->y[keySize-1] = 0x1;
			return;
		}
		else{
			cpPointZ(p,r);
			p_dbl_proj_w_fr(r);
			return;
		}
	}
	gmp_mul_mod_fr(p->z,q->z,p_E);				// E = Z1*Z2	ZZ
	gmp_add_mod(p_A,p_B,p_F);				// T=U1+U2			
	gmp_mul_mod_fr(p_F,p_F,p_G);				// TT = T^2 (F^2)	
	gmp_add_mod(p_C,p_D,p_s);				// M = C+D (S1+S2) 
	gmp_mul_mod_fr(p_A,p_B,p_tmp1);			// U1*U2 
	gmp_mul_mod_fr(p_E,p_E,p_tmp2);			// E^2 (ZZ^2)
	gmp_mul_mod_fr(dp->a,p_tmp2,p_tmp3);		// a*ZZ^2
	gmp_sub_mod(p_G,p_tmp1,p_tmp1);			//TT-U1*U2
	gmp_add_mod(p_tmp1,p_tmp3,p_tmp1);		// R = G-(A*B)+a*E^2 (TT-U1*U2+a*ZZ^2)
	gmp_mul_mod_fr(p_E,p_s,p_tmp2);			// F = ZZ*M
	gmp_mul_mod_fr(p_tmp2,p_s,p_D);			// L = F*M
	gmp_mul_mod_fr(p_D,p_D,p_tmp3);			// LL = L^2
	gmp_add_mod(p_F,p_D,p_C);				// (T+L)
	gmp_mul_mod_fr(p_C,p_C,p_D);				// (T+L)^2
	gmp_sub_mod(p_D,p_G,p_C);				// (T+L)^2 - TT
	gmp_sub_mod(p_C,p_tmp3,p_D);			// G = (T+L)^2 - TT - LL
	gmp_mul_mod_fr(p_tmp1,p_tmp1,p_tmp4);		// R^2
	shiftArrayLeft_fr_mod(p_tmp4,keySize);	// 2*R^2
	gmp_sub_mod(p_tmp4,p_D,p_tmp4);			// W = 2*R^2-G
	gmp_mul_mod_fr(p_tmp2,p_tmp4,r->x);		// F*W
	shiftArrayLeft_fr_mod(r->x,keySize);	// X = 2*F*W
	
	shiftArrayLeft_fr_mod(p_tmp4,keySize);	// 2*W
	gmp_sub_mod(p_D,p_tmp4,p_tmp4);			// G-2*W
	gmp_mul_mod_fr(p_tmp1,p_tmp4,p_s	);		// R*(G-2*W)
	shiftArrayLeft_fr_mod(p_tmp3,keySize);	// 2*LL
	gmp_sub_mod(p_s,p_tmp3,r->y);			// Y = R*(G-2*W) - 2*LL
	
	gmp_mul_mod_fr(p_tmp2,p_tmp2,p_tmp3);		// F^2
	gmp_mul_mod_fr(p_tmp3,p_tmp2,r->z);		// F^2*F
	shiftArrayLeft_fr_mod(r->z,keySize);	// 2*F^2*F
	shiftArrayLeft_fr_mod(r->z,keySize);	// 4*F^2*F
}

// Source: 2007 Bernstein–Lange.
// dbl-2007-bl
void p_dbl_proj_w_fr(point *p){

	if (compareAI(p->y,0,keySize)==1){
		for (i=0;i<keySize;i++){
			p->x[i] = 0;
			p->y[i] = 0;
			p->z[i] = 0;
		}
		p->y[keySize-1] = 0x1;
		return;
	}
	if ((compareAI(p->x,0,keySize)==1) && (compareAI(p->y,1,keySize)==1) && (compareAI(p->z,0,keySize)==1)){
		return;
	}
	gmp_mul_mod_fr(p->x,p->x,p_A);				// XX = X1^2
	gmp_mul_mod_fr(p->z,p->z,p_B);				// ZZ = Z1^2
	gmp_mul_mod_fr(dp->a,p_B,p_tmp1);			// a*ZZ
	gmp_mul_mod_fr(c3,p_A,p_tmp2);				// 3*XX
	gmp_add_mod(p_tmp1,p_tmp2,p_C);			// w = a*ZZ + 3*XX
	gmp_mul_mod_fr(p->y,p->z,p_D);				// Y1*Z1
	shiftArrayLeft_fr_mod(p_D,keySize);		// s = 2*Y1*Z1
	gmp_mul_mod_fr(p_D,p_D,p_E);				// ss = s^2
	gmp_mul_mod_fr(p_D,p_E,p->z);				// z = sss = s^2 * s
	gmp_mul_mod_fr(p->y,p_D,p_G);				// R = Y1 * s
	gmp_mul_mod_fr(p_G,p_G,p_tmp1);			// RR = R^2
	gmp_add_mod(p->x,p_G,p_tmp2);			// X1+R
	gmp_mul_mod_fr(p_tmp2,p_tmp2,p_s);			// (X1+R)^2
	gmp_sub_mod(p_s,p_A,p_tmp2); 			// (X1+R)^2 - XX
	gmp_sub_mod(p_tmp2,p_tmp1,p_s); 		// B = (X1+R)^2 - XX - RR
	gmp_mul_mod_fr(p_C,p_C,p_tmp3);			// w^2
	cpArray(p_s,p_tmp4,keySize);
	shiftArrayLeft_fr_mod(p_tmp4,keySize);	// 2*B
	gmp_sub_mod(p_tmp3,p_tmp4,p_tmp3);		// h = w^2-2*B
	gmp_mul_mod_fr(p_tmp3,p_D,p->x);			// x = h*s
	gmp_sub_mod(p_s,p_tmp3,p_tmp4);			// B-h
	gmp_mul_mod_fr(p_C,p_tmp4,p_tmp2);			// w*(B-h)
	shiftArrayLeft_fr_mod(p_tmp1,keySize);	// 2*RR
	gmp_sub_mod(p_tmp2,p_tmp1,p->y);		// w*(B-h) - 2*RR
}

void transforProjToAff_fr(point *p){
	
	// if Z > 1
	if (compareAI(p->z,1,keySize)==2){
		mod_div(p->x,p->z,p->x);
		mod_div(p->y,p->z,p->y);	
		mod_div(p->z,p->z,p->z);
	}
}

char p_mul_aff_w_bin_fr(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 msbPos;
	point *baseP;
	
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	cpPoint(baseP,public_key);
	
	for (i=msbPos;i>=0;i--){
		
		p_dbl_aff_w_fr(public_key);
		if (getIthBit(private_key,i)){
			p_add_aff_w_fr(public_key,baseP);
		}
	}
	free(baseP);
	return 0;
}

char p_mul_proj_w_bin_fr(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 msbPos;
	point *baseP;
	
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	
	cpPointZ(baseP,public_key);
	for (i=msbPos;i>=0;i--){
		p_dbl_proj_w_fr(public_key);
		if (getIthBit(private_key,i)){
			p_add_proj_w_fr(public_key,baseP,public_key);
		}
	}
	transforProjToAff_fr(public_key);
	free(baseP);
	return 0;
}

// handbook of elliptic and hyperelliptic curve cryptography
// angepasste sliding window mehtode, im buch ist der rahmen nicht fest hier schon immer 4 oder 2 bits
char p_mul_proj_w_rndSW_fr(bit_32 *private_key, point *public_key, bit_32 *rnd1, bit_32 *rnd2){
	
	bit_16 w, i, h; 
	signed int msbPos,msbPosRnd1;
	rndPermPoint(rnd1);
	cpPointZ(&arrPoint[rp[0]],public_key);
	msbPos = getMSBPos(private_key,keySize);
	msbPosRnd1 = getMSBPos(rnd2,keySize)-1;
	if (getIthBit(rnd2,msbPosRnd1) && (msbPos>=4)){
		w=4;
	}
	else {
		w=2;
	}
	msbPosRnd1--;
	msbPos = msbPos - (msbPos%w);
		
	while (msbPos >= 0){
		for (i=0;i<w;i++){
			p_dbl_proj_w_fr(public_key);	
		}
		h = getIthBits(private_key,msbPos,w);
		p_add_proj_w_fr(public_key,&arrPoint[rp[h]],public_key);
		if ((getIthBit(rnd2,msbPosRnd1)) && (msbPos>=4)){
			w=4;
		}
		else {
			w=2;
		}
		msbPosRnd1--;
		msbPos=msbPos-w;
	}
	transforProjToAff_fr(public_key);
	return 0;
}

// random scalar
char p_mul_proj_w_rndSca_fr(bit_32 *private_key, point *public_key, bit_32 rnd){
	
	bit_16 i;
	bit_16 msbPos;
	bit_32 *pk;
	pk = calloc(keySize2,sizeof(bit_32));
	point *baseP;
	
	randomize_scalar(private_key, rnd, pk);
	msbPos = getMSBPos(pk,keySize2)-1;
	
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;	
	cpPointZ(baseP,public_key);
	for (i=msbPos;i>=0;i--){
		p_dbl_proj_w_fr(public_key);
		if (getIthBit_rs(pk,i)){
			p_add_proj_w_fr(public_key,baseP,public_key);
		}	

	}
	transforProjToAff_fr(public_key);
	free(baseP);
	free(pk);
	return 0;
}

// double and add always
char p_mul_proj_w_daa_fr(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 it1 = 0, it2;
	bit_16 msbPos;
	point *baseP, *poiP;
	
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	
	poiP = calloc(2,sizeof(point));
	poiP[0].x = calloc(keySize,sizeof(bit_32));
	poiP[0].y = calloc(keySize,sizeof(bit_32));
	poiP[0].z = calloc(keySize,sizeof(bit_32));
			
	poiP[1].x = calloc(keySize,sizeof(bit_32));
	poiP[1].y = calloc(keySize,sizeof(bit_32));
	poiP[1].z = calloc(keySize,sizeof(bit_32));
	
	cpPointZ(baseP,&poiP[0]);
	cpPointZ(baseP,&poiP[1]);
	p_dbl_proj_w_fr(&poiP[1]);
	
	for (i=msbPos;i>=0;i--){
		it2 = 1-getIthBit(private_key,i);
		p_dbl_proj_w_fr(&poiP[it1]);
		p_add_proj_w_fr(&poiP[it1],baseP,&poiP[it2]);
	}
	transforProjToAff_fr(&poiP[0]);
	cpPointZ(&poiP[0],public_key);
	free(baseP);
	free(poiP);
	return 0;
}

// randomise coordiantes
// Resistance Against Differential Power Analysis For Elliptic Curve Cryptosystems
char p_mul_proj_w_rndCoor_fr(bit_32 *private_key, point *public_key, bit_32 rnd){
	
	bit_16 i;
	bit_16 msbPos;
	point *baseP;
	bit_32 *rnd1;

	rnd1 = calloc(keySize,sizeof(bit_32));
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	rnd1[keySizeM1] = rnd;
	
	cpPointZ(baseP,public_key);
	for (i=msbPos;i>=0;i--){
		// randomise coordinate
		// fibonacci-LFSR x^32+x^23+x^7+1
		rnd1[keySizeM1] = ((rnd1[keySizeM1]>>31) ^ (rnd1[keySizeM1]>>22 & 1) ^ (rnd1[keySizeM1]>>6 & 1)) | (rnd1[keySizeM1]<<1 );
		gmp_mul_mod_fr(public_key->x,rnd1,public_key->x);
		gmp_mul_mod_fr(public_key->y,rnd1,public_key->y);
		gmp_mul_mod_fr(public_key->z,rnd1,public_key->z);
		
		p_dbl_proj_w_fr(public_key);
		if (getIthBit(private_key,i)){
			p_add_proj_w_fr(public_key,baseP,public_key);
		}
	}
	transforProjToAff_fr(public_key);
	free(baseP);
	free(rnd1);
	return 0;
}

// point multiplication montgomery ladder
char p_mul_proj_w_ml_fr(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 it1 = 0, it2;
	bit_16 msbPos;
	point *baseP, *poiP;
	
	msbPos = getMSBPos(private_key,keySize);
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	
	poiP = calloc(2,sizeof(point));
	poiP[0].x = calloc(keySize,sizeof(bit_32));
	poiP[0].y = calloc(keySize,sizeof(bit_32));
	poiP[0].z = calloc(keySize,sizeof(bit_32));
			
	poiP[1].x = calloc(keySize,sizeof(bit_32));
	poiP[1].y = calloc(keySize,sizeof(bit_32));
	poiP[1].z = calloc(keySize,sizeof(bit_32));
	
	poiP[0].y[keySizeM1] = 0x1;
	cpPointZ(baseP,&poiP[1]);
	
	for (i=msbPos;i>=0;i--){
		it1 = getIthBit(private_key,i);
		it2 = 1-it1;
		p_add_proj_w_fr(&poiP[0],&poiP[1],&poiP[it2]);
		p_dbl_proj_w_fr(&poiP[it1]);
		
	}
	transforProjToAff_fr(&poiP[0]);
	cpPointZ(&poiP[0],public_key);
	free(baseP);
	free(poiP);
	return 0;
}

// point multiplication randomised skalar splitting
char p_mul_proj_w_rss_fr(bit_32 *private_key, point *public_key, bit_32 *rnd){
	
	bit_16 i;
	bit_16 msbPos1, msbPos2, carry;
	bit_32 *rnd2;
	rnd2 = calloc(keySize,sizeof(bit_32));
	point *baseP, *poiP;
	gmp_sub(private_key,rnd,rnd2,&carry,keySize);
	
	msbPos1 = getMSBPos(rnd,keySize)-1;
	msbPos2 = getMSBPos(rnd2,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	
	poiP = calloc(2,sizeof(point));
	poiP[0].x = calloc(keySize,sizeof(bit_32));
	poiP[0].y = calloc(keySize,sizeof(bit_32));
	poiP[0].z = calloc(keySize,sizeof(bit_32));
			
	poiP[1].x = calloc(keySize,sizeof(bit_32));
	poiP[1].y = calloc(keySize,sizeof(bit_32));
	poiP[1].z = calloc(keySize,sizeof(bit_32));
	
	cpPointZ(baseP,&poiP[0]);
	cpPointZ(baseP,&poiP[1]);
	
	for (i=msbPos1;i>=0;i--){
		p_dbl_proj_w_fr(&poiP[0]);
		if (getIthBit(rnd,i)){
			p_add_proj_w_fr(&poiP[0],baseP,&poiP[0]);
		}
	}
	
	for (i=msbPos2;i>=0;i--){
		p_dbl_proj_w_fr(&poiP[1]);
		if (getIthBit(rnd2,i)){
			p_add_proj_w_fr(&poiP[1],baseP,&poiP[1]);
		}
	}

	p_add_proj_w_fr(&poiP[0],&poiP[1],public_key);
	transforProjToAff_fr(public_key);
	free(baseP);
	free(poiP);
	return 0;
}


// projective plane in weierstrass with fast reduction p256
// #############################################################################################
// #############################################################################################

// ==============================================================================================================================

// ###############################################################################################
// #############################################################################################
// Edward Curve affine

// point addition in affine plane for edward curve
void p_add_aff_ed(point *p, point *q){
	// c=1
	// x=(pxqy+pyqx/(1+ px qx py qy))
	// y=(pyqy-pxqx/(1- px qx py qy))
	
	gmp_mul_mod(p->x,q->y,p_tmp1); // px * qy
	gmp_mul_mod(p->y,q->x,p_tmp2); // py * qx
	gmp_add_mod(p_tmp1,p_tmp2,p_tmp1); // px * qy + py * qx
	
	gmp_mul_mod(p->x,q->x,p_tmp2); // px * qx
	gmp_mul_mod(p->y,q->y,p_tmp3); // py * qy
	gmp_sub_mod(p_tmp2,p_tmp3,p_s); // px * qy - py * qx

	gmp_mul_mod(p_tmp2,p_tmp3,p_tmp4); // px*qx*py*qy
	gmp_mul_mod(p_tmp4,dp->d,p_tmp2); // d*px*qx*py*qy
	
	gmp_add_mod(c1,p_tmp2,p_tmp4); // 1 + d*px*qx*py*qy
	gmp_sub_mod(c1,p_tmp2,p_tmp3); // 1 - d*px*qx*py*qy
	
	mod_div(p_tmp1,p_tmp4,p->x);
	mod_div(p_s,p_tmp3,p->y);
}

// http://hyperelliptic.org/EFD/g1p/auto-twisted.html
// point addition in affine plane for twisted edward curve
void p_add_aff_ted(point *p, point *q){
	// c=1
	// a=-1
	// x=(pxqy+pyqx/(1+ px qx py qy))
	// y=(pyqy+pxqx/(1- px qx py qy))
	
	gmp_mul_mod(p->x,q->y,p_tmp1); // px * qy
	gmp_mul_mod(p->y,q->x,p_tmp2); // py * qx
	gmp_add_mod(p_tmp1,p_tmp2,p_tmp3); // px * qy + py * qx
	
	gmp_mul_mod(p->y,q->y,p_tmp1); // py * qy
	gmp_mul_mod(p->x,q->x,p_tmp2); // px * qx
	gmp_add_mod(p_tmp1,p_tmp2,p_s); // px * qy - -1py * qx

	gmp_mul_mod(p_tmp1,p_tmp2,p_tmp4); // px*qx*py*qy	
	gmp_mul_mod(p_tmp4,dp->d,p_tmp2); // d*px*qx*py*qy
	
	gmp_add_mod(c1,p_tmp2,p_tmp4); // 1 + d*px*qx*py*qy
	gmp_sub_mod(c1,p_tmp2,p_tmp1); // 1 - d*px*qx*py*qy
	
	mod_div(p_tmp3,p_tmp4,p->x);
	mod_div(p_s,p_tmp1,p->y);
}

// point multiplication in affine plane for twisted edward curve with the binary methode
char p_mul_aff_ted_bin(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 msbPos;
	point *baseP;
	
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	cpPoint(baseP,public_key);
	for (i=msbPos;i>=0;i--){
		p_add_aff_ted(public_key,public_key);
		if (getIthBit(private_key,i)){
			p_add_aff_ted(public_key,baseP);
		}
	}
	free(baseP);
	return 0;
}


// "add-2007-bl"
// 2007 Bernstein–Lange.
// http://hyperelliptic.org/EFD/g1p/auto-edwards-projective.html#addition-add-2007-bl-2
// point addition in projective plane with edward curve
void p_add_proj_ed(point *p, point *q){
	// c=1
	gmp_mul_mod(p->z,q->z,p_A); 		// A=Z1*Z2
	gmp_mul_mod(p_A,p_A,p_B); 			// B=A^2
	gmp_mul_mod(p->x,q->x,p_C); 		// C=X1*X2
	gmp_mul_mod(p->y,q->y,p_D); 		// D=Y1*Y2
	gmp_mul_mod(dp->d,p_C,p_tmp1); 		// d*C 
	gmp_mul_mod(p_tmp1,p_D,p_E); 		// E=d*C*D 
	gmp_sub_mod(p_B,p_E,p_F);			// F=B-E
	gmp_add_mod(p_B,p_E,p_G);			// G=B+E
	// X3 = A*F*((X1+Y1)*(X2+Y2)-C-D)
	gmp_add_mod(p->x,p->y,p_tmp1);		// X1+Y1
	gmp_add_mod(q->x,q->y,p_tmp2);		// X2+Y2
	gmp_mul_mod(p_tmp1,p_tmp2,p_tmp3);	// (X1+Y1)*(X2+Y2)
	gmp_sub_mod(p_tmp3,p_C,p_tmp1);		// (X1+Y1)*(X2+Y2)-C
	gmp_sub_mod(p_tmp1,p_D,p_tmp2);		// (X1+Y1)*(X2+Y2)-C-D
	gmp_mul_mod(p_A,p_F,p_tmp1);		//  A*F
	gmp_mul_mod(p_tmp1,p_tmp2,p->x);	// X3= A*F*((X1+Y1)*(X2+Y2)-C-D)
	
	// Y3 = A*G*(D-C)
	gmp_sub_mod(p_D,p_C,p_tmp1);		// (D-C)
	gmp_mul_mod(p_A,p_G,p_tmp2);		// A*G
	gmp_mul_mod(p_tmp1,p_tmp2,p_tmp3);	// A*G*(D-C)
	gmp_sub_mod(dp->p,p_tmp3,p->y);
	// Z3 = c*F*G
	gmp_mul_mod(p_F,p_G,p->z);			// Y3=c*F*G , c=1
	
}

// Edward Curve affine
// #############################################################################################
// #############################################################################################

// ==============================================================================================================================

// ###############################################################################################
// #############################################################################################
// Twisted Edward Curve projective

// http://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#addition-add-2008-bbjlp
// add-2008-bbjlp
// Source: 2008 Bernstein–Birkner–Joye–Lange–Peters http://eprint.iacr.org/2008/013 Section 6. 
// point addition in projective plane with twisted edward curve
void p_add_proj_ted(point *p, point *q, point *r){
	// c=1
	// a=-1
	gmp_mul_mod(p->z,q->z,p_A); 		// A=Z1*Z2
	gmp_mul_mod(p_A,p_A,p_B); 			// B=A^2
	gmp_mul_mod(p->x,q->x,p_C); 		// C=X1*X2
	gmp_mul_mod(p->y,q->y,p_D); 		// D=Y1*Y2
	gmp_mul_mod(dp->d,p_C,p_tmp1); 		// d*C 
	gmp_mul_mod(p_tmp1,p_D,p_E); 		// E=d*C*D 
	gmp_sub_mod(p_B,p_E,p_F);			// F=B-E
	gmp_add_mod(p_B,p_E,p_G);			// G=B+E
	// X3 = A*F*((X1+Y1)*(X2+Y2)-C-D)
	gmp_add_mod(p->x,p->y,p_tmp1);		// X1+Y1
	gmp_add_mod(q->x,q->y,p_tmp2);		// X2+Y2
	gmp_mul_mod(p_tmp1,p_tmp2,p_tmp3);	// (X1+Y1)*(X2+Y2)
	gmp_sub_mod(p_tmp3,p_C,p_tmp1);		// (X1+Y1)*(X2+Y2)-C
	gmp_sub_mod(p_tmp1,p_D,p_tmp2);		// (X1+Y1)*(X2+Y2)-C-D
	gmp_mul_mod(p_A,p_F,p_tmp1);		//  A*F
	gmp_mul_mod(p_tmp1,p_tmp2,r->x);	// X3= A*F*((X1+Y1)*(X2+Y2)-C-D)
	
	// Y3 = A*G*(D-a*C)
	gmp_add_mod(p_D,p_C,p_tmp1);		// (D-a*C), a=-1 => D+C
	gmp_mul_mod(p_A,p_G,p_tmp2);		// A*G
	gmp_mul_mod(p_tmp1,p_tmp2,r->y);	// A*G*(D-a*C)
//	gmp_sub_mod(dp->p,p_tmp3,p->y);
	// Z3 = c*F*G
	gmp_mul_mod(p_F,p_G,r->z);			// Y3=c*F*G , c=1
	
}

char p_mul_proj_ted_bin(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 msbPos;
	point *baseP;
	
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	cpPointZ(baseP,public_key);

	for (i=msbPos;i>=0;i--){
		p_add_proj_ted(public_key,public_key,public_key);
		if (getIthBit(private_key,i)){
			p_add_proj_ted(public_key,baseP,public_key);
		}
	}
	
	transforProjToAff_ed(public_key);
	free(baseP);
	return 0;
}

void transforProjToAff_ed(point *p){
		mod_div(p->x,p->z,p->x);
		mod_div(p->y,p->z,p->y);	
		mod_div(p->z,p->z,p->z);
	}

// pre compute Points for sliding window
void preComputSW_ted(){
	bit_16 i;
	point *baseP;
	
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	cpPointZ(baseP,&oArrPoint[1]);
	cpPointZ(baseP,&oArrPoint[15]);
	
	for (i=2;i<15;i++){
		p_add_proj_ted(&oArrPoint[15],&oArrPoint[1],&oArrPoint[15]);
		cpPointZ(&oArrPoint[15],&oArrPoint[i]);
	}
	// point at infinity
	oArrPoint[0].y[keySizeM1] = 0x1;
	oArrPoint[0].z[keySizeM1] = 0x1;
	p_add_proj_ted(&oArrPoint[15],&oArrPoint[1],&oArrPoint[15]);
	
}

/// handbook of elliptic and hyperelliptic curve cryptography
// angepasste sliding window mehtode, im buch ist der rahmen nicht fest hier schon immer 4 oder 2 bits
char p_mul_proj_ted_rndSW(bit_32 *private_key, point *public_key, bit_32 *rnd1, bit_32 *rnd2){
	
	bit_16 w, i, h; 
	signed int msbPos, msbPosRnd1;
	
	rndPermPoint(rnd1);
	cpPointZ(&arrPoint[rp[0]],public_key);
	msbPos = getMSBPos(private_key,keySize);
	msbPosRnd1 = getMSBPos(rnd2,keySize)-1;
	if ((getIthBit(rnd2,msbPosRnd1)) && (msbPos>=4)){
		w=4;
	}
	else {
		w=2;
	}
	msbPosRnd1--;
	msbPos = msbPos - (msbPos%w);
	
	while (msbPos >= 0){
		
		for (i=0;i<w;i++){
			p_add_proj_ted(public_key,public_key,public_key);	

		}
		h = getIthBits(private_key,msbPos,w);
		
		p_add_proj_ted(public_key,&arrPoint[rp[h]],public_key);
		if ((getIthBit(rnd2,msbPosRnd1)) && (msbPos>=4)){
			w=4;
		}
		else {
			w=2;
		}
		msbPosRnd1--;
		msbPos=msbPos-w;
	}
	transforProjToAff_ed(public_key);
	return 0;
}

char p_mul_proj_ted_rndSca(bit_32 *private_key, point *public_key, bit_32 rnd){
	
	bit_16 i;
	bit_16 msbPos;
	bit_32 *pk;
	pk = calloc(keySize2,sizeof(bit_32));
	point *baseP;
	
	randomize_scalar(private_key, rnd, pk);
	msbPos = getMSBPos(pk,keySize2)-1;
	
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;	
	cpPointZ(baseP,public_key);
	for (i=msbPos;i>=0;i--){
		p_add_proj_ted(public_key,public_key,public_key);
		if (getIthBit_rs(pk,i)){
			p_add_proj_ted(public_key,baseP,public_key);
		}	

	}
	transforProjToAff_ed(public_key);
	free(baseP);
	free(pk);
	return 0;
}

char p_mul_proj_ted_daa(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 it1 = 0, it2;
	bit_16 msbPos;
	point *baseP, *poiP;
	
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	
	poiP = calloc(2,sizeof(point));
	poiP[0].x = calloc(keySize,sizeof(bit_32));
	poiP[0].y = calloc(keySize,sizeof(bit_32));
	poiP[0].z = calloc(keySize,sizeof(bit_32));
	poiP[1].x = calloc(keySize,sizeof(bit_32));
	poiP[1].y = calloc(keySize,sizeof(bit_32));
	poiP[1].z = calloc(keySize,sizeof(bit_32));
	
	cpPointZ(baseP,&poiP[0]);
	cpPointZ(baseP,&poiP[1]);
	p_add_proj_ted(&poiP[1],&poiP[0],&poiP[1]);
	
	for (i=msbPos;i>=0;i--){
		it2 = (1-getIthBit(private_key,i));
		p_add_proj_ted(&poiP[it1],&poiP[it1],&poiP[it1]);
		p_add_proj_ted(&poiP[it1],baseP,&poiP[it2]);
	}
	
	transforProjToAff_ed(&poiP[0]);
	cpPointZ(&poiP[0],public_key);
	
	free(baseP);
	free(poiP);
	return 0;
}

// Resistance Against Differential Power Analysis For Elliptic Curve Cryptosystems
char p_mul_proj_ted_rndCoor(bit_32 *private_key, point *public_key, bit_32 rnd){

	
	bit_16 i;
	bit_16 msbPos;
	point *baseP;
	bit_32 *rnd1;

	rnd1 = calloc(keySize,sizeof(bit_32));
	msbPos = getMSBPos(private_key,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	rnd1[keySizeM1] = rnd;
	
	cpPointZ(baseP,public_key);
	for (i=msbPos;i>=0;i--){
		// randomise coordinate
		// fibonacci-LFSR x^32+x^23+x^7+1
		rnd1[keySizeM1] = ((rnd1[keySizeM1]>>31) ^ (rnd1[keySizeM1]>>22 & 1) ^ (rnd1[keySizeM1]>>6 & 1)) | (rnd1[keySizeM1]<<1 );
		gmp_mul_mod(public_key->x,rnd1,public_key->x);
		gmp_mul_mod(public_key->y,rnd1,public_key->y);
		gmp_mul_mod(public_key->z,rnd1,public_key->z);
		
		p_add_proj_ted(public_key,public_key,public_key);
		if (getIthBit(private_key,i)){
			p_add_proj_ted(public_key,baseP,public_key);
		}
	}
	transforProjToAff_ed(public_key);
	free(baseP);
	free(rnd1);
	return 0;
}

// point multiplication montgomery ladder
char p_mul_proj_ted_ml(bit_32 *private_key, point *public_key){
	
	bit_16 i;
	bit_16 it1, it2;
	bit_16 msbPos;
	point *baseP, *poiP;
	
	msbPos = getMSBPos(private_key,keySize);
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	
	poiP = calloc(2,sizeof(point));
	poiP[0].x = calloc(keySize,sizeof(bit_32));
	poiP[0].y = calloc(keySize,sizeof(bit_32));
	poiP[0].z = calloc(keySize,sizeof(bit_32));
			
	poiP[1].x = calloc(keySize,sizeof(bit_32));
	poiP[1].y = calloc(keySize,sizeof(bit_32));
	poiP[1].z = calloc(keySize,sizeof(bit_32));
	
	poiP[0].y[keySizeM1] = 0x1;
	poiP[0].z[keySizeM1] = 0x1;
	cpPointZ(baseP,&poiP[1]);
	
	for (i=msbPos;i>=0;i--){
		it1 = getIthBit(private_key,i);
		it2 = 1-it1;
		p_add_proj_ted(&poiP[0],&poiP[1],&poiP[it2]);
		p_add_proj_ted(&poiP[it1],&poiP[it1],&poiP[it1]);
		
	}
	transforProjToAff_ed(&poiP[0]);
	cpPointZ(&poiP[0],public_key);
	free(baseP);
	free(poiP);
	return 0;
}

// point multiplication randomised skalar splitting
char p_mul_proj_ted_rss(bit_32 *private_key, point *public_key, bit_32 *rnd){
	
	bit_16 i;
	signed int msbPos1, msbPos2, carry;
	bit_32 *rnd2;
	rnd2 = calloc(keySize,sizeof(bit_32));
	point *baseP, *poiP;
	gmp_sub(private_key,rnd,rnd2,&carry,keySize);
	
	msbPos1 = getMSBPos(rnd,keySize)-1;
	msbPos2 = getMSBPos(rnd2,keySize)-1;
	baseP = calloc(1,sizeof(point));
	baseP->x = dp->g_x;
	baseP->y = dp->g_y;
	baseP->z = c1;
	
	poiP = calloc(2,sizeof(point));
	poiP[0].x = calloc(keySize,sizeof(bit_32));
	poiP[0].y = calloc(keySize,sizeof(bit_32));
	poiP[0].z = calloc(keySize,sizeof(bit_32));
			
	poiP[1].x = calloc(keySize,sizeof(bit_32));
	poiP[1].y = calloc(keySize,sizeof(bit_32));
	poiP[1].z = calloc(keySize,sizeof(bit_32));
	
	cpPointZ(baseP,&poiP[0]);
	cpPointZ(baseP,&poiP[1]);
	
	for (i=msbPos1;i>=0;i--){
		p_add_proj_ted(&poiP[0],&poiP[0],&poiP[0]);
		if (getIthBit(rnd,i)){
			p_add_proj_ted(&poiP[0],baseP,&poiP[0]);
		}
	}
	for (i=msbPos2;i>=0;i--){
		p_add_proj_ted(&poiP[1],&poiP[1],&poiP[1]);
		if (getIthBit(rnd2,i)){
			p_add_proj_ted(&poiP[1],baseP,&poiP[1]);
		}
	}
	p_add_proj_ted(&poiP[0],&poiP[1],public_key);
	transforProjToAff_ed(public_key);
	free(baseP);
	free(poiP);
	return 0;
}

// Twisted Edward Curve projective
// #############################################################################################
// #############################################################################################