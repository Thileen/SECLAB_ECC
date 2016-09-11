/*
eliptic curv domain parameter 

*/


#ifndef EC_DOMAIN_PARAMETER_H
#define EC_DOMAIN_PARAMETER_H


#define BIT_SIZE 32

#define KEY_SIZE_128 128/BIT_SIZE
#define KEY_SIZE_192 192/BIT_SIZE
#define KEY_SIZE_224 224/BIT_SIZE
#define KEY_SIZE_256 256/BIT_SIZE
#define KEY_SIZE_320 320/BIT_SIZE
#define KEY_SIZE_384 384/BIT_SIZE
#define KEY_SIZE_512 512/BIT_SIZE
#define KEY_SIZE_521 (521/BIT_SIZE)+1

// define elliptic curve domain parameter
#define SECP128R1 20

#define SECP192K1 49
#define SECP256K1 50
#define SECP256R1 51

#define P_192     58
#define P_224     59
#define P_256     60
#define P_384     61
#define P_521     62

#define BRAINPOOLP256R1	70
#define BRAINPOOLP320R1 71
#define BRAINPOOLP384R1 72
#define BRAINPOOLP512R1 73

// Edard curve parameter
#define CURVE1174 90
#define ED25519   91

#define CORTEX_M3

#ifdef CYGWIN
	typedef short bit_16;
	typedef unsigned int bit_32;
	typedef long long unsigned int bit_64;
#endif
#ifdef CORTEX_M3
//	typedef uint16_t bit_16;
//	typedef uint32_t bit_32;
//	typedef uint64_t bit_64;
	typedef short bit_16;
	typedef unsigned int bit_32;
	typedef long long unsigned int bit_64;
#endif

typedef struct{
	bit_32 *x;
	bit_32 *y;
	bit_32 *z;
}point;


// 6-tuple T=(p,a,b,G,n,h)
typedef struct {
	const bit_32 key_size;
	const bit_32 *p;
	const bit_32 *a;
	const bit_32 *b;
	const bit_32 *g_x;
	const bit_32 *g_y;
    const bit_32 *n;
	const bit_16 h;
	const bit_32 *d;
} DOMAIN_PARAMETER;


#endif /*EC_DOMAIN_PARAMETER_H*/