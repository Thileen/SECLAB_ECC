/*
skalarae punkt multiplikation eines ECCs 

*/

#include "ec_domain_parameter.h"

#ifndef SECLABEC_H
#define SECLABEC_H

void shiftArrayRightBySteps(bit_32 *a, bit_16 len, bit_16 steps);
void shiftArrayLeftBySteps(bit_32 *a, bit_16 len, bit_16 steps);

static const DOMAIN_PARAMETER *dp;

void gmp_add(bit_32 *a, bit_32 *b, bit_32 *c, bit_16 *carry, bit_16 len);
void gmp_sub(bit_32 *a, bit_32 *b, bit_32 *c, bit_16 *carry, bit_16 len);

void gmp_add_mod(bit_32 *a, bit_32 *b, bit_32 *c);
void gmp_sub_mod(bit_32 *a, bit_32 *b, bit_32 *c);

void gmp_mul(bit_32 *a, bit_32 *b, bit_32 *c, bit_16 len);
void gmp_mul_mod(bit_32 *a, bit_32 *b, bit_32 *c);
//void mp_div(bit_32 *a, bit_32 *b, bit_32 *c, bit_32 *r);
void omp_div(bit_32 *a, bit_16 lenA, bit_32 *b, bit_16 lenB, bit_32 *c, bit_32 *r);

// multiplication modulo fast reduction p256
void gmp_mul_mod_fr(bit_32 *a, bit_32 *b, bit_32 *c);


void mp_add_mod(bit_32 *a, bit_32 *b, bit_32 *c);
void mp_sub(bit_32 *a, bit_32 *b, bit_32 *c);
bit_16 compareAA(bit_32 *a, bit_32 *b, bit_16 len);
bit_16 compareAI(bit_32 *a, bit_32 b, bit_16 len);

void init_barrett_reduction();
void conclude_barrett_reduction();
void br_mod(bit_32 *a, bit_32 *c);
void free_memory();
void init_memory();
void gmp_mont_inv(bit_32 *a, bit_32 *c);
void mod_div(bit_32 *a, bit_32 *b, bit_32 *c);

bit_16 init_ecc(int name);
void conclude_ecc();
// input p and q
// output p
void p_add_aff_w(point *p, point *q);
// input p 
// output p
void p_dbl_aff_w(point *p);

char p_mul_aff_w_bin(bit_32 *private_key, point *public_key);
char p_mul_aff_w_bin_fr(bit_32 *private_key, point *public_key);

void point_add_projective(point *p, point *q);
void point_double_projective(point *p);
// point multiplication affine plane
char p_mul_proj_w_bin(bit_32 *private_key, point *public_key);

// point multiplicatio affine plane with fast reduction p256
char p_mul_proj_w_bin_fr(bit_32 *private_key, point *public_key);


// edward curve
// twisted edward curve a=-1
char p_mul_aff_ted_bin(bit_32 *private_key, point *public_key);
char p_mul_proj_ted_bin(bit_32 *private_key, point *public_key);

// randomise Sliding Window
void preComputSW_w();
void preComputSW_ted();
void rndPermPoint(bit_32 *rndNum);
char p_mul_proj_w_rndSW(bit_32 *private_key, point *public_key, bit_32 *rnd1, bit_32 *rnd2);
char p_mul_proj_w_rndSW_fr(bit_32 *private_key, point *public_key, bit_32 *rnd1, bit_32 *rnd2);
char p_mul_proj_ted_rndSW(bit_32 *private_key, point *public_key, bit_32 *rnd1, bit_32 *rnd2);

// randomise scalar
void randomize_scalar(bit_32 *private_key, bit_32 rnd, bit_32 *rs_private_key);
char p_mul_proj_w_rndSca(bit_32 *private_key, point *public_key, bit_32 rnd);
char p_mul_proj_w_rndSca_fr(bit_32 *private_key, point *public_key, bit_32 rnd);
char p_mul_proj_ted_rndSca(bit_32 *private_key, point *public_key, bit_32 rnd);

// double-and_add_always
char p_mul_proj_w_daa(bit_32 *private_key, point *public_key);
char p_mul_proj_w_daa_fr(bit_32 *private_key, point *public_key);
char p_mul_proj_ted_daa(bit_32 *private_key, point *public_key);

// randomise coordinate
char p_mul_proj_w_rndCoor(bit_32 *private_key, point *public_key, bit_32 rnd);
char p_mul_proj_w_rndCoor_fr(bit_32 *private_key, point *public_key, bit_32 rnd);
char p_mul_proj_ted_rndCoor(bit_32 *private_key, point *public_key, bit_32 rnd);

// montgomery ladder
char p_mul_proj_w_ml(bit_32 *private_key, point *public_key);
char p_mul_proj_w_ml_fr(bit_32 *private_key, point *public_key);
char p_mul_proj_ted_ml(bit_32 *private_key, point *public_key);

// randomised skalar splitting
char p_mul_proj_w_rss(bit_32 *private_key, point *public_key, bit_32 *rnd);
char p_mul_proj_w_rss_fr(bit_32 *private_key, point *public_key, bit_32 *rnd);
char p_mul_proj_ted_rss(bit_32 *private_key, point *public_key, bit_32 *rnd);
#endif