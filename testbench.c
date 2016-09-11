
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include "ec_domain_parameter.h"
#include "seclabec.h"



// parameter:
// 1: length = argc
// 2: commando = argv[1]
// 3: data

// MSB...LSB
bit_32 psrnd = 0x12345678;
bit_32 test(){
	bit_32 tmp = psrnd+~psrnd;
	return (psrnd^tmp);
}

int main( int argc, char *argv[] ){
	
	int i,j;
	int shift;
    bit_32 *a, *b, *c, *d, *r;
	bit_16 carry;
	point *poiP, *poiQ;
	if (argc <= 1){
		printf("\n%d",sizeof(long long));
		printf("\n%d",sizeof(long long int));
		
		
		printf("shift: %x",0x1<<31);
		
		return -1;
	} 

	bit_32 aa, bb, cc;
	
	int com = (bit_32)strtol(argv[1], NULL, 16);
	int size = (bit_32)strtol(argv[2], NULL, 16);
	
	//printf("com=%d -- size=%d\n",com,size);

	switch(size){

		case SECP128R1:
			init_ecc(SECP128R1);
			size=4;
			//printf("secp128r1\n");
			break;
		case P_192:
			init_ecc(P_192);
			size=6;
			//printf("p_192\n");
			break;
		case P_224:
			init_ecc(P_224);
			size=7;
			//printf("p_224\n");
			break;
		case P_256:
			init_ecc(P_256);
			size=8;
			//printf("p_256\n");
			break;
		case P_384:
			init_ecc(P_384);
			size=12;
			//printf("p_384\n");
			break;
		case P_521:
			init_ecc(P_521);
			size=16;
			//printf("p_512\n");
			break;
		case BRAINPOOLP256R1:
			init_ecc(BRAINPOOLP256R1);
			size=8;
			//printf("BRAINPOOLP256R1\n");
			break;
		case CURVE1174:
			init_ecc(CURVE1174);
			size=8;
			//printf("curve1174\n");
			break;
		case ED25519:
			init_ecc(ED25519);
			size=8;
			//printf("ed25519\n");
			break;	
		default: break;
	}

	switch (com){
		case 1: //void gmp_add
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
						
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
				b[i] = (bit_32)strtoul(argv[3+size+i], NULL, 16);
			}
			
			//mp_add_mod(a,b,c);
			gmp_add(a,b,c,&carry,size);
			
			printf("%08x ",carry);
			for (i = 0; i<size;i++){
				printf("%08x ",c[i]);
			}
			
			free(a);
			free(b);
			free(c);
			break;
		case 2: // void gmp_sub
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
						
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
				b[i] = (bit_32)strtoul(argv[3+size+i], NULL, 16);
			}
			
			gmp_sub(a,b,c,&carry,size);
		//	printf("%08x ",carry);
			for (i = 0; i<size;i++){
				printf("%08x ",c[i]);
			}
			free(a);
			free(b);
			free(c);
			break;
		case 3: // bit_16 compare
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			//printf("\n %d \n",(sizeof(a)));
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
				b[i] = (bit_32)strtoul(argv[3+size+i], NULL, 16);
			}
			
			printf("%d ",compareAA(a,b,size));
			free(a);
			free(b);
			break;
		case 4: //void gmp_add_mod
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
						
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
				b[i] = (bit_32)strtoul(argv[3+size+i], NULL, 16);
			}
			gmp_add_mod(a,b,c);
			for (i = 0; i<size;i++){
				printf("%08x ",c[i]);
			}
			free(a);
			free(b);
			free(c);
			break;
		case 5: // void gmp_sub_mod
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
						
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
				b[i] = (bit_32)strtoul(argv[3+size+i], NULL, 16);
			}
			
			gmp_sub_mod(a,b,c);
			for (i = 0; i<size;i++){
				printf("%08x ",c[i]);
			}
			free(a);
			free(b);
			free(c);
			break;
		case 6: // void gmp_mul
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size*2,sizeof(bit_32));
						
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
				b[i] = (bit_32)strtoul(argv[3+size+i], NULL, 16);
			}
			
			gmp_mul(a,b,c,size);
			for (i = 0; i<(size*2);i++){
				printf("%08x ",c[i]);
			}
			free(a);
			free(b);
			free(c);
			break;
		case 7: // void gmp_sqrt
			a = calloc(size,sizeof(bit_32));
			c = calloc(size*2,sizeof(bit_32));
						
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
		//	gmp_sqrt(a,c,size);
			for (i = 0; i<(size*2);i++){
				printf("%08x ",c[i]);
			}
			free(a);
			free(c);
			break;
		case 8: // void mp_div
			a = calloc(size*2,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size*2,sizeof(bit_32));
			r = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size*2);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			for (i = 0; i<(size);i++){
				b[i] = (bit_32)strtoul(argv[3+(size*2)+i], NULL, 16);
			}
			
			omp_div(a,size*2,b,size,c,r);
			for (i = 0; i<(size*2);i++){
				printf("%08x ",c[i]);
			}
			for (i = 0; i<(size);i++){
				printf("%08x ",r[i]);
			}
			
			free(a);
			free(b);
			free(c);
			free(r);
			break;
		case 9: // void mp_inv
			a = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
		//	gmp_inv(a,c);
			for (i = 0; i<(size);i++){
				printf("%08x ",c[i]);
			}
			
			free(a);
			free(c);
			break;	
		case 10: // void mp_div_len
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			r = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			for (i = 0; i<(size);i++){
				b[i] = (bit_32)strtoul(argv[3+(size)+i], NULL, 16);
			}
			
			omp_div(a,size,b,size,c,r);
			for (i = 0; i<(size);i++){
				printf("%08x ",c[i]);
			}
			for (i = 0; i<(size);i++){
				printf("%08x ",r[i]);
			}
			
			free(a);
			free(b);
			free(c);
			free(r);
			break;
		case 11: // void mp_mul_mod
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
				b[i] = (bit_32)strtoul(argv[3+(size)+i], NULL, 16);
			}
			
			gmp_mul_mod(a,b,c);
			for (i = 0; i<(size);i++){
				printf("%08x ",c[i]);
			}
			printf("\nfast reduction\n");
			gmp_mul_mod_fr(a,b,c);
			for (i = 0; i<(size);i++){
				printf("%08x ",c[i]);
			}
			
			free(a);
			free(b);
			free(c);
			break;
		case 12: // void point_add
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			d = calloc(size,sizeof(bit_32));

			poiP = calloc(1,sizeof(point));
			poiQ = calloc(1,sizeof(point));

			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
				b[i] = (bit_32)strtoul(argv[3+(size)+i], NULL, 16);
				c[i] = (bit_32)strtoul(argv[3+(size*2)+i], NULL, 16);
				d[i] = (bit_32)strtoul(argv[3+(size*3)+i], NULL, 16);
			}
			
			//printf("%x %x %x %x\n%x %x %x %x\n",a[0],a[1],a[2],a[3],b[0],b[1],b[2],b[3]);
			//printf("a=%x\nb=%x\n",a[0],b[0]);
			
			poiP->x = a;
			poiP->y = b;
			
			poiQ->x = c;
			poiQ->y = d;
			
		//	printf("p=%x %x\nq=%x %x\n",poiP->x[0],poiP->y[0],poiQ->x[0],poiQ->y[0]);
			
			p_add_aff_w(poiP,poiQ);
			
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->x[i]);
			}
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->y[i]);
			}
			
			free(a);
			free(b);
			break;
		case 13: // void point_double
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));

			poiP = calloc(1,sizeof(point));

			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
				b[i] = (bit_32)strtoul(argv[3+(size)+i], NULL, 16);
			}
			
			//printf("%x %x %x %x\n%x %x %x %x\n",a[0],a[1],a[2],a[3],b[0],b[1],b[2],b[3]);
			//printf("a=%x\nb=%x\n",a[0],b[0]);
			
			poiP->x = a;
			poiP->y = b;
			
		//	printf("p=%x %x\n",poiP->x[0],poiP->y[0]);
			p_dbl_aff_w(poiP);
			
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->x[i]);
			}
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->y[i]);
			}
			
			free(a);
			free(b);
			break;
		case 14: // point_mul
			a = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			poiP = calloc(1,sizeof(point));
			poiP->x = calloc(size,sizeof(bit_32));
			poiP->y = calloc(size,sizeof(bit_32));
			poiP->z = calloc(size,sizeof(bit_32));
			
			p_mul_aff_w_bin(a,poiP);
			
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->x[i]);
			}
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->y[i]);
			}
			
			free(a);
			free(poiP);
			break;
		case 15: // void gmp_sqrt_mod
			a = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			poiP = calloc(1,sizeof(point));
			poiP->x = calloc(size,sizeof(bit_32));
			poiP->y = calloc(size,sizeof(bit_32));
			poiP->z = calloc(size,sizeof(bit_32));
			
			p_mul_proj_w_bin(a,poiP);
			
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->x[i]);
			}
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->y[i]);
			}
			
			free(a);
			free(poiP);
			break;
		case 16: // void shiftArrayRightBySteps
			a = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			shiftArrayRightBySteps(a, size, 32);
			for (i = 0; i<(size);i++){
				printf("%08x ",a[i]);
			}
			free(a);
			break;
		case 17: // void shiftArrayLeftBySteps
			a = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			shiftArrayLeftBySteps(a, size, 32);
			for (i = 0; i<(size);i++){
				printf("%08x ",a[i]);
			}
			free(a);
			break;
		case 18: // void init_barrett_reduction
			
			init_barrett_reduction();
			break;	
			
		case 19: // void barrett reduction
			a = calloc(size*2,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size*2);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			init_barrett_reduction();
			br_mod(a, c);
			conclude_barrett_reduction();
			for (i = 0; i<(size);i++){
				printf("%08x ",c[i]);
			}
			free(a);
			break;	
		case 20: // inverse Handbook of elliptic and hyperelliptic cc
			a = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
		//	hmp_inv(a, c);
			
			for (i = 0; i<(size);i++){
				printf("%08x ",c[i]);
			}
			free(a);
			break;	
		case 21: // inverse Handbook of elliptic and hyperelliptic cc
			a = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			gmp_mont_inv(a, c);
			
			for (i = 0; i<(size);i++){
				printf("%08x ",c[i]);
			}
			free(a);
			break;		
		
		case 22: // mod div
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
				b[i] = (bit_32)strtoul(argv[3+(size)+i], NULL, 16);
			}
			
			mod_div(a,b, c);
			
			for (i = 0; i<(size);i++){
				printf("%08x ",c[i]);
			}
			free(a);
			break;			
		case 23: // twisted edward curve 
			a = calloc(size,sizeof(bit_32));
			
			poiP = calloc(1,sizeof(point));
			poiP->x = calloc(size,sizeof(bit_32));
			poiP->y = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			p_mul_aff_ted_bin(a, poiP);
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->x[i]);
			}
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->y[i]);
			}
			free(a);
			free(poiP);
			break;
		case 24: // mul projective 
			a = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			
			poiP = calloc(1,sizeof(point));
			poiP->x = calloc(size,sizeof(bit_32));
			poiP->y = calloc(size,sizeof(bit_32));
			poiP->z = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			p_mul_proj_ted_bin(a, poiP);
			
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->x[i]);
			}
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->y[i]);
			}
			free(a);
			free(c);
			free(poiP);
			break;
		case 25: // randomise sliding window weierstrass 
			a = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			
			poiP = calloc(1,sizeof(point));
			poiP->x = calloc(size,sizeof(bit_32));
			poiP->y = calloc(size,sizeof(bit_32));
			poiP->z = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			preComputSW_w();
			c[0] = 0;
			c[1] = 0xabcdef01;
			c[2] = 0xfda32463;
			c[4] = 0x98873662;
			
			p_mul_proj_w_rndSW(a, poiP, c, c);
			
			printf("\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->z[i]);
			}
			
			free(a);
			free(c);
			free(poiP);
			break;
		case 26: // random scalar 
			a = calloc(size,sizeof(bit_32));
			b = calloc(1,sizeof(bit_32));
			c = calloc((size*2),sizeof(bit_32));
			b[0]=0x10000000;
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			poiP = calloc(1,sizeof(point));
			poiP->x = calloc(size,sizeof(bit_32));
			poiP->y = calloc(size,sizeof(bit_32));
			poiP->z = calloc(size,sizeof(bit_32));
			
//			randomize_scalar(a, b[0],c);
			p_mul_proj_w_rndSca(a,poiP,b[0]);
			
			printf("\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->z[i]);
			}
			free(a);
			free(b);
			free(c);
			free(poiP);
			break;
		case 27: // double-and-add-always 
			a = calloc(size,sizeof(bit_32));
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			poiP = calloc(2,sizeof(point));
			poiP[0].x = calloc(size,sizeof(bit_32));
			poiP[0].y = calloc(size,sizeof(bit_32));
			poiP[0].z = calloc(size,sizeof(bit_32));
			
			poiP[1].x = calloc(size,sizeof(bit_32));
			poiP[1].y = calloc(size,sizeof(bit_32));
			poiP[1].z = calloc(size,sizeof(bit_32));
			
			p_mul_proj_w_daa(a, poiP);
			
			printf("\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			free(a);
			free(poiP);
			break;		
		case 28: // randomise coordinate 
			a = calloc(size,sizeof(bit_32));
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			poiP = calloc(2,sizeof(point));
			poiP[0].x = calloc(size,sizeof(bit_32));
			poiP[0].y = calloc(size,sizeof(bit_32));
			poiP[0].z = calloc(size,sizeof(bit_32));
			
			poiP[1].x = calloc(size,sizeof(bit_32));
			poiP[1].y = calloc(size,sizeof(bit_32));
			poiP[1].z = calloc(size,sizeof(bit_32));
			
			p_mul_proj_w_rndCoor(a, poiP,test);
			
			printf("\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			free(a);
			free(poiP);
			break;	
			
		case 29: // random sliding window  twisted edward
			a = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			
			poiP = calloc(1,sizeof(point));
			poiP->x = calloc(size,sizeof(bit_32));
			poiP->y = calloc(size,sizeof(bit_32));
			poiP->z = calloc(size,sizeof(bit_32));
			
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			preComputSW_ted();
		
			c[0] = 0;
			c[1] = 0xabcdef01;
			c[2] = 0xfda32463;
			c[4] = 0x98873662;
			
			p_mul_proj_ted_rndSW(a, poiP, c, c);
			
			printf("\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->z[i]);
			}
			
			free(a);
			free(c);
			free(poiP);
			break;
		case 30: // random scalar twisted edward kurve 
			a = calloc(size,sizeof(bit_32));
			b = calloc(1,sizeof(bit_32));
			c = calloc((size*2),sizeof(bit_32));
			b[0]=0x10000000;
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			poiP = calloc(1,sizeof(point));
			poiP->x = calloc(size,sizeof(bit_32));
			poiP->y = calloc(size,sizeof(bit_32));
			poiP->z = calloc(size,sizeof(bit_32));
			
			p_mul_proj_ted_rndSca(a,poiP,b[0]);
			
			printf("\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP->z[i]);
			}
			free(a);
			free(b);
			free(c);
			free(poiP);
			break;
		case 31: // double-and-add-always twisted edward kurve 
			a = calloc(size,sizeof(bit_32));
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			poiP = calloc(2,sizeof(point));
			poiP[0].x = calloc(size,sizeof(bit_32));
			poiP[0].y = calloc(size,sizeof(bit_32));
			poiP[0].z = calloc(size,sizeof(bit_32));
			
			poiP[1].x = calloc(size,sizeof(bit_32));
			poiP[1].y = calloc(size,sizeof(bit_32));
			poiP[1].z = calloc(size,sizeof(bit_32));
			
			p_mul_proj_ted_daa(a, poiP);
			
			printf("\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			free(a);
			free(poiP);
			break;		
		case 32: // randomise coordinate twisted edward curve
			a = calloc(size,sizeof(bit_32));
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			poiP = calloc(2,sizeof(point));
			poiP[0].x = calloc(size,sizeof(bit_32));
			poiP[0].y = calloc(size,sizeof(bit_32));
			poiP[0].z = calloc(size,sizeof(bit_32));
			
			poiP[1].x = calloc(size,sizeof(bit_32));
			poiP[1].y = calloc(size,sizeof(bit_32));
			poiP[1].z = calloc(size,sizeof(bit_32));
			
			p_mul_proj_ted_rndCoor(a, poiP,test);
			
			printf("\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			free(a);
			free(poiP);
			break;
			
		case 33: // randomise coordinate twisted edward curve
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			c[0] = test();
			c[1] = test();
			c[2] = test();
			c[3] = test();
			c[4] = test();
			c[5] = test();
			c[6] = test();
			c[7] = test();
			
			for (i=0;i<size;i++){
				if (a[i]!=0){
					b[i] = a[i]-1;
				}
			}
			
			
			preComputSW_w();
			
			poiP = calloc(2,sizeof(point));
			poiP[0].x = calloc(size,sizeof(bit_32));
			poiP[0].y = calloc(size,sizeof(bit_32));
			poiP[0].z = calloc(size,sizeof(bit_32));
			
			poiP[1].x = calloc(size,sizeof(bit_32));
			poiP[1].y = calloc(size,sizeof(bit_32));
			poiP[1].z = calloc(size,sizeof(bit_32));
			
			p_mul_aff_w_bin(a,&poiP[0]);;
			printf("\n\n Weierstrass affine plane binary method for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			
			p_mul_proj_w_bin(a,&poiP[0]);
			printf("\n\n Weierstrass projective plane binary method for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_w_daa(a, &poiP[0]);
			printf("\n\n Weierstrass projective plane double and add always for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_w_ml(a,&poiP[0]);
			printf("\n\n Weierstrass projective plane montgomery-ladder for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_w_rndSca(a,&poiP[0],c[0]);
			printf("\n\n Weierstrass projective plane randomise scalar for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_w_rss(a,&poiP[0], b);
			printf("\n\n Weierstrass projective plane randomise scalar splitting for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}

			p_mul_proj_w_rndSW(a, poiP, c, c);
			printf("\n\n Weierstrass projective plane randomise sliding window for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_w_rndCoor(a, &poiP[0],test);
			printf("\n\n Weierstrass projective plane randomise coordinate for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			free(a);
			free(b);
			free(c);
			free(poiP);
			break;
			
		case 34: // randomise coordinate twisted edward curve
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			c[0] = test();
			c[1] = test();
			c[2] = test();
			c[3] = test();
			c[4] = test();
			c[5] = test();
			c[6] = test();
			c[7] = test();
			
			for (i=0;i<size;i++){
				if (a[i]!=0){
					b[i] = a[i]-1;
				}
			}
			
			
			preComputSW_w();
			
			poiP = calloc(2,sizeof(point));
			poiP[0].x = calloc(size,sizeof(bit_32));
			poiP[0].y = calloc(size,sizeof(bit_32));
			poiP[0].z = calloc(size,sizeof(bit_32));
			
			poiP[1].x = calloc(size,sizeof(bit_32));
			poiP[1].y = calloc(size,sizeof(bit_32));
			poiP[1].z = calloc(size,sizeof(bit_32));
			
			p_mul_aff_w_bin_fr(a,&poiP[0]);;
			printf("\n\n Weierstrass affine plane binary method for P-256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			
			p_mul_proj_w_bin_fr(a,&poiP[0]);
			printf("\n\n Weierstrass projective plane binary method for P-256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_w_daa_fr(a, &poiP[0]);
			printf("\n\n Weierstrass projective plane double and add always for P-256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_w_ml_fr(a,&poiP[0]);
			printf("\n\n Weierstrass projective plane montgomery-ladder for P-256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_w_rndSca_fr(a,&poiP[0],c[0]);
			printf("\n\n Weierstrass projective plane randomise scalar for P-256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_w_rss_fr(a,&poiP[0], b);
			printf("\n\n Weierstrass projective plane randomise scalar splitting for P-256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}

			p_mul_proj_w_rndSW_fr(a, poiP, c, c);
			printf("\n\n Weierstrass projective plane randomise sliding window for P-256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_w_rndCoor_fr(a, &poiP[0],test);
			printf("\n\n Weierstrass projective plane randomise coordinate for P-256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			free(a);
			free(b);
			free(c);
			free(poiP);
			break;
			
		case 35: // randomise coordinate twisted edward curve
			a = calloc(size,sizeof(bit_32));
			b = calloc(size,sizeof(bit_32));
			c = calloc(size,sizeof(bit_32));
			for (i = 0; i<(size);i++){
				a[i] = (bit_32)strtoul(argv[3+i], NULL, 16);
			}
			
			c[0] = test();
			c[1] = test();
			c[2] = test();
			c[3] = test();
			c[4] = test();
			c[5] = test();
			c[6] = test();
			c[7] = test();
			
			for (i=0;i<size;i++){
				if (a[i]!=0){
					b[i] = a[i]-2;
				}
			}
			
			preComputSW_ted();
			
			poiP = calloc(2,sizeof(point));
			poiP[0].x = calloc(size,sizeof(bit_32));
			poiP[0].y = calloc(size,sizeof(bit_32));
			poiP[0].z = calloc(size,sizeof(bit_32));
			
			poiP[1].x = calloc(size,sizeof(bit_32));
			poiP[1].y = calloc(size,sizeof(bit_32));
			poiP[1].z = calloc(size,sizeof(bit_32));
			
			p_mul_aff_ted_bin(a,&poiP[0]);;
			printf("\n\n TED affine plane binary method for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			
			p_mul_proj_ted_bin(a,&poiP[0]);
			printf("\n\n TED projective plane binary method for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_ted_daa(a, &poiP[0]);
			printf("\n\n TED projective plane double and add always for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_ted_ml(a,&poiP[0]);
			printf("\n\n TED projective plane montgomery-ladder for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_ted_rndSca(a,&poiP[0],c[0]);
			printf("\n\n TED projective plane randomise scalar for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_ted_rss(a,&poiP[0], b);
			printf("\n\n TED projective plane randomise scalar splitting for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}

			p_mul_proj_ted_rndSW(a, poiP, c, c);
			printf("\n\n TED projective plane randomise sliding window for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			p_mul_proj_ted_rndCoor(a, &poiP[0],test);
			printf("\n\n TED projective plane randomise coordinate for brainpool256\nx=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].x[i]);
			}
			printf("\ny=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].y[i]);
			}
			printf("\nz=");
			for (i = 0; i<(size);i++){
				printf("%08x ",poiP[0].z[i]);
			}
			
			free(a);
			free(b);
			free(c);
			free(poiP);
			break;
	}
	return 0;
}