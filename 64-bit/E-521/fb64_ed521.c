/* Test program for ed521 scalar point multiplication
 The following information are taken from http://safecurves.cr.yp.to/base.html
 Curve Equation :: x^2+y^2 = 1-376014x^2y^2 where d=-376014
 Uses Daniel J. Bernstein et al.s point multiplication method from Curve41417 paper
 Cache safety thanks to ed25519
 Fully Tested and debugged
 We have replaced gmul() by TMV_product()
 We have added rdtscp() to measure the clock cycles count as suggested by Paoloni; http://www.intel.com.tr/content/dam/www/public/us/en/documents/white-papers/ia-32-ia-64-benchmark-code-execution-paper.pdf
 We have amended scr(), gmuli(), gsqr(), and gsqr2() according to modulus p
 gcc -Wall -O3 fb64_ed521.c -o fb64_ed521.exe */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#define WINDOW 5 //6 //4

#define BIT_LENGTH 519

#if WINDOW==4
#define V 3
#define E (int)ceil((double)BIT_LENGTH/(WINDOW*V))
#define D E*V
#define L D*WINDOW
#endif

#if WINDOW==5
#define V 3
#define E (int)ceil((double)BIT_LENGTH/(WINDOW*V))
#define D E*V
#define L D*WINDOW
#endif

#if WINDOW==6
#define V 3
#define E (int)ceil((double)BIT_LENGTH/(WINDOW*V))
#define D E*V
#define L D*WINDOW
#endif

#define M (1<<(WINDOW-1))
#define LIMBS 9
#define CURVE_PARAMETER 376014

//#define TEST  /* define to multiply by the group order */

typedef __int128 type128;
typedef int64_t type64;


static const type64 bot58bits = 0x3ffffffffffffff;
static const type64 bot57bits = 0x1ffffffffffffff;

__inline__ uint64_t rdtsc(){
   uint32_t lo, hi;
   __asm__ __volatile__ ("xorl %%eax,%%eax \n       cpuid" ::: "%rax", "%rbx", "%rcx", "%rdx");
   __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
   return (uint64_t)hi << 32 | lo;
}

__inline__ uint64_t rdtscp(){
   uint32_t lo, hi;
   __asm__ __volatile__ ("rdtscp": "=a" (lo), "=d" (hi) :: "%rcx");
   __asm__ __volatile__ ("cpuid"::: "%rax", "%rbx", "%rcx", "%rdx");
   return (uint64_t)hi << 32 | lo;
}

//////////////////////////////////////////////////////////////////////////////
// w=x+y
void gadd(const type64 x[], const type64 y[], type64 w[]){
	w[0]=x[0]+y[0];
	w[1]=x[1]+y[1];
	w[2]=x[2]+y[2];
	w[3]=x[3]+y[3];
	w[4]=x[4]+y[4];
	w[5]=x[5]+y[5];
	w[6]=x[6]+y[6];
	w[7]=x[7]+y[7];
	w[8]=x[8]+y[8];
}

// w=x-y
void gsub(const type64 x[], const type64 y[], type64 w[]){
	w[0]=x[0]-y[0];
	w[1]=x[1]-y[1];
	w[2]=x[2]-y[2];
	w[3]=x[3]-y[3];
	w[4]=x[4]-y[4];
	w[5]=x[5]-y[5];
	w[6]=x[6]-y[6];
	w[7]=x[7]-y[7];
	w[8]=x[8]-y[8];
}

// w-=x
/*void gdec(const type64 x[], type64 w[]){
	w[0]-=x[0];
	w[1]-=x[1];
	w[2]-=x[2];
	w[3]-=x[3];
	w[4]-=x[4];
	w[5]-=x[5];
	w[6]-=x[6];
	w[7]-=x[7];
	w[8]-=x[8];
}*/

// w=x
void gcopy(const type64 x[], type64 w[]){
	w[0]=x[0];
	w[1]=x[1];
	w[2]=x[2];
	w[3]=x[3];
	w[4]=x[4];
	w[5]=x[5];
	w[6]=x[6];
	w[7]=x[7];
	w[8]=x[8];
}

//w=w-x-y
void gtsb(const type64 x[], const type64 y[], type64 w[]){
	w[0]-=x[0]+y[0];
	w[1]-=x[1]+y[1];
	w[2]-=x[2]+y[2];
	w[3]-=x[3]+y[3];
	w[4]-=x[4]+y[4];
	w[5]-=x[5]+y[5];
	w[6]-=x[6]+y[6];
	w[7]-=x[7]+y[7];
	w[8]-=x[8]+y[8];
}

// reduce w - Short Coefficient Reduction
void scr(type64 w[]){
	type64 t0,t1;
	t0=w[0]&bot58bits;

	t1=w[1]+(w[0]>>58);             w[1]=t1&bot58bits;
	t1=w[2]+(t1>>58);           	w[2]=t1&bot58bits;
	t1=w[3]+(t1>>58);           	w[3]=t1&bot58bits;
	t1=w[4]+(t1>>58);           	w[4]=t1&bot58bits;
	t1=w[5]+(t1>>58);           	w[5]=t1&bot58bits;
	t1=w[6]+(t1>>58);           	w[6]=t1&bot58bits;
	t1=w[7]+(t1>>58);           	w[7]=t1&bot58bits;
	t1=w[8]+(t1>>58);           	w[8]=t1&bot57bits;				//w[8] is 57-bit

	w[0]=t0+(t1>>57);				//w[0] can be at most 59-bit
}

// multiply w by a constant, w*=i especially used for multiplication by the curve constant

void gmuli(type64 w[],int i){
	type128 t;

	t=(type128)w[0]*i;                  w[0]=((type64)t)&bot58bits;
	t=(type128)w[1]*i+(t>>58);      	w[1]=((type64)t)&bot58bits;
	t=(type128)w[2]*i+(t>>58);      	w[2]=((type64)t)&bot58bits;
	t=(type128)w[3]*i+(t>>58);      	w[3]=((type64)t)&bot58bits;
	t=(type128)w[4]*i+(t>>58);          w[4]=((type64)t)&bot58bits;
	t=(type128)w[5]*i+(t>>58);      	w[5]=((type64)t)&bot58bits;
	t=(type128)w[6]*i+(t>>58);      	w[6]=((type64)t)&bot58bits;
	t=(type128)w[7]*i+(t>>58);      	w[7]=((type64)t)&bot58bits;
	t=(type128)w[8]*i+(t>>58);      	w[8]=((type64)t)&bot57bits;				//w[8] is 57-bit

	w[0]+=(type64)(t>>57);					//modulus p
}

// z=x^2
// For modulus p the changes are commented
void gsqr(const type64 x[], type64 z[]){
	type128 t0,t1,t2;

	t1=2*((type128)x[0]*x[8]+(type128)x[1]*x[7]+(type128)x[2]*x[6]+(type128)x[3]*x[5])+(type128)x[4]*x[4];
	t0=((type64) t1)&bot57bits;				//z[8] is 57-bit
	t2=4*((type128)x[1]*x[8]+(type128)x[2]*x[7]+(type128)x[3]*x[6]+(type128)x[4]*x[5])+(type128)x[0]*x[0]+(t1>>57);			//modulus p
	z[0]=((type64) t2)&bot58bits;
	t1=4*((type128)x[2]*x[8]+(type128)x[3]*x[7]+(type128)x[4]*x[6])+2*((type128)x[0]*x[1]+(type128)x[5]*x[5])+(t2>>58);
	z[1]=((type64) t1)&bot58bits;
	t2=4*((type128)x[3]*x[8]+(type128)x[4]*x[7]+(type128)x[5]*x[6])+2*(type128)x[0]*x[2]+(type128)x[1]*x[1]+(t1>>58);
	z[2]=((type64) t2)&bot58bits;
	t1=4*((type128)x[4]*x[8]+(type128)x[5]*x[7])+2*((type128)x[0]*x[3]+(type128)x[1]*x[2]+(type128)x[6]*x[6])+(t2>>58);
	z[3]=((type64) t1)&bot58bits;
	t2=4*((type128)x[5]*x[8]+(type128)x[6]*x[7])+2*((type128)x[0]*x[4]+(type128)x[1]*x[3])+(type128)x[2]*x[2]+(t1>>58);
	z[4]=((type64) t2)&bot58bits;
	t1=4*(type128)x[6]*x[8]+2*((type128)x[0]*x[5]+(type128)x[1]*x[4]+(type128)x[2]*x[3]+(type128)x[7]*x[7])+(t2>>58);
	z[5]=((type64) t1)&bot58bits;
	t2=4*(type128)x[7]*x[8]+2*((type128)x[0]*x[6]+(type128)x[1]*x[5]+(type128)x[2]*x[4])+(type128)x[3]*x[3]+(t1>>58);
	z[6]=((type64) t2)&bot58bits;
	t1=2*((type128)x[0]*x[7]+(type128)x[1]*x[6]+(type128)x[2]*x[5]+(type128)x[3]*x[4]+(type128)x[8]*x[8])+(t2>>58);
	z[7]=((type64) t1)&bot58bits;
	t0+=(t1>>58);
	z[8]=((type64)t0)&bot57bits;			//z[8] is 57-bit
	z[0]+=(type64)(t0>>57);					//modulus p
}


// z=2x^2
// For modulus p the changes are commented
void gsqr2(const type64 x[], type64 z[]){
	type128 t0,t1,t2;

	t1=4*((type128)x[0]*x[8]+(type128)x[1]*x[7]+(type128)x[2]*x[6]+(type128)x[3]*x[5])   +(type128)x[4]*(2*x[4]);
	t0=((type64) t1)&bot57bits;				//z[8] is 57-bit
	t2=8*((type128)x[1]*x[8]+(type128)x[2]*x[7]+(type128)x[3]*x[6]+(type128)x[4]*x[5])   +2*((type128)x[0]*x[0])+(t1>>57);		//modulus p
	z[0]=((type64) t2)&bot58bits;
	t1=8*((type128)x[2]*x[8]+(type128)x[3]*x[7]+(type128)x[4]*x[6])     +4*((type128)x[0]*x[1]+(type128)x[5]*x[5])  +(t2>>58);
	z[1]=((type64) t1)&bot58bits;
	t2=8*((type128)x[3]*x[8]+(type128)x[4]*x[7]+(type128)x[5]*x[6])     +4*(type128)x[0]*x[2]+(type128)x[1]*(2*x[1])+(t1>>58);
	z[2]=((type64) t2)&bot58bits;
	t1=8*((type128)x[4]*x[8]+(type128)x[5]*x[7])      +4*((type128)x[0]*x[3]+(type128)x[1]*x[2]+(type128)x[6]*x[6])  +(t2>>58);
	z[3]=((type64) t1)&bot58bits;
	t2=8*((type128)x[5]*x[8]+(type128)x[6]*x[7])  +4*((type128)x[0]*x[4]+(type128)x[1]*x[3])  +(type128)x[2]*(2*x[2])  +(t1>>58);
	z[4]=((type64) t2)&bot58bits;
	t1=8*(type128)x[6]*x[8]  +4*((type128)x[0]*x[5]+(type128)x[1]*x[4]+(type128)x[2]*x[3]+(type128)x[7]*x[7])  +(t2>>58);
	z[5]=((type64) t1)&bot58bits;
	t2=8*(type128)x[7]*x[8]  +4*((type128)x[0]*x[6]+(type128)x[1]*x[5]+(type128)x[2]*x[4])  +(type128)x[3]*(2*x[3])  +(t1>>58);
	z[6]=((type64) t2)&bot58bits;
	t1=4*((type128)x[0]*x[7]+(type128)x[1]*x[6]+(type128)x[2]*x[5]+(type128)x[3]*x[4]+(type128)x[8]*x[8])  +(t2>>58);
	z[7]=((type64) t1)&bot58bits;
	t0+=(t1>>58);
	z[8]=((type64)t0)&bot57bits;			//z[8] is 57-bit
	z[0]+=(type64)(t0>>57);					//modulus p
}


// Hybrid version
/*Residue multiplication modulo 521-bit Mersenne prime
XY = Z (mod pow(2,521) - 1) where X and Y comprise of 9-limb
X[i] and Y[i] are unsigned 58-bit where i=1,2,3,4,5,6,7
X[0], Y[0], and Z[0] are at most unsigned 59-bit
X[8], Y[8], and Z[8] are unsigned 57-bit*/
void TMV_product(const type64 X[], const type64 Y[], type64 Z[]){
	type64 T1[5]={0,0,0,2*X[6],2*X[5]}, T6[5];
	type64 T5=2*X[8], c=2*X[7];

	T1[0]=Y[0]-Y[3];		T1[1]=Y[1]-Y[4];		T1[2]=Y[2]-Y[5];
	type128 R0=(type128)X[3]*T1[0]+(type128)X[2]*T1[1]+(type128)X[1]*T1[2];
	type128 R1 = (type128)X[4]*T1[0]+(type128)X[3]*T1[1]+(type128)X[2]*T1[2];
	type128 R2 = (type128)X[5]*T1[0]+(type128)X[4]*T1[1]+(type128)X[3]*T1[2];

	T1[0]=Y[3]-Y[6];		T1[1]=Y[4]-Y[7];		T1[2]=Y[5]-Y[8];
	type128 R6 = (type128)T1[3]*T1[0]+(type128)T1[4]*T1[1]+(type128)(2*X[4])*T1[2];
	type128 R7 = (type128)c*T1[0]+(type128)T1[3]*T1[1]+(type128)T1[4]*T1[2];
	type128 R8 = (type128)T5*T1[0]+(type128)c*T1[1]+(type128)T1[3]*T1[2];

	T1[0]=Y[0]-Y[6];		T1[1]=Y[1]-Y[7];		T1[2]=Y[2]-Y[8];
	type128 R3 = (type128)X[0]*T1[0]+(type128)T5*T1[1]+(type128)c*T1[2];
	type128 R4 = (type128)X[1]*T1[0]+(type128)X[0]*T1[1]+(type128)T5*T1[2];
	type128 R5 = (type128)X[2]*T1[0]+(type128)X[1]*T1[1]+(type128)X[0]*T1[2];

	T6[0]=X[4]+X[1];		T6[1]=X[5]+X[2];		T6[2]=X[6]+X[3];
	T6[3]=X[7]+X[4];		T6[4]=X[8]+X[5];
	T1[0]=T6[0]+c;			T1[1]=T6[1]+T5;			T1[2]=T6[2]+X[0];
	T1[3]=T6[3]+X[1];		T1[4]=T6[4]+X[2];
	T6[0]+=T1[0];			T6[1]+=T1[1];			T6[2]+=T1[2];
	T6[3]+=T1[3];			T6[4]+=T1[4];
	T5=T1[2]+X[6];

	type128 S = ((type128)T1[2]*Y[2])+((type128)T1[3]*Y[1])+((type128)T1[4]*Y[0]) - R2 - R5;
	c = ((type64)S)&bot57bits;
	S = (((type128)T6[0]*Y[8])+((type128)T6[1]*Y[7])+((type128)T6[2]*Y[6]) + R3 + R6) + (S>>57);
	Z[0] = ((type64)S)&bot58bits;
	S = (((type128)T6[1]*Y[8])+((type128)T6[2]*Y[7])+((type128)T6[3]*Y[6]) + R4 + R7) + (S>>58);
	Z[1] = ((type64)S)&bot58bits;
	S = (((type128)T6[2]*Y[8])+((type128)T6[3]*Y[7])+((type128)T6[4]*Y[6]) + R5 + R8) + (S>>58);
	Z[2] = ((type64)S)&bot58bits;
	S = (((type128)T6[3]*Y[5])+((type128)T6[4]*Y[4])+((type128)T5*Y[3]) + R0 - R6) + (S>>58);
	Z[3] = ((type64)S)&bot58bits;
	S = (((type128)T6[4]*Y[5])+((type128)T5*Y[4])+((type128)T1[0]*Y[3]) + R1 - R7) + (S>>58);
	Z[4] = ((type64)S)&bot58bits;
	S = (((type128)T5*Y[5])+((type128)T1[0]*Y[4])+((type128)T1[1]*Y[3]) + R2 - R8) + (S>>58);
	Z[5] = ((type64)S)&bot58bits;
	S = (((type128)T1[0]*Y[2])+((type128)T1[1]*Y[1])+((type128)T1[2]*Y[0]) - R0 - R3) + (S>>58);
	Z[6] = ((type64)S)&bot58bits;
	S = (((type128)T1[1]*Y[2])+((type128)T1[2]*Y[1])+((type128)T1[3]*Y[0]) - R1 - R4) + (S>>58);
	Z[7] = ((type64)S)&bot58bits;
	c += ((type64)(S>>58));
	Z[8] = c&bot57bits;
	Z[0]+= (c>>57);

}


// Inverse x = 1/x = x^(p-2) mod p
void ginv(type64 x[]){
	type64 x127[LIMBS], w[LIMBS], t[LIMBS], z[LIMBS];
	gsqr(x,x127);       // x127=x^2
	TMV_product(x127,x,t);     // t=x^3
	gsqr(t,x127);       // x127=x^6
	TMV_product(x127,x,w);     // w=x^7
	gsqr(w,x127);       //
	gsqr(x127,t);
	gsqr(t,x127);       // x127=x^56
	gcopy(x127,t);		// t=x^56
	TMV_product(w,t,x127);     // x127=x^63
	gsqr(x127,t);
	TMV_product(t,x,x127);     // x127=x^127

	gsqr(x127,t);
	TMV_product(t,x,z);        // z=x^255

	gcopy(z,w);
	for (int i=0;i<4;i++)
	{
		gsqr(z,t);
		gsqr(t,z);
	}
	TMV_product(z,w,t);        // z=z16

	gcopy(t,w);
	for (int i=0;i<8;i++)
	{
		gsqr(t,z);
		gsqr(z,t);
	}
	TMV_product(t,w,z);        // z=z32

	gcopy(z,w);
	for (int i=0;i<16;i++)
	{
		gsqr(z,t);
		gsqr(t,z);
	}
	TMV_product(z,w,t);        // z=z64

	gcopy(t,w);
	for (int i=0;i<32;i++)
	{
		gsqr(t,z);
		gsqr(z,t);
	}
	TMV_product(t,w,z);        // z=z128

	gcopy(z,w);
	for (int i=0;i<64;i++)
	{
		gsqr(z,t);
		gsqr(t,z);
	}
	TMV_product(z,w,t);        // z=z256

	gcopy(t,w);
	for (int i=0;i<128;i++)
	{
		gsqr(t,z);
		gsqr(z,t);
	}
	TMV_product(t,w,z);      // z=z512

	gsqr(z,t);
	gsqr(t,z);
	gsqr(z,t);
	gsqr(t,z);
	gsqr(z,t);
	gsqr(t,z);
	gsqr(z,t);
	TMV_product(t,x127,z);
	gsqr(z,t);
	gsqr(t,z);
	TMV_product(z,x,t);
	gcopy(t,x);
}

////////////////////////////////////////////////////////////////////////////////////////////
// Point Structure such the points are in extended projective coordinate
typedef struct{
    type64 x[LIMBS];
    type64 y[LIMBS];
    type64 z[LIMBS];
    type64 t[LIMBS];
} ECp;

//P=0, Setting the point as neutral element
void inf(ECp *P){
	for (int i=0; i<LIMBS; i++)
		P->x[i]=P->y[i]=P->z[i]=P->t[i]=0;
	P->y[0]=P->z[0]=1;
}

// Initialise P
void init(type64 *x,type64 *y, ECp *P){
	for (int i=0; i<LIMBS; i++){
		P->x[i]=x[i];
		P->y[i]=y[i];
		P->z[i]=0;
	}
	P->z[0]=1;
	TMV_product(x,y,P->t);
}

//P=Q
void copy(ECp *Q,ECp *P){
	for (int i=0; i<LIMBS; i++){
		P->x[i]=Q->x[i];
		P->y[i]=Q->y[i];
		P->z[i]=Q->z[i];
		P->t[i]=Q->t[i];
	}
}

// P=-Q
void neg(ECp *Q,ECp *P){
	for (int i=0; i<LIMBS; i++){
		P->x[i]=-Q->x[i];
		P->y[i]=Q->y[i];
		P->z[i]=Q->z[i];
		P->t[i]=-Q->t[i];
	}
}


// P=2P, when P->t is not required
void doubl_2(ECp *P){
	type64 a[LIMBS], b[LIMBS], c[LIMBS], e[LIMBS], f[LIMBS], g[LIMBS], h[LIMBS];
	gsqr(P->x,a);
	gsqr(P->y,b);
	gsqr2(P->z,c);
	gadd(P->x,P->y,g); gsqr(g,e); gtsb(a,b,e); //gdec(a,e); gdec(b,e);
	gadd(a,b,g);
	gsub(g,c,f);
	gsub(a,b,h);
	TMV_product(e,f,P->x);
	TMV_product(g,h,P->y);
	TMV_product(f,g,P->z);
}

// P=2P
void doubl_3(ECp *P){
	type64 a[LIMBS], b[LIMBS], c[LIMBS], e[LIMBS], f[LIMBS], g[LIMBS], h[LIMBS];
	gsqr(P->x,a);
	gsqr(P->y,b);
	gsqr2(P->z,c);
	gadd(P->x,P->y,g); gsqr(g,e); gtsb(a,b,e); //gdec(a,e); gdec(b,e);
	gadd(a,b,g);
	gsub(g,c,f);
	gsub(a,b,h);
	TMV_product(e,f,P->x);
	TMV_product(g,h,P->y);
	TMV_product(f,g,P->z);
	TMV_product(e,h,P->t);
}

//P+=Q, when P->t is not required
void add_2(ECp *Q,ECp *P){
	type64 a[LIMBS], b[LIMBS], c[LIMBS], d[LIMBS], e[LIMBS], f[LIMBS], g[LIMBS], h[LIMBS];
	TMV_product(P->x,Q->x,a);
	TMV_product(P->y,Q->y,b);
	TMV_product(P->t,Q->t,c);
	TMV_product(P->z,Q->z,d);
	gadd(d,c,f);  /* reversed sign as d is negative */
	gsub(d,c,g);
	gsub(b,a,h);
	gadd(P->x,P->y,c); gadd(Q->x,Q->y,d); TMV_product(c,d,e); gtsb(a,b,e); //gdec(a,e); gdec(b,e);
	TMV_product(e,f,P->x);
	TMV_product(g,h,P->y);
	TMV_product(f,g,P->z);
}

//P+=Q, when P->t is required
void add_3(ECp *Q,ECp *P){
	type64 a[LIMBS], b[LIMBS], c[LIMBS], d[LIMBS], e[LIMBS], f[LIMBS], g[LIMBS], h[LIMBS];
	TMV_product(P->x,Q->x,a);
	TMV_product(P->y,Q->y,b);
	TMV_product(P->t,Q->t,c);
	TMV_product(P->z,Q->z,d);
	gadd(d,c,f);  /* reversed sign as d is negative */
	gsub(d,c,g);
	gsub(b,a,h);
	gadd(P->x,P->y,c); gadd(Q->x,Q->y,d); TMV_product(c,d,e); gtsb(a,b,e); //gdec(a,e); gdec(b,e);
	TMV_product(e,f,P->x);
	TMV_product(g,h,P->y);
	TMV_product(f,g,P->z);
	TMV_product(e,h,P->t);
}

//From Extended projective to Affine coordinate
void norm(ECp *P){
	type64 w[LIMBS], t[LIMBS];
	gcopy(P->z,w);
	ginv(w);
	TMV_product(P->x, w, t);        //t is in intermediate form
	scr(t);                         //t is a unique residue in modulo p
	gcopy(t,P->x);
	TMV_product(P->y,w,t);          //t is in intermediate form
	scr(t);                         //t is a unique residue in modulo p
	gcopy(t,P->y);
	/*TMV_product(P->z,w,t);          //t is in intermediate form
	scr(t);                         //t is a unique residue in modulo p
	gcopy(t,P->z);
	TMV_product(P->t,w,t);          //t is in intermediate form
	scr(t);                         //t is a unique residue in modulo p
	gcopy(t,P->t);*/
}

void output(ECp *P){
	int i;
	for(i=0; i<LIMBS; i++)
		printf("x[%d] = %lx\n", i, P->x[i]);
	puts("");
	for(i=0; i<LIMBS; i++)
		printf("y[%d] = %lx\n", i, P->y[i]);
	puts("");


}

/* Off-line Pre-computation by taking points in extended projective co-ordinate */
void Offline_precomp(ECp *P, ECp W[][V]){
    int32_t i, j, ti;
	ECp Q;
	ECp exp_d[WINDOW-1], exp_e[V];

	copy(P, &Q);
	copy(P, &exp_e[0]);
	for(i=1; i<V; i++){
        for(j=0; j<E-1; j++)
            doubl_2(&Q);            //Only for doubling
        doubl_3(&Q);                //Later will be used for addition
        copy(&Q, &exp_e[i]);
    }
    copy(P, &Q);
    for(i=1; i<WINDOW; i++){
        for(j=0; j<D-1; j++)
            doubl_2(&Q);            //Only for doubling
        doubl_3(&Q);                //Later will be used for addition
        copy(&Q, &exp_d[i-1]);
    }

    for(i=0; i<M; i++){
        for(j=0; j<V; j++){
            ti=i;                           //computation to be preformed on the value of i
            copy(&exp_e[j], &W[i][j]);
            for(uint32_t jj=0; jj<WINDOW-1; jj++){      //For the bits of i
                if(ti&1){                       //checking Least Significant Bit (LSB)
                    copy(&exp_d[jj], &Q);
                    if(j){                      //for j>0, extra doubling required
                        for(uint32_t k=0; k<(j*E -1); k++)
                            doubl_2(&Q);
                        doubl_3(&Q);
                    }
                    gmuli(Q.t, CURVE_PARAMETER);
                    add_3(&Q, &W[i][j]);
                }
                ti>>=1;
            }
            gmuli(W[i][j].t, CURVE_PARAMETER);
        }
    }
}

//Constant time table look-up - borrowed from ed25519

void fe_cmov(type64 f[],type64 g[],int ib){
  f[0]^=(f[0]^g[0])&(type64)ib;
  f[1]^=(f[1]^g[1])&(type64)ib;
  f[2]^=(f[2]^g[2])&(type64)ib;
  f[3]^=(f[3]^g[3])&(type64)ib;
  f[4]^=(f[4]^g[4])&(type64)ib;
  f[5]^=(f[5]^g[5])&(type64)ib;
  f[6]^=(f[6]^g[6])&(type64)ib;
  f[7]^=(f[7]^g[7])&(type64)ib;
  f[8]^=(f[8]^g[8])&(type64)ib;
}

static void cmov(ECp *w, ECp *u, int b){
  fe_cmov(w->x,u->x,b);
  fe_cmov(w->y,u->y,b);
  fe_cmov(w->z,u->z,b);
  fe_cmov(w->t,u->t,b);
}

// return 1 if b==c, no branching
static int equal(int b, int c){
	int x=b^c;
	x-=1;  // if x=0, x now -1
	return (x>>31);
}

static void ct_select(ECp *T, ECp W[][V], const uint32_t col, const int32_t s, const uint32_t b){
  ECp MT;

  cmov(T,&W[0][col],equal(b,0));  // conditional move
  cmov(T,&W[1][col],equal(b,1));
  cmov(T,&W[2][col],equal(b,2));
  cmov(T,&W[3][col],equal(b,3));
  cmov(T,&W[4][col],equal(b,4));
  cmov(T,&W[5][col],equal(b,5));
  cmov(T,&W[6][col],equal(b,6));
  cmov(T,&W[7][col],equal(b,7));
#if WINDOW>4
  cmov(T,&W[8][col],equal(b,8));
  cmov(T,&W[9][col],equal(b,9));
  cmov(T,&W[10][col],equal(b,10));
  cmov(T,&W[11][col],equal(b,11));
  cmov(T,&W[12][col],equal(b,12));
  cmov(T,&W[13][col],equal(b,13));
  cmov(T,&W[14][col],equal(b,14));
  cmov(T,&W[15][col],equal(b,15));
#if WINDOW>5
  cmov(T,&W[16][col],equal(b,16));
  cmov(T,&W[17][col],equal(b,17));
  cmov(T,&W[18][col],equal(b,18));
  cmov(T,&W[19][col],equal(b,19));
  cmov(T,&W[20][col],equal(b,20));
  cmov(T,&W[21][col],equal(b,21));
  cmov(T,&W[22][col],equal(b,22));
  cmov(T,&W[23][col],equal(b,23));
  cmov(T,&W[24][col],equal(b,24));
  cmov(T,&W[25][col],equal(b,25));
  cmov(T,&W[26][col],equal(b,26));
  cmov(T,&W[27][col],equal(b,27));
  cmov(T,&W[28][col],equal(b,28));
  cmov(T,&W[29][col],equal(b,29));
  cmov(T,&W[30][col],equal(b,30));
  cmov(T,&W[31][col],equal(b,31));
#endif
#endif
  neg(T,&MT);       // minus t
  cmov(T,&MT,s>>1);     //because s={-1, 1}, therefore, we shift only by two bits
}


// Scalar Multiplication - scalar is 519-bit
void scalar_mul(int32_t S[][E], uint32_t T[][E], ECp W[][V], ECp *P){
    int32_t i, j;
    ECp Q;

    ct_select(&Q, W, 0, S[0][E-1], T[0][E-1]);
    //P+=Q where P is the neutral element
    TMV_product(Q.x, Q.z, P->x);
    TMV_product(Q.y, Q.z, P->y);
    gsqr(Q.z, P->z);
    TMV_product(Q.x, Q.y, P->t);

	for(i=1; i<(V-1); i++){
        ct_select(&Q, W, i, S[i][E-1], T[i][E-1]);
        add_3(&Q, P);
	}
	ct_select(&Q, W, i, S[i][E-1], T[i][E-1]);
    add_2(&Q, P);

    for(i=E-2; i>=0; i--){
        doubl_3(P);
        for(j=0; j<(V-1); j++){
            ct_select(&Q, W, j, S[j][i], T[j][i]);
            add_3(&Q, P);
        }
        ct_select(&Q, W, j, S[j][i], T[j][i]);
        add_2(&Q, P);
    }
    norm(P);

}


/////////////////////////////////////////////////////////////////////////////////////////////////

int main(){
	uint64_t bef, aft, min_ccycle=100000000;
	int i, j,lpz=10000;
	type64 xs[LIMBS], ys[LIMBS];
	ECp P, PreT[M][V];

/* Base point on Twisted Edwards Curve (from SafeCurves Website) */

	xs[0]=0x2A940A2F19BA6CLL;           	xs[1]=0x3EC4CD920E2A8CLL;           	xs[2]=0x1D568FC99C6059DLL;
	xs[3]=0x3331C90D2C6BA52LL;          	xs[4]=0xC6203913F6ECC5LL;           	xs[5]=0x1B2063B22FCF270LL;
	xs[6]=0x2878A3BFD9F42FCLL;      	    xs[7]=0x6277E432C8A5ACLL;           	xs[8]=0x752CB45C48648BLL;

	ys[0]=0xcLL;        ys[1]=0x0LL;        ys[2]=0x0LL;
	ys[3]=0x0LL;        ys[4]=0x0LL;        ys[5]=0x0LL;
	ys[6]=0x0LL;    	ys[7]=0x0LL;    	ys[8]=0x0LL;

	/*We have chosen the following at random for testing
    Private key =Scalar=1194302085452638630674850501692388541378855954395238638681997303850335875726918918577735190749226078935456026504296836368946206823467968167314247278775950881
    The mLSB-set of the above private key for window width-4 is
    mLSB-set={-1, -1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1,
    -1, 1, -1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, -1, -1, -1,
    -1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 0, 0, -1, 0, 1, -1, 0, 0, 1, 1, 1, 1, 0, 1, 0, -1, 1, 1, 0, 0, -1,
    0, 0, 0, -1, -1, 0, 1, 0, 0, 0, -1, 0, 1, 0, 1, -1, 1, 0, 0, -1, -1, 0, -1, 1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, 1, 0, 0, 0, 1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, -1,
    0, -1, 1, -1, -1, 0, 1, -1, -1, 0, -1, 1, -1, 0, 1, 0, -1, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, -1, 0, 1, 0, 1, -1, 1, 1, 0, -1, 1, 1, 0, 0, -1, 0, 0, 0, 0,
    -1, -1, -1, 1, -1, -1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, -1, 1, -1, -1, 0, 1, 1, 0, -1, -1, 1, 1, 0, 1, -1, 0, -1, 0, 0, -1, 0, -1, 1, 0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 0,
    -1, 1, 1, -1, -1, 1, -1, 0, -1, -1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, -1, 0, -1, -1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, -1, 0, 0, -1, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, 1, -1, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 1, 0, 1, 1, -1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, -1, -1, 1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, -1, 0, 1, 0, -1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1, 0, 1, -1, 0,
    0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 1, -1, 1, 0, 0, 0, 0, 0, -1, 1, 1, -1, -1, 1, -1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, -1, 1, 0, 0, -1, 1, 0, 0, -1, -1, 1, -1, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 1, 1, 0, 1, -1, 0, 1, 0, 1, 1, 0, 0, -1, 1, 0, 0, 1, 0, 0, -1, 0, 1, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0}
    The mLSB-set of the above private key for window width-5 is
    mLSB-set={-1, -1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1,
    -1, 1, -1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, -1, -1, -1,
    -1, 1, 1, 1, 1, 1, 1, 0, -1, 0, -1, 0, 0, 0, 0, 1, 1, 0, 1, -1, 0, 0, -1, 0, 0, 1, 0, -1, 0, -1, 1, 0, -1, 1, 0, 0, -1, 0, -1, 0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 0, -1, 0, 1, 0, -1,
    0, 0, 0, -1, -1, 0, -1, -1, -1, 1, 1, -1, -1, 1, -1, 0, -1, -1, 1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
    0, 1, 0, 0, -1, -1, 0, -1, 0, -1, 0, 1, 0, 1, -1, 1, 0, -1, 1, 1, 1, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, -1, -1, 1, 0, 1, 1, 0, 1, 0, 0, -1, 0, 1, -1, 1, 1, 0, -1, 1, 0, 1, -1, 0, 0,
    0, -1, -1, 1, 0, -1, -1, 0, 0, 0, -1, 0, 1, -1, 1, 0, 1, 1, -1, 1, 0, 0, -1, -1, 0, 1, -1, -1, 0, 1, -1, -1, -1, 0, 1, -1, 0, 0, -1, 0, 0, -1, 0, -1, -1, 1, 0, 1, 0, 0, 1, 0, -1, 0,
    -1, 0, -1, -1, 0, 0, 1, 0, 1, -1, 1, 0, -1, 0, 1, 0, 0, -1, 1, -1, 0, -1, -1, 1, 1, 1, 0, 0, -1, 1, 0, 0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 1, 1, 0, -1, 0, -1, 1, 0, -1, -1, -1, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, -1, 1, -1, 1, 1, 1, 1, -1, 0, 0, -1, 0, 0, -1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, -1, 0, -1, 1, 1, 1, 1, 0, 1, -1, 0, -1, -1, 1, -1, 0, -1, 0,
    1, 0, 1, -1, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, -1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0, -1, 1, -1, 0, 0, -1, -1, -1, 1, 1, 0, 0, 0, -1, -1, -1,
    0, 0, 0, 0, 1, 1, 1, -1, 0, 0, 0, -1, 0, -1, 0, 0, 0, -1, 1, -1, 0, -1, -1, 0, -1, 1, 0, -1, -1, -1, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0}
    The mLSB-set of the above private key for window width-6 is
    mLSB-set={-1, -1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1,
    1, -1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, 0, 0, -1, -1, 0, 0, -1, -1, 1, 0, 0, 0, 0, 1, 1, -1,
    0, 0, 1, -1, 0, 1, -1, 0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 1, 1, 1, 0, 0, -1, 0, -1, -1, 1, -1, 1, 0, 0, -1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 1, 0, -1, 1, -1, -1, -1, 0, 0, 0, 0, 1, 0, 1, 0, 0,
    0, 0, -1, 0, -1, 1, -1, 0, -1, 1, 0, -1, 1, -1, -1, -1, -1, 1, 0, 0, 0, 1, 0, 1, 1, -1, 0, 1, 0, 0, 0, 0, -1, 0, 1, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 1, 0, 1, 0, 0, -1, 1, -1, -1, 0, -1,
    0, 1, 0, -1, 0, -1, 1, -1, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, -1, 1, 1, 0,
    1, 0, 0, -1, 0, 0, -1, 1, 0, -1, 1, 1, 0, -1, -1, 0, 1, 1, 0, 0, -1, 1, -1, 0, 0, -1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 0, 1, 1, 0, 0, 1, -1, -1, -1, 0, 0, -1, 1, 1, 1, 0, 0, 0,
    1, 0, -1, 0, -1, 1, 0, -1, 0, 1, 0, -1, 1, 0, 0, 0, 0, 1, 0, -1, -1, 0, 0, 0, 1, -1, 1, 1, -1, 1, 0, 0, -1, -1, 1, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 0, 0, 1, -1, 0,
    1, 0, 0, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, -1, 0, -1, 1, -1, 0, 0, 1, 1, -1, 0, 1, 0, 0, -1, 0, 1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0,
    1, 0, -1, 0, 1, 0, 1, -1, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 1, 0, 0, 0, -1, 1, 0, 0, -1, 1, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0, -1, 1, 1, 1, 0, 0, 1, 0, 0,
    0, -1, -1, 1, -1, 0, -1, 0, -1, -1, 1}

    */

    uint32_t TT[V][E];         //TT is the table, T, of the absolute values
    int32_t  SS[V][E];         //SS is the table, s, of (first row) signs

#ifndef TEST
#if WINDOW==6
    TT[0][0]=2; TT[0][1]=2; TT[0][2]=3; TT[0][3]=19; TT[0][4]=14; TT[0][5]=0; TT[0][6]=9; TT[0][7]=25; TT[0][8]=7; TT[0][9]=16; TT[0][10]=22; TT[0][11]=14; TT[0][12]=14; TT[0][13]=13; TT[0][14]=15; TT[0][15]=9; TT[0][16]=12; TT[0][17]=16; TT[0][18]=1; TT[0][19]=31; TT[0][20]=8; TT[0][21]=27; TT[0][22]=13; TT[0][23]=22; TT[0][24]=18; TT[0][25]=28; TT[0][26]=13; TT[0][27]=5; TT[0][28]=16;
    TT[1][0]=5; TT[1][1]=4; TT[1][2]=2; TT[1][3]=14; TT[1][4]=7; TT[1][5]=25; TT[1][6]=3; TT[1][7]=28; TT[1][8]=4; TT[1][9]=7; TT[1][10]=18; TT[1][11]=3; TT[1][12]=7; TT[1][13]=9; TT[1][14]=27; TT[1][15]=21; TT[1][16]=14; TT[1][17]=0; TT[1][18]=19; TT[1][19]=24; TT[1][20]=11; TT[1][21]=30; TT[1][22]=18; TT[1][23]=2; TT[1][24]=5; TT[1][25]=20; TT[1][26]=0; TT[1][27]=0; TT[1][28]=5;
    TT[2][0]=15; TT[2][1]=0; TT[2][2]=9; TT[2][3]=21; TT[2][4]=21; TT[2][5]=29; TT[2][6]=5; TT[2][7]=8; TT[2][8]=10; TT[2][9]=28; TT[2][10]=20; TT[2][11]=21; TT[2][12]=30; TT[2][13]=11; TT[2][14]=8; TT[2][15]=18; TT[2][16]=12; TT[2][17]=0; TT[2][18]=7; TT[2][19]=24; TT[2][20]=23; TT[2][21]=29; TT[2][22]=27; TT[2][23]=12; TT[2][24]=27; TT[2][25]=5; TT[2][26]=18; TT[2][27]=23; TT[2][28]=21;

    SS[0][0]=-1; SS[0][1]=-1; SS[0][2]=-1; SS[0][3]=-1; SS[0][4]=1; SS[0][5]=-1; SS[0][6]=-1; SS[0][7]=-1; SS[0][8]=1; SS[0][9]=1; SS[0][10]=1; SS[0][11]=1; SS[0][12]=-1; SS[0][13]=1; SS[0][14]=1; SS[0][15]=-1; SS[0][16]=1; SS[0][17]=1; SS[0][18]=1; SS[0][19]=-1; SS[0][20]=-1; SS[0][21]=1; SS[0][22]=-1; SS[0][23]=1; SS[0][24]=-1; SS[0][25]=-1; SS[0][26]=1; SS[0][27]=1; SS[0][28]=1;
    SS[1][0]=-1; SS[1][1]=-1; SS[1][2]=-1; SS[1][3]=1; SS[1][4]=1; SS[1][5]=1; SS[1][6]=1; SS[1][7]=-1; SS[1][8]=1; SS[1][9]=-1; SS[1][10]=1; SS[1][11]=-1; SS[1][12]=-1; SS[1][13]=1; SS[1][14]=-1; SS[1][15]=1; SS[1][16]=1; SS[1][17]=-1; SS[1][18]=-1; SS[1][19]=1; SS[1][20]=-1; SS[1][21]=1; SS[1][22]=-1; SS[1][23]=-1; SS[1][24]=-1; SS[1][25]=-1; SS[1][26]=-1; SS[1][27]=-1; SS[1][28]=1;
    SS[2][0]=1; SS[2][1]=-1; SS[2][2]=-1; SS[2][3]=1; SS[2][4]=-1; SS[2][5]=-1; SS[2][6]=-1; SS[2][7]=-1; SS[2][8]=1; SS[2][9]=-1; SS[2][10]=1; SS[2][11]=1; SS[2][12]=1; SS[2][13]=1; SS[2][14]=-1; SS[2][15]=1; SS[2][16]=1; SS[2][17]=-1; SS[2][18]=-1; SS[2][19]=-1; SS[2][20]=-1; SS[2][21]=1; SS[2][22]=-1; SS[2][23]=-1; SS[2][24]=-1; SS[2][25]=1; SS[2][26]=-1; SS[2][27]=-1; SS[2][28]=1;
#endif
#if WINDOW==5
    TT[0][0]=8; TT[0][1]=5; TT[0][2]=10; TT[0][3]=15; TT[0][4]=8; TT[0][5]=14; TT[0][6]=4; TT[0][7]=10; TT[0][8]=1; TT[0][9]=15; TT[0][10]=0; TT[0][11]=15; TT[0][12]=15; TT[0][13]=6; TT[0][14]=0; TT[0][15]=7; TT[0][16]=2; TT[0][17]=6; TT[0][18]=11; TT[0][19]=10; TT[0][20]=5; TT[0][21]=4; TT[0][22]=15; TT[0][23]=3; TT[0][24]=12; TT[0][25]=5; TT[0][26]=13; TT[0][27]=4; TT[0][28]=4; TT[0][29]=1; TT[0][30]=2; TT[0][31]=7; TT[0][32]=14; TT[0][33]=0; TT[0][34]=10;
    TT[1][0]=11; TT[1][1]=1; TT[1][2]=6; TT[1][3]=8; TT[1][4]=5; TT[1][5]=11; TT[1][6]=5; TT[1][7]=2; TT[1][8]=11; TT[1][9]=14; TT[1][10]=7; TT[1][11]=0; TT[1][12]=7; TT[1][13]=2; TT[1][14]=12; TT[1][15]=14; TT[1][16]=11; TT[1][17]=5; TT[1][18]=4; TT[1][19]=13; TT[1][20]=11; TT[1][21]=11; TT[1][22]=11; TT[1][23]=13; TT[1][24]=3; TT[1][25]=3; TT[1][26]=1; TT[1][27]=9; TT[1][28]=8; TT[1][29]=11; TT[1][30]=5; TT[1][31]=7; TT[1][32]=7; TT[1][33]=7; TT[1][34]=13;
    TT[2][0]=14; TT[2][1]=14; TT[2][2]=14; TT[2][3]=2; TT[2][4]=0; TT[2][5]=4; TT[2][6]=10; TT[2][7]=3; TT[2][8]=13; TT[2][9]=6; TT[2][10]=2; TT[2][11]=6; TT[2][12]=8; TT[2][13]=15; TT[2][14]=10; TT[2][15]=2; TT[2][16]=10; TT[2][17]=8; TT[2][18]=3; TT[2][19]=11; TT[2][20]=9; TT[2][21]=4; TT[2][22]=11; TT[2][23]=12; TT[2][24]=8; TT[2][25]=2; TT[2][26]=4; TT[2][27]=10; TT[2][28]=6; TT[2][29]=14; TT[2][30]=5; TT[2][31]=6; TT[2][32]=5; TT[2][33]=0; TT[2][34]=7;

    SS[0][0]=-1; SS[0][1]=-1; SS[0][2]=-1; SS[0][3]=-1; SS[0][4]=1; SS[0][5]=-1; SS[0][6]=-1; SS[0][7]=-1; SS[0][8]=1; SS[0][9]=1; SS[0][10]=1; SS[0][11]=1; SS[0][12]=-1; SS[0][13]=1; SS[0][14]=1; SS[0][15]=-1; SS[0][16]=1; SS[0][17]=1; SS[0][18]=1; SS[0][19]=-1; SS[0][20]=-1; SS[0][21]=1; SS[0][22]=-1; SS[0][23]=1; SS[0][24]=-1; SS[0][25]=-1; SS[0][26]=1; SS[0][27]=1; SS[0][28]=1; SS[0][29]=-1; SS[0][30]=-1; SS[0][31]=-1; SS[0][32]=1; SS[0][33]=1; SS[0][34]=1;
    SS[1][0]=1; SS[1][1]=-1; SS[1][2]=1; SS[1][3]=-1; SS[1][4]=1; SS[1][5]=-1; SS[1][6]=-1; SS[1][7]=1; SS[1][8]=-1; SS[1][9]=1; SS[1][10]=1; SS[1][11]=-1; SS[1][12]=-1; SS[1][13]=1; SS[1][14]=-1; SS[1][15]=1; SS[1][16]=-1; SS[1][17]=-1; SS[1][18]=-1; SS[1][19]=-1; SS[1][20]=-1; SS[1][21]=-1; SS[1][22]=1; SS[1][23]=1; SS[1][24]=-1; SS[1][25]=-1; SS[1][26]=1; SS[1][27]=-1; SS[1][28]=-1; SS[1][29]=-1; SS[1][30]=-1; SS[1][31]=1; SS[1][32]=-1; SS[1][33]=1; SS[1][34]=1;
    SS[2][0]=1; SS[2][1]=1; SS[2][2]=-1; SS[2][3]=1; SS[2][4]=1; SS[2][5]=-1; SS[2][6]=-1; SS[2][7]=-1; SS[2][8]=-1; SS[2][9]=1; SS[2][10]=-1; SS[2][11]=-1; SS[2][12]=-1; SS[2][13]=1; SS[2][14]=-1; SS[2][15]=-1; SS[2][16]=-1; SS[2][17]=-1; SS[2][18]=1; SS[2][19]=-1; SS[2][20]=1; SS[2][21]=1; SS[2][22]=-1; SS[2][23]=-1; SS[2][24]=-1; SS[2][25]=-1; SS[2][26]=-1; SS[2][27]=-1; SS[2][28]=-1; SS[2][29]=1; SS[2][30]=1; SS[2][31]=1; SS[2][32]=1; SS[2][33]=1; SS[2][34]=1;
#endif
#if WINDOW==4
    TT[0][0]=4; TT[0][1]=2; TT[0][2]=7; TT[0][3]=6; TT[0][4]=7; TT[0][5]=3; TT[0][6]=2; TT[0][7]=0; TT[0][8]=1; TT[0][9]=1; TT[0][10]=7; TT[0][11]=3; TT[0][12]=4; TT[0][13]=3; TT[0][14]=0; TT[0][15]=5; TT[0][16]=3; TT[0][17]=7; TT[0][18]=2; TT[0][19]=4; TT[0][20]=1; TT[0][21]=0; TT[0][22]=6; TT[0][23]=2; TT[0][24]=3; TT[0][25]=3; TT[0][26]=4; TT[0][27]=3; TT[0][28]=2; TT[0][29]=4; TT[0][30]=2; TT[0][31]=7; TT[0][32]=2; TT[0][33]=7; TT[0][34]=0; TT[0][35]=7; TT[0][36]=7; TT[0][37]=1; TT[0][38]=2; TT[0][39]=0; TT[0][40]=1; TT[0][41]=3; TT[0][42]=4; TT[0][43]=7;
    TT[1][0]=3; TT[1][1]=0; TT[1][2]=4; TT[1][3]=2; TT[1][4]=6; TT[1][5]=4; TT[1][6]=7; TT[1][7]=0; TT[1][8]=3; TT[1][9]=0; TT[1][10]=2; TT[1][11]=0; TT[1][12]=7; TT[1][13]=7; TT[1][14]=6; TT[1][15]=6; TT[1][16]=6; TT[1][17]=7; TT[1][18]=7; TT[1][19]=1; TT[1][20]=7; TT[1][21]=2; TT[1][22]=0; TT[1][23]=4; TT[1][24]=6; TT[1][25]=0; TT[1][26]=0; TT[1][27]=2; TT[1][28]=0; TT[1][29]=7; TT[1][30]=2; TT[1][31]=5; TT[1][32]=1; TT[1][33]=0; TT[1][34]=7; TT[1][35]=5; TT[1][36]=3; TT[1][37]=3; TT[1][38]=4; TT[1][39]=7; TT[1][40]=1; TT[1][41]=1; TT[1][42]=4; TT[1][43]=5;
    TT[2][0]=7; TT[2][1]=5; TT[2][2]=2; TT[2][3]=1; TT[2][4]=2; TT[2][5]=5; TT[2][6]=1; TT[2][7]=2; TT[2][8]=0; TT[2][9]=2; TT[2][10]=0; TT[2][11]=0; TT[2][12]=7; TT[2][13]=4; TT[2][14]=0; TT[2][15]=6; TT[2][16]=4; TT[2][17]=0; TT[2][18]=6; TT[2][19]=2; TT[2][20]=5; TT[2][21]=6; TT[2][22]=0; TT[2][23]=0; TT[2][24]=5; TT[2][25]=4; TT[2][26]=1; TT[2][27]=0; TT[2][28]=7; TT[2][29]=0; TT[2][30]=3; TT[2][31]=5; TT[2][32]=3; TT[2][33]=7; TT[2][34]=6; TT[2][35]=5; TT[2][36]=5; TT[2][37]=1; TT[2][38]=2; TT[2][39]=0; TT[2][40]=1; TT[2][41]=2; TT[2][42]=0; TT[2][43]=0;

    SS[0][0]=-1; SS[0][1]=-1; SS[0][2]=-1; SS[0][3]=-1; SS[0][4]=1; SS[0][5]=-1; SS[0][6]=-1; SS[0][7]=-1; SS[0][8]=1; SS[0][9]=1; SS[0][10]=1; SS[0][11]=1; SS[0][12]=-1; SS[0][13]=1; SS[0][14]=1; SS[0][15]=-1; SS[0][16]=1; SS[0][17]=1; SS[0][18]=1; SS[0][19]=-1; SS[0][20]=-1; SS[0][21]=1; SS[0][22]=-1; SS[0][23]=1; SS[0][24]=-1; SS[0][25]=-1; SS[0][26]=1; SS[0][27]=1; SS[0][28]=1; SS[0][29]=-1; SS[0][30]=-1; SS[0][31]=-1; SS[0][32]=1; SS[0][33]=1; SS[0][34]=1; SS[0][35]=1; SS[0][36]=-1; SS[0][37]=1; SS[0][38]=-1; SS[0][39]=1; SS[0][40]=-1; SS[0][41]=-1; SS[0][42]=1; SS[0][43]=-1;
    SS[1][0]=1; SS[1][1]=1; SS[1][2]=-1; SS[1][3]=-1; SS[1][4]=1; SS[1][5]=-1; SS[1][6]=1; SS[1][7]=-1; SS[1][8]=-1; SS[1][9]=-1; SS[1][10]=-1; SS[1][11]=-1; SS[1][12]=-1; SS[1][13]=1; SS[1][14]=1; SS[1][15]=-1; SS[1][16]=-1; SS[1][17]=1; SS[1][18]=-1; SS[1][19]=-1; SS[1][20]=-1; SS[1][21]=-1; SS[1][22]=1; SS[1][23]=-1; SS[1][24]=1; SS[1][25]=1; SS[1][26]=1; SS[1][27]=1; SS[1][28]=-1; SS[1][29]=1; SS[1][30]=1; SS[1][31]=-1; SS[1][32]=-1; SS[1][33]=-1; SS[1][34]=-1; SS[1][35]=1; SS[1][36]=-1; SS[1][37]=-1; SS[1][38]=-1; SS[1][39]=1; SS[1][40]=-1; SS[1][41]=-1; SS[1][42]=-1; SS[1][43]=-1;
    SS[2][0]=1; SS[2][1]=-1; SS[2][2]=1; SS[2][3]=1; SS[2][4]=-1; SS[2][5]=-1; SS[2][6]=-1; SS[2][7]=-1; SS[2][8]=-1; SS[2][9]=-1; SS[2][10]=-1; SS[2][11]=1; SS[2][12]=1; SS[2][13]=1; SS[2][14]=1; SS[2][15]=1; SS[2][16]=-1; SS[2][17]=1; SS[2][18]=1; SS[2][19]=-1; SS[2][20]=1; SS[2][21]=1; SS[2][22]=1; SS[2][23]=1; SS[2][24]=-1; SS[2][25]=1; SS[2][26]=-1; SS[2][27]=1; SS[2][28]=1; SS[2][29]=1; SS[2][30]=1; SS[2][31]=-1; SS[2][32]=1; SS[2][33]=1; SS[2][34]=-1; SS[2][35]=-1; SS[2][36]=1; SS[2][37]=1; SS[2][38]=-1; SS[2][39]=-1; SS[2][40]=-1; SS[2][41]=1; SS[2][42]=-1; SS[2][43]=1;
#endif

#else
#if WINDOW==6
    TT[0][0]=3; TT[0][1]=2; TT[0][2]=3; TT[0][3]=0; TT[0][4]=1; TT[0][5]=3; TT[0][6]=2; TT[0][7]=0; TT[0][8]=3; TT[0][9]=1; TT[0][10]=3; TT[0][11]=3; TT[0][12]=1; TT[0][13]=2; TT[0][14]=1; TT[0][15]=0; TT[0][16]=3; TT[0][17]=3; TT[0][18]=1; TT[0][19]=1; TT[0][20]=0; TT[0][21]=2; TT[0][22]=0; TT[0][23]=1; TT[0][24]=1; TT[0][25]=1; TT[0][26]=0; TT[0][27]=2; TT[0][28]=3;
    TT[1][0]=1; TT[1][1]=0; TT[1][2]=1; TT[1][3]=1; TT[1][4]=2; TT[1][5]=0; TT[1][6]=1; TT[1][7]=3; TT[1][8]=3; TT[1][9]=0; TT[1][10]=0; TT[1][11]=0; TT[1][12]=1; TT[1][13]=2; TT[1][14]=1; TT[1][15]=1; TT[1][16]=2; TT[1][17]=2; TT[1][18]=2; TT[1][19]=0; TT[1][20]=3; TT[1][21]=0; TT[1][22]=0; TT[1][23]=0; TT[1][24]=2; TT[1][25]=0; TT[1][26]=1; TT[1][27]=3; TT[1][28]=2;
    TT[2][0]=3; TT[2][1]=0; TT[2][2]=1; TT[2][3]=3; TT[2][4]=0; TT[2][5]=1; TT[2][6]=3; TT[2][7]=1; TT[2][8]=3; TT[2][9]=0; TT[2][10]=0; TT[2][11]=2; TT[2][12]=2; TT[2][13]=0; TT[2][14]=3; TT[2][15]=2; TT[2][16]=3; TT[2][17]=1; TT[2][18]=3; TT[2][19]=2; TT[2][20]=1; TT[2][21]=2; TT[2][22]=1; TT[2][23]=1; TT[2][24]=3; TT[2][25]=1; TT[2][26]=18; TT[2][27]=17; TT[2][28]=16;

    SS[0][0]=1; SS[0][1]=-1; SS[0][2]=1; SS[0][3]=-1; SS[0][4]=1; SS[0][5]=1; SS[0][6]=-1; SS[0][7]=1; SS[0][8]=-1; SS[0][9]=1; SS[0][10]=1; SS[0][11]=-1; SS[0][12]=-1; SS[0][13]=-1; SS[0][14]=-1; SS[0][15]=-1; SS[0][16]=-1; SS[0][17]=-1; SS[0][18]=1; SS[0][19]=1; SS[0][20]=-1; SS[0][21]=-1; SS[0][22]=-1; SS[0][23]=1; SS[0][24]=-1; SS[0][25]=1; SS[0][26]=-1; SS[0][27]=1; SS[0][28]=1;
    SS[1][0]=1; SS[1][1]=1; SS[1][2]=1; SS[1][3]=-1; SS[1][4]=1; SS[1][5]=-1; SS[1][6]=1; SS[1][7]=1; SS[1][8]=-1; SS[1][9]=-1; SS[1][10]=-1; SS[1][11]=-1; SS[1][12]=1; SS[1][13]=-1; SS[1][14]=-1; SS[1][15]=1; SS[1][16]=-1; SS[1][17]=-1; SS[1][18]=-1; SS[1][19]=1; SS[1][20]=-1; SS[1][21]=1; SS[1][22]=-1; SS[1][23]=1; SS[1][24]=1; SS[1][25]=1; SS[1][26]=-1; SS[1][27]=-1; SS[1][28]=-1;
    SS[2][0]=-1; SS[2][1]=-1; SS[2][2]=-1; SS[2][3]=1; SS[2][4]=-1; SS[2][5]=1; SS[2][6]=-1; SS[2][7]=1; SS[2][8]=-1; SS[2][9]=-1; SS[2][10]=-1; SS[2][11]=1; SS[2][12]=-1; SS[2][13]=1; SS[2][14]=1; SS[2][15]=1; SS[2][16]=1; SS[2][17]=1; SS[2][18]=-1; SS[2][19]=-1; SS[2][20]=-1; SS[2][21]=1; SS[2][22]=1; SS[2][23]=1; SS[2][24]=1; SS[2][25]=-1; SS[2][26]=-1; SS[2][27]=-1; SS[2][28]=1;
#endif
#if WINDOW==5
    TT[0][0]=2; TT[0][1]=3; TT[0][2]=1; TT[0][3]=0; TT[0][4]=0; TT[0][5]=1; TT[0][6]=3; TT[0][7]=1; TT[0][8]=0; TT[0][9]=2; TT[0][10]=1; TT[0][11]=1; TT[0][12]=3; TT[0][13]=2; TT[0][14]=0; TT[0][15]=0; TT[0][16]=0; TT[0][17]=3; TT[0][18]=2; TT[0][19]=0; TT[0][20]=2; TT[0][21]=2; TT[0][22]=2; TT[0][23]=1; TT[0][24]=0; TT[0][25]=3; TT[0][26]=0; TT[0][27]=0; TT[0][28]=2; TT[0][29]=2; TT[0][30]=0; TT[0][31]=3; TT[0][32]=3; TT[0][33]=2; TT[0][34]=2;
    TT[1][0]=0; TT[1][1]=2; TT[1][2]=3; TT[1][3]=1; TT[1][4]=2; TT[1][5]=1; TT[1][6]=2; TT[1][7]=3; TT[1][8]=3; TT[1][9]=3; TT[1][10]=1; TT[1][11]=2; TT[1][12]=3; TT[1][13]=0; TT[1][14]=0; TT[1][15]=0; TT[1][16]=0; TT[1][17]=0; TT[1][18]=0; TT[1][19]=1; TT[1][20]=0; TT[1][21]=1; TT[1][22]=0; TT[1][23]=0; TT[1][24]=0; TT[1][25]=1; TT[1][26]=0; TT[1][27]=1; TT[1][28]=0; TT[1][29]=0; TT[1][30]=0; TT[1][31]=0; TT[1][32]=1; TT[1][33]=0; TT[1][34]=1;
    TT[2][0]=1; TT[2][1]=1; TT[2][2]=0; TT[2][3]=0; TT[2][4]=1; TT[2][5]=1; TT[2][6]=1; TT[2][7]=1; TT[2][8]=0; TT[2][9]=1; TT[2][10]=1; TT[2][11]=1; TT[2][12]=0; TT[2][13]=1; TT[2][14]=0; TT[2][15]=1; TT[2][16]=1; TT[2][17]=1; TT[2][18]=0; TT[2][19]=0; TT[2][20]=1; TT[2][21]=1; TT[2][22]=1; TT[2][23]=0; TT[2][24]=0; TT[2][25]=0; TT[2][26]=1; TT[2][27]=1; TT[2][28]=0; TT[2][29]=8; TT[2][30]=0; TT[2][31]=0; TT[2][32]=1; TT[2][33]=1; TT[2][34]=1;

    SS[0][0]=1; SS[0][1]=-1; SS[0][2]=1; SS[0][3]=-1; SS[0][4]=1; SS[0][5]=1; SS[0][6]=-1; SS[0][7]=1; SS[0][8]=-1; SS[0][9]=1; SS[0][10]=1; SS[0][11]=-1; SS[0][12]=-1; SS[0][13]=-1; SS[0][14]=-1; SS[0][15]=-1; SS[0][16]=-1; SS[0][17]=-1; SS[0][18]=1; SS[0][19]=1; SS[0][20]=-1; SS[0][21]=-1; SS[0][22]=-1; SS[0][23]=1; SS[0][24]=-1; SS[0][25]=1; SS[0][26]=-1; SS[0][27]=1; SS[0][28]=1; SS[0][29]=1; SS[0][30]=1; SS[0][31]=1; SS[0][32]=-1; SS[0][33]=1; SS[0][34]=-1;
    SS[1][0]=1; SS[1][1]=1; SS[1][2]=-1; SS[1][3]=-1; SS[1][4]=-1; SS[1][5]=-1; SS[1][6]=1; SS[1][7]=-1; SS[1][8]=-1; SS[1][9]=1; SS[1][10]=-1; SS[1][11]=-1; SS[1][12]=-1; SS[1][13]=1; SS[1][14]=-1; SS[1][15]=1; SS[1][16]=-1; SS[1][17]=1; SS[1][18]=1; SS[1][19]=1; SS[1][20]=-1; SS[1][21]=-1; SS[1][22]=-1; SS[1][23]=-1; SS[1][24]=-1; SS[1][25]=-1; SS[1][26]=1; SS[1][27]=-1; SS[1][28]=1; SS[1][29]=-1; SS[1][30]=1; SS[1][31]=-1; SS[1][32]=-1; SS[1][33]=-1; SS[1][34]=1;
    SS[2][0]=-1; SS[2][1]=1; SS[2][2]=1; SS[2][3]=1; SS[2][4]=1; SS[2][5]=1; SS[2][6]=-1; SS[2][7]=-1; SS[2][8]=-1; SS[2][9]=1; SS[2][10]=1; SS[2][11]=1; SS[2][12]=1; SS[2][13]=-1; SS[2][14]=-1; SS[2][15]=-1; SS[2][16]=1; SS[2][17]=-1; SS[2][18]=1; SS[2][19]=-1; SS[2][20]=1; SS[2][21]=1; SS[2][22]=-1; SS[2][23]=-1; SS[2][24]=1; SS[2][25]=-1; SS[2][26]=1; SS[2][27]=1; SS[2][28]=-1; SS[2][29]=1; SS[2][30]=-1; SS[2][31]=1; SS[2][32]=-1; SS[2][33]=-1; SS[2][34]=1;
#endif
#if WINDOW==4
    TT[0][0]=0; TT[0][1]=0; TT[0][2]=0; TT[0][3]=0; TT[0][4]=1; TT[0][5]=1; TT[0][6]=1; TT[0][7]=0; TT[0][8]=0; TT[0][9]=0; TT[0][10]=1; TT[0][11]=0; TT[0][12]=1; TT[0][13]=1; TT[0][14]=0; TT[0][15]=1; TT[0][16]=1; TT[0][17]=1; TT[0][18]=0; TT[0][19]=0; TT[0][20]=1; TT[0][21]=0; TT[0][22]=0; TT[0][23]=0; TT[0][24]=0; TT[0][25]=0; TT[0][26]=0; TT[0][27]=1; TT[0][28]=0; TT[0][29]=1; TT[0][30]=1; TT[0][31]=1; TT[0][32]=1; TT[0][33]=1; TT[0][34]=1; TT[0][35]=1; TT[0][36]=1; TT[0][37]=1; TT[0][38]=0; TT[0][39]=0; TT[0][40]=1; TT[0][41]=0; TT[0][42]=1; TT[0][43]=0;
    TT[1][0]=1; TT[1][1]=0; TT[1][2]=0; TT[1][3]=1; TT[1][4]=0; TT[1][5]=0; TT[1][6]=1; TT[1][7]=1; TT[1][8]=1; TT[1][9]=1; TT[1][10]=1; TT[1][11]=0; TT[1][12]=1; TT[1][13]=0; TT[1][14]=1; TT[1][15]=1; TT[1][16]=0; TT[1][17]=0; TT[1][18]=0; TT[1][19]=1; TT[1][20]=1; TT[1][21]=0; TT[1][22]=0; TT[1][23]=0; TT[1][24]=0; TT[1][25]=1; TT[1][26]=1; TT[1][27]=1; TT[1][28]=0; TT[1][29]=0; TT[1][30]=0; TT[1][31]=1; TT[1][32]=0; TT[1][33]=0; TT[1][34]=1; TT[1][35]=0; TT[1][36]=0; TT[1][37]=0; TT[1][38]=0; TT[1][39]=0; TT[1][40]=1; TT[1][41]=0; TT[1][42]=0; TT[1][43]=1;
    TT[2][0]=1; TT[2][1]=0; TT[2][2]=1; TT[2][3]=0; TT[2][4]=1; TT[2][5]=0; TT[2][6]=0; TT[2][7]=1; TT[2][8]=1; TT[2][9]=0; TT[2][10]=1; TT[2][11]=1; TT[2][12]=0; TT[2][13]=1; TT[2][14]=1; TT[2][15]=1; TT[2][16]=1; TT[2][17]=1; TT[2][18]=1; TT[2][19]=0; TT[2][20]=1; TT[2][21]=0; TT[2][22]=0; TT[2][23]=1; TT[2][24]=0; TT[2][25]=0; TT[2][26]=1; TT[2][27]=1; TT[2][28]=1; TT[2][29]=1; TT[2][30]=1; TT[2][31]=0; TT[2][32]=1; TT[2][33]=0; TT[2][34]=0; TT[2][35]=4; TT[2][36]=1; TT[2][37]=0; TT[2][38]=1; TT[2][39]=1; TT[2][40]=0; TT[2][41]=0; TT[2][42]=0; TT[2][43]=0;

    SS[0][0]=1; SS[0][1]=-1; SS[0][2]=1; SS[0][3]=-1; SS[0][4]=1; SS[0][5]=1; SS[0][6]=-1; SS[0][7]=1; SS[0][8]=-1; SS[0][9]=1; SS[0][10]=1; SS[0][11]=-1; SS[0][12]=-1; SS[0][13]=-1; SS[0][14]=-1; SS[0][15]=-1; SS[0][16]=-1; SS[0][17]=-1; SS[0][18]=1; SS[0][19]=1; SS[0][20]=-1; SS[0][21]=-1; SS[0][22]=-1; SS[0][23]=1; SS[0][24]=-1; SS[0][25]=1; SS[0][26]=-1; SS[0][27]=1; SS[0][28]=1; SS[0][29]=1; SS[0][30]=1; SS[0][31]=1; SS[0][32]=-1; SS[0][33]=1; SS[0][34]=-1; SS[0][35]=1; SS[0][36]=1; SS[0][37]=-1; SS[0][38]=-1; SS[0][39]=-1; SS[0][40]=-1; SS[0][41]=1; SS[0][42]=-1; SS[0][43]=-1;
    SS[1][0]=1; SS[1][1]=-1; SS[1][2]=-1; SS[1][3]=-1; SS[1][4]=1; SS[1][5]=-1; SS[1][6]=1; SS[1][7]=-1; SS[1][8]=1; SS[1][9]=1; SS[1][10]=1; SS[1][11]=-1; SS[1][12]=-1; SS[1][13]=-1; SS[1][14]=-1; SS[1][15]=-1; SS[1][16]=-1; SS[1][17]=1; SS[1][18]=-1; SS[1][19]=1; SS[1][20]=-1; SS[1][21]=1; SS[1][22]=-1; SS[1][23]=-1; SS[1][24]=-1; SS[1][25]=1; SS[1][26]=-1; SS[1][27]=1; SS[1][28]=1; SS[1][29]=1; SS[1][30]=1; SS[1][31]=1; SS[1][32]=-1; SS[1][33]=-1; SS[1][34]=-1; SS[1][35]=1; SS[1][36]=1; SS[1][37]=1; SS[1][38]=1; SS[1][39]=-1; SS[1][40]=-1; SS[1][41]=-1; SS[1][42]=1; SS[1][43]=-1;
    SS[2][0]=1; SS[2][1]=-1; SS[2][2]=1; SS[2][3]=1; SS[2][4]=-1; SS[2][5]=-1; SS[2][6]=1; SS[2][7]=-1; SS[2][8]=1; SS[2][9]=1; SS[2][10]=-1; SS[2][11]=1; SS[2][12]=-1; SS[2][13]=1; SS[2][14]=-1; SS[2][15]=-1; SS[2][16]=-1; SS[2][17]=1; SS[2][18]=-1; SS[2][19]=-1; SS[2][20]=-1; SS[2][21]=1; SS[2][22]=1; SS[2][23]=-1; SS[2][24]=-1; SS[2][25]=-1; SS[2][26]=1; SS[2][27]=1; SS[2][28]=-1; SS[2][29]=1; SS[2][30]=1; SS[2][31]=1; SS[2][32]=1; SS[2][33]=-1; SS[2][34]=1; SS[2][35]=1; SS[2][36]=1; SS[2][37]=1; SS[2][38]=1; SS[2][39]=-1; SS[2][40]=-1; SS[2][41]=1; SS[2][42]=-1; SS[2][43]=1;
#endif

#endif

    init(xs, ys, &P);
    Offline_precomp(&P, PreT);

	for(i=0;i<10;i++){
		bef=rdtsc();
		for (j=0; j<lpz; j++){
			inf(&P);               //Initializing as neutral element
			scalar_mul(SS, TT, PreT, &P);
		}
		aft=rdtscp();
		if((aft-bef)/(lpz) < min_ccycle)
				min_ccycle= (aft-bef)/(lpz);
	}
	printf("Window width :: %u\n%s%lu\n", WINDOW, "The minimum clock cycles count is :: ", min_ccycle);
	output(&P);
}	
