/* Test program for ed521 scalar point multiplication
 Uses Bernstein et al.s point multiplication method from Curve41417 paper
 Cache safety thanks to ed25519
 Fully Tested and debugged
 We have replaced gmul() by TMV_product()
 We have added rdtscp() to measure the clock cycles count as suggested by Paoloni; http://www.intel.com.tr/content/dam/www/public/us/en/documents/white-papers/ia-32-ia-64-benchmark-code-execution-paper.pdf
 We have amended scr(), gmuli(), gsqr(), and gsqr2() according to modulus p
 gcc -Wall -O3 vb64_ed521.c -o vb64_ed521.exe */

#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>

#define WINDOW 4  //5 //6

#if WINDOW==4
#define PANES 131
#endif

#if WINDOW==5
#define PANES 105
#endif

#if WINDOW==6
#define PANES 87
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
void gdec(const type64 x[], type64 w[]){
	w[0]-=x[0];
	w[1]-=x[1];
	w[2]-=x[2];
	w[3]-=x[3];
	w[4]-=x[4];
	w[5]-=x[5];
	w[6]-=x[6];
	w[7]-=x[7];
	w[8]-=x[8];
}

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

// w*=2
void gmul2(type64 w[]){
	w[0]*=2;
	w[1]*=2;
	w[2]*=2;
	w[3]*=2;
	w[4]*=2;
	w[5]*=2;
	w[6]*=2;
	w[7]*=2;
	w[8]*=2;
}

// w-=2*x
void gsb2(const type64 x[], type64 w[]){
	w[0]-=2*x[0];
	w[1]-=2*x[1];
	w[2]-=2*x[2];
	w[3]-=2*x[3];
	w[4]-=2*x[4];
	w[5]-=2*x[5];
	w[6]-=2*x[6];
	w[7]-=2*x[7];
	w[8]-=2*x[8];
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


// P=2P, where P->z=1
void double_1(ECp *P){
	type64 a[LIMBS], b[LIMBS], e[LIMBS], f[LIMBS], g[LIMBS], h[LIMBS];
	gsqr(P->x,a);
	gsqr(P->y,b);
	gcopy(P->t,e);
	gmul2(e);
	gadd(a,b,g);
	gcopy(g,f); f[0]-=2;
	gsub(a,b,h);
	TMV_product(e,f,P->x);
	TMV_product(g,h,P->y);
	gsqr(g,P->z);
	gsb2(g,P->z);
	TMV_product(e,h,P->t);
	scr(P->z);
}

// P=2P, when P->t is not required
void double_2(ECp *P){
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
void double_3(ECp *P){
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

//P+=Q, when P->t is required
void add_1(ECp *Q,ECp *P){
	type64 a[LIMBS], b[LIMBS], c[LIMBS], d[LIMBS], e[LIMBS], f[LIMBS], g[LIMBS], h[LIMBS];
	TMV_product(P->x,Q->x,a);
	TMV_product(P->y,Q->y,b);
	TMV_product(P->t,Q->t,c);
	gadd(P->z,c,f);  /* reversed sign as d is negative */
	gsub(P->z,c,g);
	gsub(b,a,h);
	gadd(P->x,P->y,c); gadd(Q->x,Q->y,d); TMV_product(c,d,e); gtsb(a,b,e); //gdec(a,e); gdec(b,e);
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

//P=0, Setting the point as neutral element
void inf(ECp *P){
	for (int i=0; i<LIMBS; i++)
		P->x[i]=P->y[i]=P->z[i]=P->t[i]=0;
	P->y[0]=P->z[0]=1;
}

// Initialise P
void init(type64 *x,type64 *y,ECp *P){
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


// Pre-computation
void precomp(ECp *P,ECp W[]){
	inf(&W[0]);
	copy(P,&W[1]); gmuli(W[1].t,CURVE_PARAMETER);
	copy(P,&W[2]); double_1(&W[2]);
	copy(&W[2],&W[3]); add_1(&W[1],&W[3]);
	copy(&W[2],&W[4]); double_3(&W[4]);
	copy(&W[4],&W[5]); add_1(&W[1],&W[5]);
	copy(&W[3],&W[6]); double_3(&W[6]);
	copy(&W[6],&W[7]); add_1(&W[1],&W[7]);
	copy(&W[4],&W[8]); double_3(&W[8]);

#if WINDOW>4
	copy(&W[8],&W[9]); add_1(&W[1],&W[9]);
	copy(&W[5],&W[10]); double_3(&W[10]);
	copy(&W[10],&W[11]); add_1(&W[1],&W[11]);
	copy(&W[6],&W[12]); double_3(&W[12]);
	copy(&W[12],&W[13]); add_1(&W[1],&W[13]);
	copy(&W[7],&W[14]); double_3(&W[14]);
	copy(&W[14],&W[15]); add_1(&W[1],&W[15]);
	copy(&W[8],&W[16]); double_3(&W[16]);

#if WINDOW>5
	copy(&W[16],&W[17]); add_1(&W[1],&W[17]);
	copy(&W[9],&W[18]); double_3(&W[18]);
	copy(&W[18],&W[19]); add_1(&W[1],&W[19]);
	copy(&W[10],&W[20]); double_3(&W[20]);
	copy(&W[20],&W[21]); add_1(&W[1],&W[21]);
	copy(&W[11],&W[22]); double_3(&W[22]);
	copy(&W[22],&W[23]); add_1(&W[1],&W[23]);
	copy(&W[12],&W[24]); double_3(&W[24]);
	copy(&W[24],&W[25]); add_1(&W[1],&W[25]);
	copy(&W[13],&W[26]); double_3(&W[26]);
	copy(&W[26],&W[27]); add_1(&W[1],&W[27]);
	copy(&W[14],&W[28]); double_3(&W[28]);
	copy(&W[28],&W[29]); add_1(&W[1],&W[29]);
	copy(&W[15],&W[30]); double_3(&W[30]);
	copy(&W[30],&W[31]); add_1(&W[1],&W[31]);
	copy(&W[16],&W[32]); double_3(&W[32]);
#endif
#endif
/* pre-multiply t parameter by curve constant/parameter */

	gmuli(W[2].t,CURVE_PARAMETER);
	gmuli(W[3].t,CURVE_PARAMETER);
	gmuli(W[4].t,CURVE_PARAMETER);
	gmuli(W[5].t,CURVE_PARAMETER);
	gmuli(W[6].t,CURVE_PARAMETER);
	gmuli(W[7].t,CURVE_PARAMETER);
	gmuli(W[8].t,CURVE_PARAMETER);
#if WINDOW>4
	gmuli(W[9].t,CURVE_PARAMETER);
	gmuli(W[10].t,CURVE_PARAMETER);
	gmuli(W[11].t,CURVE_PARAMETER);
	gmuli(W[12].t,CURVE_PARAMETER);
	gmuli(W[13].t,CURVE_PARAMETER);
	gmuli(W[14].t,CURVE_PARAMETER);
	gmuli(W[15].t,CURVE_PARAMETER);
	gmuli(W[16].t,CURVE_PARAMETER);
#if WINDOW>5
	gmuli(W[17].t,CURVE_PARAMETER);
	gmuli(W[18].t,CURVE_PARAMETER);
	gmuli(W[19].t,CURVE_PARAMETER);
	gmuli(W[20].t,CURVE_PARAMETER);
	gmuli(W[21].t,CURVE_PARAMETER);
	gmuli(W[22].t,CURVE_PARAMETER);
	gmuli(W[23].t,CURVE_PARAMETER);
	gmuli(W[24].t,CURVE_PARAMETER);
	gmuli(W[25].t,CURVE_PARAMETER);
	gmuli(W[26].t,CURVE_PARAMETER);
	gmuli(W[27].t,CURVE_PARAMETER);
	gmuli(W[28].t,CURVE_PARAMETER);
	gmuli(W[29].t,CURVE_PARAMETER);
	gmuli(W[30].t,CURVE_PARAMETER);
	gmuli(W[31].t,CURVE_PARAMETER);
	gmuli(W[32].t,CURVE_PARAMETER);
#endif
#endif
}

// Windows of width 4-6

void window(ECp *Q,ECp *P){
	double_2(P);
	double_2(P);
	double_2(P);
#if WINDOW>4
	double_2(P);
#if WINDOW>5
	double_2(P);
#endif
#endif
	double_3(P);
	add_2(Q,P);
}

// Printing the Affine coordinates of the computed new point on the screen
void output(ECp *P){
    int32_t i;
    for(i=0; i<LIMBS; i++)
        printf("x[%d]= %lx\n", i, P->x[i]);

	puts("");
	for(i=0; i<LIMBS; i++)
        printf("y[%d]= %lx\n", i, P->y[i]);

    puts("");
}

//Constant time table look-up - borrowed from ed25519

void fe_cmov(type64 f[],type64 g[],int ib){
  type64 b=-ib;
  f[0]^=(f[0]^g[0])&b;
  f[1]^=(f[1]^g[1])&b;
  f[2]^=(f[2]^g[2])&b;
  f[3]^=(f[3]^g[3])&b;
  f[4]^=(f[4]^g[4])&b;
  f[5]^=(f[5]^g[5])&b;
  f[6]^=(f[6]^g[6])&b;
  f[7]^=(f[7]^g[7])&b;
  f[8]^=(f[8]^g[8])&b;
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
	return ((x>>31)&1);
}

static void ct_select(ECp *T, ECp W[], int b){
  ECp MT;
  int m=b>>31;
  int babs=(b^m)-m;

  cmov(T,&W[0],equal(babs,0));  // conditional move
  cmov(T,&W[1],equal(babs,1));
  cmov(T,&W[2],equal(babs,2));
  cmov(T,&W[3],equal(babs,3));
  cmov(T,&W[4],equal(babs,4));
  cmov(T,&W[5],equal(babs,5));
  cmov(T,&W[6],equal(babs,6));
  cmov(T,&W[7],equal(babs,7));
  cmov(T,&W[8],equal(babs,8));
#if WINDOW>4
  cmov(T,&W[9],equal(babs,9));
  cmov(T,&W[10],equal(babs,10));
  cmov(T,&W[11],equal(babs,11));
  cmov(T,&W[12],equal(babs,12));
  cmov(T,&W[13],equal(babs,13));
  cmov(T,&W[14],equal(babs,14));
  cmov(T,&W[15],equal(babs,15));
  cmov(T,&W[16],equal(babs,16));
#if WINDOW>5
  cmov(T,&W[17],equal(babs,17));
  cmov(T,&W[18],equal(babs,18));
  cmov(T,&W[19],equal(babs,19));
  cmov(T,&W[20],equal(babs,20));
  cmov(T,&W[21],equal(babs,21));
  cmov(T,&W[22],equal(babs,22));
  cmov(T,&W[23],equal(babs,23));
  cmov(T,&W[24],equal(babs,24));
  cmov(T,&W[25],equal(babs,25));
  cmov(T,&W[26],equal(babs,26));
  cmov(T,&W[27],equal(babs,27));
  cmov(T,&W[28],equal(babs,28));
  cmov(T,&W[29],equal(babs,29));
  cmov(T,&W[30],equal(babs,30));
  cmov(T,&W[31],equal(babs,31));
  cmov(T,&W[32],equal(babs,32));
#endif
#endif
  neg(T,&MT);  // minus t
  cmov(T,&MT,m&1);
}


// Scalar/Point Multiplication - scalar is 519-bit
void scalar_mul(int *w, ECp *P){
	ECp W[M+1], Q;

	precomp(P, W);
	copy(&W[w[PANES-1]], P);
	for (int i=PANES-2; i>=0; i--){
		ct_select(&Q, W, w[i]);
		window(&Q, P);
	}
	norm(P);
}


///////////////////////////////////////////////////////////////////////////////////
int main(){
	uint64_t bef, aft, min_ccycle=100000000;
	int w[PANES], i, j,lpz=10000;
	type64 xs[LIMBS], ys[LIMBS];
	ECp P;

/* Base point on Twisted Edwards Curve (from SafeCurves Website) */

	xs[0]=0x2A940A2F19BA6CLL;           	xs[1]=0x3EC4CD920E2A8CLL;           	xs[2]=0x1D568FC99C6059DLL;
	xs[3]=0x3331C90D2C6BA52LL;          	xs[4]=0xC6203913F6ECC5LL;           	xs[5]=0x1B2063B22FCF270LL;
	xs[6]=0x2878A3BFD9F42FCLL;      	    xs[7]=0x6277E432C8A5ACLL;           	xs[8]=0x752CB45C48648BLL;

	ys[0]=0xcLL;        ys[1]=0x0LL;        ys[2]=0x0LL;
	ys[3]=0x0LL;        ys[4]=0x0LL;        ys[5]=0x0LL;
	ys[6]=0x0LL;    	ys[7]=0x0LL;    	ys[8]=0x0LL;



// The same 519-bit random scalar of Scott is used here
#ifndef TEST
#if WINDOW==4
w[0]= 2; w[1]= 0; w[2]= -4; w[3]= -7; w[4]= 1; w[5]= 7; w[6]= -3; w[7]= 5; w[8]= -4; w[9]= -1; w[10]= -7; w[11]= -2; w[12]= 4; w[13]= -5; w[14]= -3;
w[15]= 7; w[16]= -3; w[17]= -2; w[18]= 6; w[19]= 4; w[20]= 1; w[21]= -1; w[22]= -1; w[23]= 3; w[24]= 0; w[25]= -4; w[26]= 7; w[27]= 7; w[28]= -5; w[29]= 5;
w[30]= -5; w[31]= -3; w[32]= -6; w[33]= -7; w[34]= -7; w[35]= -7; w[36]= 4; w[37]= -8; w[38]= -7; w[39]= -5; w[40]= -7; w[41]= 2; w[42]= -5; w[43]= 6; w[44]= -3;
w[45]= 0; w[46]= -7; w[47]= 2; w[48]= 7; w[49]= 0; w[50]= -4; w[51]= -6; w[52]= 3; w[53]= 4; w[54]= -8; w[55]= 2; w[56]= 3; w[57]= -1; w[58]= 5; w[59]= -8;
w[60]= -2; w[61]= 6; w[62]= 0; w[63]= -6; w[64]= 5; w[65]= -6; w[66]= 1; w[67]= 4; w[68]= -7; w[69]= -1; w[70]= -1; w[71]= -7; w[72]= 6; w[73]= -5; w[74]= 5;
w[75]= 3; w[76]= -5; w[77]= 5; w[78]= 1; w[79]= 6; w[80]= -6; w[81]= -8; w[82]= -3; w[83]= 1; w[84]= -5; w[85]= -8; w[86]= 1; w[87]= -6; w[88]= -8; w[89]= -2;
w[90]= 4; w[91]= 3; w[92]= -6; w[93]= 1; w[94]= -2; w[95]= 0; w[96]= -2; w[97]= -3; w[98]= -3; w[99]= 5; w[100]= 0; w[101]= 0; w[102]= -1; w[103]= 4; w[104]= 2;
w[105]= 5; w[106]= 0; w[107]= 4; w[108]= 5; w[109]= -8; w[110]= -1; w[111]= -6; w[112]= -1; w[113]= -1; w[114]= -6; w[115]= -2; w[116]= 6; w[117]= -8; w[118]= -3; w[119]= 2;
w[120]= 1; w[121]= 2; w[122]= -3; w[123]= -7; w[124]= 6; w[125]= -8; w[126]= -2; w[127]= -2; w[128]= -1; w[129]= -8; w[130]= 1;
#endif
#if WINDOW==5
w[0]= 2; w[1]= 0; w[2]= 3; w[3]= 1; w[4]= -9; w[5]= 7; w[6]= -15; w[7]= -2; w[8]= -7; w[9]= -1; w[10]= 13; w[11]= -7; w[12]= -9; w[13]= 15; w[14]= -9;
w[15]= 9; w[16]= -15; w[17]= -8; w[18]= 12; w[19]= 0; w[20]= 12; w[21]= -5; w[22]= 14; w[23]= 9; w[24]= 11; w[25]= 14; w[26]= 2; w[27]= -15; w[28]= -7; w[29]= 2;
w[30]= 2; w[31]= -11; w[32]= -7; w[33]= -7; w[34]= -9; w[35]= -5; w[36]= -16; w[37]= 13; w[38]= -4; w[39]= 1; w[40]= -4; w[41]= -11; w[42]= -15; w[43]= -15; w[44]= -14;
w[45]= -6; w[46]= -12; w[47]= -15; w[48]= -2; w[49]= 3; w[50]= 8; w[51]= 9; w[52]= 10; w[53]= 0; w[54]= 5; w[55]= -3; w[56]= 15; w[57]= 12; w[58]= 13; w[59]= 9;
w[60]= -13; w[61]= 6; w[62]= 5; w[63]= 12; w[64]= -6; w[65]= 4; w[66]= 3; w[67]= -10; w[68]= 8; w[69]= -16; w[70]= -1; w[71]= -5; w[72]= -12; w[73]= -14; w[74]= 3;
w[75]= -4; w[76]= 0; w[77]= 7; w[78]= -13; w[79]= 10; w[80]= 0; w[81]= -8; w[82]= -16; w[83]= 5; w[84]= 5; w[85]= 0; w[86]= -11; w[87]= -15; w[88]= -1; w[89]= -11;
w[90]= -4; w[91]= -12; w[92]= -2; w[93]= 3; w[94]= -14; w[95]= 4; w[96]= 1; w[97]= 9; w[98]= 3; w[99]= 11; w[100]= -8; w[101]= 15; w[102]= -5; w[103]= -16; w[104]= 1;

#endif
#if WINDOW==6
w[0]= 2; w[1]= -16; w[2]= 9; w[3]= 28; w[4]= 13; w[5]= -15; w[6]= 15; w[7]= -10; w[8]= -12; w[9]= -13; w[10]= 23; w[11]= -9; w[12]= 6; w[13]= 5; w[14]= -17;
w[15]= 12; w[16]= 0; w[17]= 27; w[18]= -9; w[19]= 19; w[20]= 11; w[21]= -25; w[22]= 9; w[23]= -30; w[24]= 4; w[25]= -30; w[26]= 11; w[27]= 6; w[28]= 27; w[29]= -11;
w[30]= 16; w[31]= 6; w[32]= 7; w[33]= -16; w[34]= -22; w[35]= 17; w[36]= 24; w[37]= 12; w[38]= 15; w[39]= -31; w[40]= 30; w[41]= 1; w[42]= 10; w[43]= -23; w[44]= 1;
w[45]= -27; w[46]= -17; w[47]= -28; w[48]= -10; w[49]= 19; w[50]= -13; w[51]= 19; w[52]= -31; w[53]= -22; w[54]= 8; w[55]= 3; w[56]= -5; w[57]= 2; w[58]= -6; w[59]= -10;
w[60]= -12; w[61]= -23; w[62]= -31; w[63]= 0; w[64]= 14; w[65]= -13; w[66]= 5; w[67]= 0; w[68]= -1; w[69]= 9; w[70]= 5; w[71]= 16; w[72]= 5; w[73]= -6; w[74]= -22;
w[75]= -4; w[76]= 26; w[77]= 23; w[78]= 8; w[79]= 7; w[80]= -31; w[81]= -11; w[82]= 25; w[83]= -31; w[84]= 30; w[85]= -5; w[86]= 8;
#endif

// Group Order - for testing
#else

#if WINDOW==6
w[0]= -21; w[1]= -10; w[2]= 1; w[3]= 6; w[4]= -11; w[5]= 24; w[6]= 3; w[7]= 9; w[8]= -22; w[9]= 4; w[10]= 20; w[11]= 17; w[12]= 31; w[13]= -4; w[14]= -23; w[15]= -25; w[16]= 23; w[17]= 17;
w[18]= 12; w[19]= -10; w[20]= -4; w[21]= 20; w[22]= -16; w[23]= 16; w[24]= 5; w[25]= -5; w[26]= -24; w[27]= 24; w[28]= -17; w[29]= -29;
w[30]= -20; w[31]= 14; w[32]= -9; w[33]= 24; w[34]= 8; w[35]= -1; w[36]= 7; w[37]= 29; w[38]= -28; w[39]= -14; w[40]= -9; w[41]= 23;
w[42]= 17; w[43]= -1; w[44]= 0; w[45]= 0; w[46]= 0; w[47]= 0; w[48]= 0; w[49]= 0; w[50]= 0; w[51]= 0; w[52]= 0; w[53]= 0; w[54]= 0; w[55]= 0; w[56]= 0; w[57]= 0; w[58]= 0; w[59]= 0;
w[60]= 0; w[61]= 0; w[62]= 0; w[63]= 0; w[64]= 0; w[65]= 0; w[66]= 0; w[67]= 0; w[68]= 0; w[69]= 0; w[70]= 0; w[71]= 0; w[72]= 0; w[73]= 0; w[74]= 0; w[75]= 0; w[76]= 0; w[77]= 0;
w[78]= 0; w[79]= 0; w[80]= 0; w[81]= 0; w[82]= 0; w[83]= 0; w[84]= 0; w[85]= 0; w[86]= 8;
#endif

#if WINDOW==5

w[0]= 11; w[1]= 11; w[2]= 3; w[3]= -16; w[4]= -14; w[5]= -5; w[6]= -8; w[7]= 7; w[8]= 4; w[9]= -15; w[10]= -5; w[11]= 2; w[12]= -12; w[13]= 3; w[14]= -3;
w[15]= 4; w[16]= 15; w[17]= -12; w[18]= 7; w[19]= 13; w[20]= 5; w[21]= 2; w[22]= 3; w[23]= -5; w[24]= -4; w[25]= 8; w[26]= 1; w[27]= -2; w[28]= -12; w[29]= 3;
w[30]= -5; w[31]= -16; w[32]= -1; w[33]= -5; w[34]= 12; w[35]= -15; w[36]= 12; w[37]= -5; w[38]= -3; w[39]= -1; w[40]= 6; w[41]= 4; w[42]= -1; w[43]= 14; w[44]= -12;
w[45]= 4; w[46]= -7; w[47]= -7; w[48]= -9; w[49]= 14; w[50]= 5; w[51]= -6; w[52]= 0; w[53]= 0; w[54]= 0; w[55]= 0; w[56]= 0; w[57]= 0; w[58]= 0; w[59]= 0;
w[60]= 0; w[61]= 0; w[62]= 0; w[63]= 0; w[64]= 0; w[65]= 0; w[66]= 0; w[67]= 0; w[68]= 0; w[69]= 0; w[70]= 0; w[71]= 0; w[72]= 0; w[73]= 0; w[74]= 0;
w[75]= 0; w[76]= 0; w[77]= 0; w[78]= 0; w[79]= 0; w[80]= 0; w[81]= 0; w[82]= 0; w[83]= 0; w[84]= 0; w[85]= 0; w[86]= 0; w[87]= 0; w[88]= 0; w[89]= 0;
w[90]= 0; w[91]= 0; w[92]= 0; w[93]= 0; w[94]= 0; w[95]= 0; w[96]= 0; w[97]= 0; w[98]= 0; w[99]= 0; w[100]= 0; w[101]= 0; w[102]= 0; w[103]= -16;
w[104]= 1;

#endif

#if WINDOW==4
w[0]= -5; w[1]= 7; w[2]= -3; w[3]= 1; w[4]= -8; w[5]= 2; w[6]= 5; w[7]= -1; w[8]= 6; w[9]= 3; w[10]= 4; w[11]= 2; w[12]= -6; w[13]= -1; w[14]= 1;
w[15]= 4; w[16]= 5; w[17]= 4; w[18]= -1; w[19]= 2; w[20]= -1; w[21]= -7; w[22]= -5; w[23]= -6; w[24]= 7; w[25]= 5; w[26]= 4; w[27]= -4; w[28]= -7; w[29]= -2;
w[30]= -4; w[31]= 0; w[32]= 5; w[33]= 0; w[34]= -1; w[35]= 4; w[36]= 5; w[37]= -4; w[38]= -1; w[39]= -8; w[40]= -1; w[41]= 6; w[42]= -1; w[43]= -5; w[44]= -7;
w[45]= -4; w[46]= 7; w[47]= 3; w[48]= 7; w[49]= -1; w[50]= 6; w[51]= -8; w[52]= -3; w[53]= 0; w[54]= 7; w[55]= 4; w[56]= 7; w[57]= 4; w[58]= 6; w[59]= -4;
w[60]= 7; w[61]= -5; w[62]= 6; w[63]= 1; w[64]= -3; w[65]= 0; w[66]= 0; w[67]= 0; w[68]= 0; w[69]= 0; w[70]= 0; w[71]= 0; w[72]= 0; w[73]= 0; w[74]= 0;
w[75]= 0; w[76]= 0; w[77]= 0; w[78]= 0; w[79]= 0; w[80]= 0; w[81]= 0; w[82]= 0; w[83]= 0; w[84]= 0; w[85]= 0; w[86]= 0; w[87]= 0; w[88]= 0; w[89]= 0;
w[90]= 0; w[91]= 0; w[92]= 0; w[93]= 0; w[94]= 0; w[95]= 0; w[96]= 0; w[97]= 0; w[98]= 0; w[99]= 0; w[100]= 0; w[101]= 0; w[102]= 0; w[103]= 0; w[104]= 0;
w[105]= 0; w[106]= 0; w[107]= 0; w[108]= 0; w[109]= 0; w[110]= 0; w[111]= 0; w[112]= 0; w[113]= 0; w[114]= 0; w[115]= 0; w[116]= 0; w[117]= 0; w[118]= 0; w[119]= 0;
w[120]= 0; w[121]= 0; w[122]= 0; w[123]= 0; w[124]= 0; w[125]= 0; w[126]= 0; w[127]= 0; w[128]= 0; w[129]= -8; w[130]= 1;
#endif
#endif

	for(i=0; i<10; i++){
		bef=rdtsc();
		for (j=0; j<lpz; j++){
			init(xs,ys,&P);
			scalar_mul(w,&P);
		}
		aft=rdtscp();
		if((aft-bef)/(lpz) < min_ccycle)
				min_ccycle= (aft-bef)/(lpz);
	}
	printf("Window width :: %u \nThe minimum clock cycles count is :: %lu\n", WINDOW, min_ccycle);
	output(&P);

}
