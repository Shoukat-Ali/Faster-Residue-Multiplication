/* Test program for ws521 scalar point multiplication
 This is NIST standard Weierstrass curve p-521
 Fully Tested and debugged
 Uses constant time method described by Bos et al. - http://eprint.iacr.org/2014/130
 Cache safety thanks to ed25519
 We have replaced gmul() by TMV_product()
 We have added rdtscp() to measure the clock cycles count as suggested by Paoloni; http://www.intel.com.tr/content/dam/www/public/us/en/documents/white-papers/ia-32-ia-64-benchmark-code-execution-paper.pdf
 We have amended scr() and gsqr() according to modulus p
 gcc -Wall -O3 vb64_ws521.c -o vb64_ws521.exe */

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

////////////////////////////////////////////////////////////////////////////////////////
// w=1
void gone(type64 w[]){
	w[0]=1;
	w[1]=0;
	w[2]=0;
	w[3]=0;
	w[4]=0;
	w[5]=0;
	w[6]=0;
	w[7]=0;
	w[8]=0;
}

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

// fused operations
// w+=x
void ginc(const type64 x[], type64 w[]){
	w[0]+=x[0];
	w[1]+=x[1];
	w[2]+=x[2];
	w[3]+=x[3];
	w[4]+=x[4];
	w[5]+=x[5];
	w[6]+=x[6];
	w[7]+=x[7];
	w[8]+=x[8];
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

//w=w-2x-y
void gtsb2(const type64 x[], const type64 y[], type64 w[]){
	w[0]-=2*x[0]+y[0];
	w[1]-=2*x[1]+y[1];
	w[2]-=2*x[2]+y[2];
	w[3]-=2*x[3]+y[3];
	w[4]-=2*x[4]+y[4];
	w[5]-=2*x[5]+y[5];
	w[6]-=2*x[6]+y[6];
	w[7]-=2*x[7]+y[7];
	w[8]-=2*x[8]+y[8];
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

// w*=3
void gmul3(type64 w[]){
	w[0]*=3;
	w[1]*=3;
	w[2]*=3;
	w[3]*=3;
	w[4]*=3;
	w[5]*=3;
	w[6]*=3;
	w[7]*=3;
	w[8]*=3;
}

// w=2(x-y)
void g2sb(const type64 x[], const type64 y[], type64 w[]){
	w[0]=2*(x[0]-y[0]);
	w[1]=2*(x[1]-y[1]);
	w[2]=2*(x[2]-y[2]);
	w[3]=2*(x[3]-y[3]);
	w[4]=2*(x[4]-y[4]);
	w[5]=2*(x[5]-y[5]);
	w[6]=2*(x[6]-y[6]);
	w[7]=2*(x[7]-y[7]);
	w[8]=2*(x[8]-y[8]);
}

// w=3(x+y)
void g3ad(const type64 x[], const type64 y[], type64 w[]){
	w[0]=3*(x[0]+y[0]);
	w[1]=3*(x[1]+y[1]);
	w[2]=3*(x[2]+y[2]);
	w[3]=3*(x[3]+y[3]);
	w[4]=3*(x[4]+y[4]);
	w[5]=3*(x[5]+y[5]);
	w[6]=3*(x[6]+y[6]);
	w[7]=3*(x[7]+y[7]);
	w[8]=3*(x[8]+y[8]);
}


// w=4*w-x
void g4sb(const type64 x[], type64 w[]){
	w[0]=4*w[0]-x[0];
	w[1]=4*w[1]-x[1];
	w[2]=4*w[2]-x[2];
	w[3]=4*w[3]-x[3];
	w[4]=4*w[4]-x[4];
	w[5]=4*w[5]-x[5];
	w[6]=4*w[6]-x[6];
	w[7]=4*w[7]-x[7];
	w[8]=4*w[8]-x[8];
}

// w*=4
void gmul4(type64 w[]){
	w[0]*=4;
	w[1]*=4;
	w[2]*=4;
	w[3]*=4;
	w[4]*=4;
	w[5]*=4;
	w[6]*=4;
	w[7]*=4;
	w[8]*=4;
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

// w-=8*x
void gsb8(const type64 x[], type64 w[]){
	w[0]-=8*x[0];
	w[1]-=8*x[1];
	w[2]-=8*x[2];
	w[3]-=8*x[3];
	w[4]-=8*x[4];
	w[5]-=8*x[5];
	w[6]-=8*x[6];
	w[7]-=8*x[7];
	w[8]-=8*x[8];
}

// reduce w - Short Coefficient Reduction
// Changes are made according to the modulus p

void scr(type64 w[]){
	type64 t0,t1;
	t0=w[0]&bot58bits;

	t1=w[1]+(w[0]>>58);             w[1]=t1&bot58bits;
	t1=w[2]+(t1>>58);               w[2]=t1&bot58bits;
	t1=w[3]+(t1>>58);           	w[3]=t1&bot58bits;
	t1=w[4]+(t1>>58);           	w[4]=t1&bot58bits;
	t1=w[5]+(t1>>58);           	w[5]=t1&bot58bits;
	t1=w[6]+(t1>>58);           	w[6]=t1&bot58bits;
	t1=w[7]+(t1>>58);           	w[7]=t1&bot58bits;
	t1=w[8]+(t1>>58);           	w[8]=t1&bot57bits;				//w[8] is 57-bit

	w[0]=t0+(t1>>57);           // w[0] can be at most 59-bit
}

// z=x^2
// For modulus p the changes are commented

void gsqr(const type64 x[], type64 z[]){
	type128 t0,t1,t2;

	t1=2*((type128)x[0]*x[8]+(type128)x[1]*x[7]+(type128)x[2]*x[6]+(type128)x[3]*x[5])+(type128)x[4]*x[4];
	t0=((type64) t1)&bot57bits;			//z[8] is 57-bit
	t2=4*((type128)x[1]*x[8]+(type128)x[2]*x[7]+(type128)x[3]*x[6]+(type128)x[4]*x[5])+(type128)x[0]*x[0]+(t1>>57);		//modulus p
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
	z[8]=((type64)t0)&bot57bits;		//z[8] is 57-bit
	z[0]+=(type64)(t0>>57);				//modulus p
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


//Fermat's Little Theorem
// Inverse x = 1/x = x^(p-2) mod p
// 13 muls, 520 sqrs
//
void ginv(type64 x[]){
	type64 x127[9],w[9],t[9],z[9];
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
	TMV_product(z,w,t);		// z=z256

	gcopy(t,w);
	for (int i=0;i<128;i++)
	{
		gsqr(t,z);
		gsqr(z,t);
	}
	TMV_product(t,w,z);		// z=z512

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

////////////////////////////////////////////////////////////////////////////////////////////////////
// Point Structure
typedef struct{
    type64 x[LIMBS];
    type64 y[LIMBS];
    type64 z[LIMBS];
    int inf;
} ECp;


// Set P to point at infinity
void inf(ECp *P){
	P->inf=1;
}

// Initializing P to an Affine Point where P is in Jacobian coordinate
void init(const type64 x[], const type64 y[], ECp *P){
	for (int i=0; i<LIMBS; i++){
		P->x[i]=x[i];
		P->y[i]=y[i];
		P->z[i]=0;
	}
	P->z[0]=1;
	P->inf=0;
}

// P=Q
void copy(ECp *Q,ECp *P){
	for (int i=0; i<LIMBS; i++){
		P->x[i]=Q->x[i];
		P->y[i]=Q->y[i];
		P->z[i]=Q->z[i];
	}
	P->inf=Q->inf;
}

// P=-Q, Setting Negation of a Point
void neg(ECp *Q,ECp *P){
	for (int i=0; i<LIMBS; i++){
		P->x[i]=Q->x[i];
		P->y[i]=-Q->y[i];
		P->z[i]=Q->z[i];
	}
	P->inf=Q->inf;
}

// From Jacobian to Affine coordinate
void norm(ECp *P){
	type64 iz2[LIMBS], iz3[LIMBS], w[LIMBS], t[LIMBS];
	if (P->inf) return;
	gcopy(P->z,w);
	ginv(w);
	gsqr(w,iz2);
	TMV_product(iz2,w,iz3);
	TMV_product(P->x,iz2,t);        //t is in intermediate form
	scr(t);                         //t is a unique residue in modulo p
	gcopy(t,P->x);
	TMV_product(P->y,iz3,t);        //t is in intermediate form
	scr(t);                         //t is a unique residue in modulo p
	gcopy(t,P->y);
	gone(P->z);
}

// P=2P, Point doubling where P is Affine
void doubl_1(ECp *P){
    type64 r0[9],r1[9],r2[9],r3[9];

	gcopy(P->y, P->z);
	gmul2(P->z);
	scr(P->z);

	gsqr(P->x,r0);
	gsqr(P->y,r1);
	gsqr(r1,r2);
	ginc(P->x,r1);
	gsqr(r1,r3);
	gtsb(r0,r2,r3);
	gmul2(r3);
	r0[0]-=1;
	gmul3(r0);
	gsqr(r0,P->x);
	gsb2(r3,P->x);
	scr(P->x);
	gdec(P->x,r3);
	TMV_product(r3,r0,P->y);
	gsb8(r2,P->y);
	scr(P->y);

}

// P=2P, Point doubling where P is in Jacobian coordinate
void doubl_2(ECp *P){
	type64 r0[9],r1[9],r2[9],r3[9];

	gsqr(P->z,r0);
	gsqr(P->y,r1);
	TMV_product(P->x,r1,r2);
	gadd(P->y,P->z,r3);
	g3ad(P->x,r0,P->y);
	gdec(r0,P->x);
	TMV_product(P->x,P->y,P->z);
	gsqr(P->z,P->x);
	gsb8(r2,P->x);
	scr(P->x);
	g4sb(P->x,r2);
	TMV_product(P->z,r2,P->y);
	gsqr(r1,r2);
	gsb8(r2,P->y);
	scr(P->y);
	gsqr(r3,P->z);
	gtsb(r0,r1,P->z);
	scr(P->z);
}


// P+=Q, Mixed Point addition where Q is Jacobian and P is Affine but the resultant P will be Jacobian
void mixed_add_1(ECp *Q, ECp *P){
    type64 r0[9],r1[9],r2[9];

	gsqr(Q->z,r0);
	TMV_product(P->x,r0,r1);
	TMV_product(P->y,r0,r2);
	TMV_product(r2,Q->z,P->x);
	gdec(Q->x,r1);
	gsqr(r1,r2);
	gadd(Q->z,r1,P->y);
	gsqr(P->y, P->z);
	gtsb(r0,r2,P->z);
	scr(P->z);

	gmul4(r2);
	TMV_product(r2,r1,r0);
	g2sb(P->x,Q->y,r1);
	gsqr(r1,P->x);
	TMV_product(Q->x,r2,P->y);
	gtsb2(P->y,r0,P->x);
	scr(P->x);

	gcopy(P->y,r2);
	gdec(P->x,r2);
	TMV_product(r1,r2,P->y);
	TMV_product(r0,Q->y,r1);
	gsb2(r1,P->y);
	scr(P->y);

}

// P+=Q, Mixed Point addition where P is Jacobian and Q is Affine
void mixed_add_2(ECp *Q, ECp *P){
	type64 r0[9],r1[9],r2[9],r3[9],r4[9];

	gsqr(P->z,r3);  // not if its one
	TMV_product(Q->x,r3,r4); // ditto
	TMV_product(Q->y,r3,r2); //ditto
	TMV_product(r2,P->z,r1); // not if its one

	gdec(P->x,r4);
	g2sb(r1,P->y,r2);
	gsqr(r4,r0);
	gadd(P->z,r4,r1);
	gsqr(r1,P->z);
	gtsb(r3,r0,P->z);
	scr(P->z);

	gmul4(r0);
	TMV_product(r4,r0,r1);
	TMV_product(r1,P->y,r3);
 	TMV_product(r0,P->x,P->y);
	gsqr(r2,P->x);
	gtsb2(P->y,r1,P->x);
	scr(P->x);

	gsub(P->y,P->x,r0);
	TMV_product(r0,r2,P->y);
	gsb2(r3,P->y);
	scr(P->y);
}

// P+=Q, Point addition where P and Q are Jacobian
void Jac_add(ECp *Q,ECp *P){
	type64 z1z1[9],z2z2[9],u1[9],u2[9],s1[9],s2[9],h[9],i[9],j[9];

	gsqr(P->z,z1z1);
	gsqr(Q->z,z2z2);    // Q->z=1, z2z2=1
	TMV_product(P->x,z2z2,u1);
	TMV_product(Q->x,z1z1,u2);
	TMV_product(P->y,Q->z,h);
	TMV_product(h,z2z2,s1);
	TMV_product(Q->y,P->z,h);
	TMV_product(h,z1z1,s2);
	gsub(u2,u1,h);
	gcopy(h,j); gmul2(j); gsqr(j,i);

	TMV_product(h,i,j);
	TMV_product(u1,i,u2);
	g2sb(s2,s1,u1);
	gsqr(u1,i);
	gsub(i,j,P->x);
	gsb2(u2,P->x);
	scr(P->x);

	TMV_product(s1,j,i);
	gsub(u2,P->x,j);
	TMV_product(j,u1,P->y);
	gsb2(i,P->y);
	scr(P->y);

	gadd(P->z,Q->z,i);
	gsqr(i,j);
	gtsb(z1z1,z2z2,j);
	TMV_product(h,j,P->z);
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

/* Normalise all of P[i] using one inversion - Montgomery's trick
 Assume P[0] is already Affine */
void multi_norm(ECp P[]){
	int i;
	type64 t1[LIMBS], t2[LIMBS], t3[LIMBS], w[M][LIMBS];
	gone(w[1]);   // 0->1
	gcopy(P[1].z,w[2]); // 0-1, 1-2
	for (i=3;i<M;i++)   // 2-3
		TMV_product(w[i-1],P[i-1].z,w[i]);

	TMV_product(w[M-1],P[M-1].z,t1);
	ginv(t1);

	gcopy(P[M-1].z,t2);
	TMV_product(w[M-1],t1,t3);
	gcopy(t3,w[M-1]);

	for (i=M-2;;i--)
	{
		if (i==1)  // 0-1
		{
			TMV_product(t1,t2,w[1]); //0-1
			break;
		}
		TMV_product(w[i],t2,t3);
		TMV_product(t3,t1,w[i]);
		TMV_product(t2,P[i].z,t3);
		gcopy(t3,t2);
	}

    for (i=1;i<M;i++)
    {
		gone(P[i].z);
		gsqr(w[i],t1);
		TMV_product(P[i].x,t1,t2);
		gcopy(t2,P[i].x);
		TMV_product(t1,w[i],t2);
		TMV_product(P[i].y,t2,t1);
		gcopy(t1,P[i].y);
	}

}


// Pre-computation
void precomp(ECp *P, ECp W[]){
	ECp Q;
	copy(P,&Q);
	doubl_1(&Q);        //because P is affine
	copy(P,&W[0]);
	copy(&W[0],&W[1]);
	mixed_add_1(&Q,&W[1]);//Q is Jacobian and W[1] is Affine

	/*for (int i=2;i<M;i++){
		copy(&W[i-1],&W[i]);
		Jac_add(&Q,&W[i]);
	}
	//The loop is slower than the following manual unrolling */
	copy(&W[1], &W[2]);
	Jac_add(&Q, &W[2]);
	copy(&W[2], &W[3]);
	Jac_add(&Q, &W[3]);
	copy(&W[3], &W[4]);
	Jac_add(&Q, &W[4]);
	copy(&W[4], &W[5]);
	Jac_add(&Q, &W[5]);
	copy(&W[5], &W[6]);
	Jac_add(&Q, &W[6]);
	copy(&W[6], &W[7]);
	Jac_add(&Q, &W[7]);
#if WINDOW>4
	copy(&W[7], &W[8]);
	Jac_add(&Q, &W[8]);
	copy(&W[8], &W[9]);
	Jac_add(&Q, &W[9]);
	copy(&W[9], &W[10]);
	Jac_add(&Q, &W[10]);
	copy(&W[10], &W[11]);
	Jac_add(&Q, &W[11]);
	copy(&W[11], &W[12]);
	Jac_add(&Q, &W[12]);
	copy(&W[12], &W[13]);
	Jac_add(&Q, &W[13]);
	copy(&W[13], &W[14]);
	Jac_add(&Q, &W[14]);
	copy(&W[14], &W[15]);
	Jac_add(&Q, &W[15]);
#if WINDOW>5
	copy(&W[15], &W[16]);
	Jac_add(&Q, &W[16]);
	copy(&W[16], &W[17]);
	Jac_add(&Q, &W[17]);
	copy(&W[17], &W[18]);
	Jac_add(&Q, &W[18]);
	copy(&W[18], &W[19]);
	Jac_add(&Q, &W[19]);
	copy(&W[19], &W[20]);
	Jac_add(&Q, &W[20]);
	copy(&W[20], &W[21]);
	Jac_add(&Q, &W[21]);
	copy(&W[21], &W[22]);
	Jac_add(&Q, &W[22]);
	copy(&W[22], &W[23]);
	Jac_add(&Q, &W[23]);
	copy(&W[23], &W[24]);
	Jac_add(&Q, &W[24]);
	copy(&W[24], &W[25]);
	Jac_add(&Q, &W[25]);
	copy(&W[25], &W[26]);
	Jac_add(&Q, &W[26]);
	copy(&W[26], &W[27]);
	Jac_add(&Q, &W[27]);
	copy(&W[27], &W[28]);
	Jac_add(&Q, &W[28]);
	copy(&W[28], &W[29]);
	Jac_add(&Q, &W[29]);
	copy(&W[29], &W[30]);
	Jac_add(&Q, &W[30]);
	copy(&W[30], &W[31]);
	Jac_add(&Q, &W[31]);
#endif
#endif
}

// Windows of width 4-6
void window(ECp *Q,ECp *P){
	doubl_2(P);
	doubl_2(P);
	doubl_2(P);
	doubl_2(P);
#if WINDOW>4
	doubl_2(P);
#if WINDOW>5
	doubl_2(P);
#endif
#endif
	mixed_add_2(Q,P);     //Q is Affine and P is Jacobian

}


//Constant-time table look-up - borrowed from ed25519
void fe_cmov(type64 f[], type64 g[], int ib){
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

  babs=(babs-1)/2;

  cmov(T,&W[0],equal(babs,0));  // conditional move
  cmov(T,&W[1],equal(babs,1));
  cmov(T,&W[2],equal(babs,2));
  cmov(T,&W[3],equal(babs,3));
  cmov(T,&W[4],equal(babs,4));
  cmov(T,&W[5],equal(babs,5));
  cmov(T,&W[6],equal(babs,6));
  cmov(T,&W[7],equal(babs,7));
#if WINDOW>4
  cmov(T,&W[8],equal(babs,8));
  cmov(T,&W[9],equal(babs,9));
  cmov(T,&W[10],equal(babs,10));
  cmov(T,&W[11],equal(babs,11));
  cmov(T,&W[12],equal(babs,12));
  cmov(T,&W[13],equal(babs,13));
  cmov(T,&W[14],equal(babs,14));
  cmov(T,&W[15],equal(babs,15));
#if WINDOW>5
  cmov(T,&W[16],equal(babs,16));
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
#endif
#endif
  neg(T,&MT);  // minus t
  cmov(T,&MT,m&1);
}

// Scalar/Point Multiplication - exponent is 521
void scalar_mul(int *w,ECp *P){
	ECp W[M], Q;

	precomp(P, W);           //All the precomputed table is in Jacobian coordinate
	multi_norm(W);          //The precomputed table is converted to Affine coordinate where Z-coordinate is 1

	copy(&W[(w[PANES-1]-1)/2],P);
	for (int i=PANES-2;i>=0;i--){
		ct_select(&Q, W, w[i]);
		window(&Q, P);
	}
	norm(P);

}

///////////////////////////////////////////////////////////////////////////////////////////////////

int main(){
	uint64_t bef,aft, min_ccycle=100000000;
	int i, j, w[PANES], lpz=10000;
	ECp P;
	type64 xs[LIMBS], ys[LIMBS];

/* Base point on NIST Curve */

	xs[0]=0x17E7E31C2E5BD66LL;              	xs[1]=0x22CF0615A90A6FELL;                  	xs[2]=0x127A2FFA8DE334LL;
	xs[3]=0x1DFBF9D64A3F877LL;              	xs[4]=0x6B4D3DBAA14B5ELL;                   	xs[5]=0x14FED487E0A2BD8LL;
	xs[6]=0x15B4429C6481390LL;              	xs[7]=0x3A73678FB2D988ELL;                  	xs[8]=0xC6858E06B70404LL;

	ys[0]=0xBE94769FD16650LL;               	ys[1]=0x31C21A89CB09022LL;                  	ys[2]=0x39013FAD0761353LL;
	ys[3]=0x2657BD099031542LL;              	ys[4]=0x3273E662C97EE72LL;                  	ys[5]=0x1E6D11A05EBEF45LL;
	ys[6]=0x3D1BD998F544495LL;              	ys[7]=0x3001172297ED0B1LL;                  	ys[8]=0x11839296A789A3BLL;

// The same odd random scalar of Scott is used here
#ifndef TEST
#if WINDOW==6
w[0]=13; w[1]=29; w[2]=-25; w[3]=-39; w[4]=-55; w[5]=53; w[6]=-35; w[7]=63; w[8]=-53; w[9]=-9; w[10]=43; w[11]=-15; w[12]=61; w[13]=-63; w[14]=-33;
w[15]=-13; w[16]=33; w[17]=-47; w[18]=-33; w[19]=-7; w[20]=-25; w[21]=21; w[22]=-53; w[23]=-35; w[24]=-39; w[25]=-25; w[26]=-23; w[27]=-63; w[28]=-59;
w[29]=-39; w[30]=45; w[31]=-5; w[32]=13; w[33]=-11; w[34]=7; w[35]=63; w[36]=27; w[37]=-5; w[38]=-41; w[39]=61; w[40]=-31; w[41]=-17; w[42]=23; w[43]=-39;
w[44]=15; w[45]=27; w[46]=-27; w[47]=55; w[48]=41; w[49]=-13; w[50]=59; w[51]=-41; w[52]=31; w[53]=41; w[54]=7; w[55]=3; w[56]=59; w[57]=-63; w[58]=59; w[59]=53;
w[60]=-13; w[61]=-23; w[62]=33; w[63]=63; w[64]=13; w[65]=-13; w[66]=-59; w[67]=1; w[68]=-1; w[69]=9; w[70]=-59; w[71]=17; w[72]=-59; w[73]=59; w[74]=41; w[75]=59;
w[76]=25; w[77]=-41; w[78]=9; w[79]=7; w[80]=-31; w[81]=-11; w[82]=25; w[83]=33; w[84]=29; w[85]=59; w[86]=15;
#endif
#if WINDOW==5
w[0]=-19; w[1]=27; w[2]=-3; w[3]=-27; w[4]=-25; w[5]=-27; w[6]=21; w[7]=27; w[8]=25; w[9]=-1; w[10]=3; w[11]=-5; w[12]=11; w[13]=3; w[14]=19;
w[15]=-17; w[16]=1; w[17]=-17; w[18]=19; w[19]=-31; w[20]=-25; w[21]=19; w[22]=-25; w[23]=-3; w[24]=7; w[25]=9; w[26]=13; w[27]=1; w[28]=-25;
w[29]=-19; w[30]=7; w[31]=-15; w[32]=-29; w[33]=1; w[34]=-31; w[35]=-19; w[36]=13; w[37]=23; w[38]=19; w[39]=9; w[40]=13; w[41]=3; w[42]=31;
w[43]=23; w[44]=13; w[45]=23; w[46]=-27; w[47]=31; w[48]=1; w[49]=-3; w[50]=-5; w[51]=-21; w[52]=7; w[53]=7; w[54]=-5; w[55]=-21; w[56]=-5;
w[57]=-17; w[58]=27; w[59]=-7; w[60]=27; w[61]=15; w[62]=25; w[63]=-21; w[64]=27; w[65]=3; w[66]=-29; w[67]=23; w[68]=-25; w[69]=-15; w[70]=-1;
w[71]=27; w[72]=19; w[73]=-15; w[74]=-29; w[75]=29; w[76]=-1; w[77]=7; w[78]=19; w[79]=-23; w[80]=-31; w[81]=25; w[82]=-17; w[83]=5; w[84]=-27;
w[85]=1; w[86]=-11; w[87]=-15; w[88]=-1; w[89]=21; w[90]=27; w[91]=19; w[92]=-3; w[93]=-29; w[94]=19; w[95]=3; w[96]=1; w[97]=9; w[98]=3; w[99]=-21;
w[100]=-7; w[101]=15; w[102]=27; w[103]=-1; w[104]=1;
#endif
#if WINDOW==4
w[0]=-3; w[1]=5; w[2]=7; w[3]=-9; w[4]=-13; w[5]=-9; w[6]=-7; w[7]=1; w[8]=13; w[9]=13; w[10]=9; w[11]=15; w[12]=-5; w[13]=9; w[14]=-3; w[15]=-5;
w[16]=-9; w[17]=-3; w[18]=13; w[19]=-9; w[20]=-15; w[21]=15; w[22]=-7; w[23]=-3; w[24]=-15; w[25]=-9; w[26]=-11; w[27]=15; w[28]=-15; w[29]=-1;
w[30]=-9; w[31]=3; w[32]=5; w[33]=-5; w[34]=1; w[35]=-9; w[36]=9; w[37]=9; w[38]=-7; w[39]=-7; w[40]=-13; w[41]=-15; w[42]=-11; w[43]=-15; w[44]=-9;
w[45]=-3; w[46]=-1; w[47]=-1; w[48]=-3; w[49]=5; w[50]=-3; w[51]=-9; w[52]=13; w[53]=15; w[54]=11; w[55]=-3; w[56]=-1; w[57]=7; w[58]=1; w[59]=15; w[60]=-15;
w[61]=11; w[62]=-5; w[63]=7; w[64]=-11; w[65]=-9; w[66]=-1; w[67]=-3; w[68]=7; w[69]=-11; w[70]=11; w[71]=13; w[72]=-7; w[73]=-1; w[74]=-3; w[75]=11; w[76]=15;
w[77]=-11; w[78]=15; w[79]=-11; w[80]=11; w[81]=-9; w[82]=-3; w[83]=1; w[84]=11; w[85]=-9; w[86]=-15; w[87]=11; w[88]=7; w[89]=13; w[90]=3; w[91]=-13; w[92]=-5;
w[93]=-15; w[94]=15; w[95]=15; w[96]=-3; w[97]=-3; w[98]=-3; w[99]=-11; w[100]=-15; w[101]=1; w[102]=15; w[103]=-13; w[104]=3; w[105]=-11; w[106]=-15; w[107]=5;
w[108]=-11; w[109]=-7; w[110]=15; w[111]=-7; w[112]=-1; w[113]=15; w[114]=9; w[115]=13; w[116]=-11; w[117]=-7; w[118]=13; w[119]=1; w[120]=-15; w[121]=3; w[122]=-3;
w[123]=9; w[124]=-11; w[125]=9; w[126]=13; w[127]=-3; w[128]=15; w[129]=-1; w[130]=1;
#endif
#else
#if WINDOW==6
w[0]=-55; w[1]=-47; w[2]=-57; w[3]=15; w[4]=-47; w[5]=59; w[6]=49; w[7]=45; w[8]=47; w[9]=45; w[10]=43; w[11]=43; w[12]=7; w[13]=49; w[14]=-39;
w[15]=-29; w[16]=-7; w[17]=-25; w[18]=29; w[19]=45; w[20]=-5; w[21]=1; w[22]=29; w[23]=41; w[24]=-55; w[25]=29; w[26]=-49; w[27]=19; w[28]=-63;
w[29]=-15; w[30]=61; w[31]=31; w[32]=43; w[33]=25; w[34]=57; w[35]=11; w[36]=-1; w[37]=-49; w[38]=57; w[39]=-31; w[40]=-57; w[41]=7; w[42]=-27;
w[43]=63; w[44]=63; w[45]=63; w[46]=63; w[47]=63; w[48]=63; w[49]=63; w[50]=63; w[51]=63; w[52]=63; w[53]=63; w[54]=63; w[55]=63; w[56]=63; w[57]=63;
w[58]=63; w[59]=63; w[60]=63; w[61]=63; w[62]=63; w[63]=63; w[64]=63; w[65]=63; w[66]=63; w[67]=63; w[68]=63; w[69]=63; w[70]=63; w[71]=63; w[72]=63;
w[73]=63; w[74]=63; w[75]=63; w[76]=63; w[77]=63; w[78]=63; w[79]=63; w[80]=63; w[81]=63; w[82]=63; w[83]=63; w[84]=63; w[85]=63; w[86]=31;
#endif
#if WINDOW==5
w[0]=-23; w[1]=1; w[2]=-7; w[3]=17; w[4]=-13; w[5]=-23; w[6]=27; w[7]=3; w[8]=23; w[9]=29; w[10]=-5; w[11]=23; w[12]=11; w[13]=-9; w[14]=-1; w[15]=-23;
w[16]=-3; w[17]=-19; w[18]=3; w[19]=17; w[20]=-5; w[21]=5; w[22]=-9; w[23]=23; w[24]=27; w[25]=-31; w[26]=21; w[27]=-21; w[28]=-5; w[29]=-27; w[30]=-3;
w[31]=-1; w[32]=-23; w[33]=-21; w[34]=-31; w[35]=-7; w[36]=29; w[37]=31; w[38]=13; w[39]=-19; w[40]=-9; w[41]=29; w[42]=-21; w[43]=31; w[44]=27; w[45]=-31;
w[46]=-1; w[47]=-15; w[48]=-25; w[49]=-19; w[50]=-11; w[51]=21; w[52]=31; w[53]=31; w[54]=31; w[55]=31; w[56]=31; w[57]=31; w[58]=31; w[59]=31; w[60]=31;
w[61]=31; w[62]=31; w[63]=31; w[64]=31; w[65]=31; w[66]=31; w[67]=31; w[68]=31; w[69]=31; w[70]=31; w[71]=31; w[72]=31; w[73]=31; w[74]=31; w[75]=31; w[76]=31;
w[77]=31; w[78]=31; w[79]=31; w[80]=31; w[81]=31; w[82]=31; w[83]=31; w[84]=31; w[85]=31; w[86]=31; w[87]=31; w[88]=31; w[89]=31; w[90]=31; w[91]=31; w[92]=31;
w[93]=31; w[94]=31; w[95]=31; w[96]=31; w[97]=31; w[98]=31; w[99]=31; w[100]=31; w[101]=31; w[102]=31; w[103]=31; w[104]=1;
#endif
#if WINDOW==4
w[0]=-7; w[1]=-15; w[2]=-11; w[3]=-9; w[4]=9; w[5]=3; w[6]=1; w[7]=-7; w[8]=15; w[9]=1; w[10]=7; w[11]=11; w[12]=-1; w[13]=7; w[14]=11; w[15]=-5;
w[16]=-1; w[17]=11; w[18]=-9; w[19]=-11; w[20]=13; w[21]=9; w[22]=-7; w[23]=-7; w[24]=9; w[25]=11; w[26]=-7; w[27]=13; w[28]=5; w[29]=11; w[30]=11;
w[31]=-13; w[32]=1; w[33]=13; w[34]=-11; w[35]=11; w[36]=-7; w[37]=1; w[38]=7; w[39]=-1; w[40]=-7; w[41]=5; w[42]=-15; w[43]=-15; w[44]=-3; w[45]=13;
w[46]=15; w[47]=7; w[48]=-5; w[49]=-9; w[50]=7; w[51]=9; w[52]=-1; w[53]=3; w[54]=15; w[55]=11; w[56]=-13; w[57]=9; w[58]=-9; w[59]=-7; w[60]=-9;
w[61]=9; w[62]=1; w[63]=-11; w[64]=11; w[65]=15; w[66]=15; w[67]=15; w[68]=15; w[69]=15; w[70]=15; w[71]=15; w[72]=15; w[73]=15; w[74]=15; w[75]=15;
w[76]=15; w[77]=15; w[78]=15; w[79]=15; w[80]=15; w[81]=15; w[82]=15; w[83]=15; w[84]=15; w[85]=15; w[86]=15; w[87]=15; w[88]=15; w[89]=15; w[90]=15;
w[91]=15; w[92]=15; w[93]=15; w[94]=15; w[95]=15; w[96]=15; w[97]=15; w[98]=15; w[99]=15; w[100]=15; w[101]=15; w[102]=15; w[103]=15; w[104]=15; w[105]=15;
w[106]=15; w[107]=15; w[108]=15; w[109]=15; w[110]=15; w[111]=15; w[112]=15; w[113]=15; w[114]=15; w[115]=15; w[116]=15; w[117]=15; w[118]=15; w[119]=15;
w[120]=15; w[121]=15; w[122]=15; w[123]=15; w[124]=15; w[125]=15; w[126]=15; w[127]=15; w[128]=15; w[129]=15; w[130]=1;
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
	printf("Window width :: %u \nThe minimum mean clock cycles count is :: %lu\n", WINDOW, min_ccycle);
	output(&P);

}
