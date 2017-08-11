/*The prime order of the group, ECsGord, is taken from FIPS PUB 186-4
ECsGord=6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449
bitlength(ECsGord)= t=521
We have implemented the Algorithm 7, Fixed-base scalar multiplication, proposed by J. W. Bos et al.
The Point operation formulas are selected from Explicit Formulas Database https://www.hyperelliptic.org/EFD/index.html
We are counting the clock cycles of the Evaluation Stage of Algorithm 7 by using the fastest explicit formulas for (Jacobian) doubling and (mixed) addition.
The program is tested 10 times for 10000 iteration in order to obtain consistent cycle count
gcc -Wall -O3 fb64_ws521.c -o fb64_ws521.exe */


#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>


#define WINDOW 4 //6 //5

#define BIT_LENGTH 521

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

///////////////////////////////////////////////////////////////////////////////////////////

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

// w=2*w+x
void g2ad(const type64 x[], type64 w[]){
	w[0]=2*w[0]+x[0];
	w[1]=2*w[1]+x[1];
	w[2]=2*w[2]+x[2];
	w[3]=2*w[3]+x[3];
	w[4]=2*w[4]+x[4];
	w[5]=2*w[5]+x[5];
	w[6]=2*w[6]+x[6];
	w[7]=2*w[7]+x[7];
	w[8]=2*w[8]+x[8];
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


////////////////////////////////////////////////////////////////////////////////////////////
// Jacobian Point Structure
typedef struct{
    type64 x[LIMBS];
    type64 y[LIMBS];
    type64 z[LIMBS];
    type64 inf;
} ECp_Jac;

// Affine Point Structure
typedef struct{
    type64 x[LIMBS];
    type64 y[LIMBS];
} ECp_Aff;


// P=0, where P is Jacobian
void inf(ECp_Jac *P){
	P->inf=1;
}

// Initialise Jacobian point P as neutral element
void init(ECp_Jac *P){
	for(int i=0;i<LIMBS;i++){
		P->x[i]=0;
		P->y[i]=0;
		P->z[i]=0;
	}
	P->x[0]=1;
	P->y[0]=1;
	P->inf=0;
}


// P=Q where P is Jacobian and Q is Affine
void Aff_copy_Jac(const ECp_Aff *Q, ECp_Jac *P){
	for (int i=0;i<LIMBS;i++){
		P->x[i]=Q->x[i];
		P->y[i]=Q->y[i];
		P->z[i]=0;
	}
	P->z[0]=1;
	P->inf=0;
}

// P=Q where both P and Q are Jacobian
void Jac_copy(const ECp_Jac *Q, ECp_Jac *P){
	for (int i=0;i<LIMBS;i++){
		P->x[i]=Q->x[i];
		P->y[i]=Q->y[i];
		P->z[i]=Q->z[i];
	}
	P->inf=Q->inf;
}

//Affine P=-Q
void neg_Aff(ECp_Aff *Q, ECp_Aff *P){
    for (int32_t i=0;i<LIMBS;i++){
		P->x[i]=Q->x[i];
		P->y[i]=-Q->y[i];
	}
}


// From Jacobian to Affine coordinate
void norm(ECp_Jac *P){
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

// Printing the Affine coordinates of the computed new point on the screen
void output(ECp_Jac *P){
    int32_t i;
        for(i=0; i<LIMBS; i++)
        printf("x[%d]= %"PRIX64"\n", i, P->x[i]);

	puts("");
	for(i=0; i<LIMBS; i++)
        printf("y[%d]= %"PRIX64"\n", i, P->y[i]);

    puts("");

}

// P=2P, Point doubling where P is in Jacobian coordinate
void doubl_Jac(ECp_Jac *P){
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

// P+=Q, where P and Q are Jacobian
void add_Jac(ECp_Jac *Q, ECp_Jac *P){
	type64 z1z1[LIMBS], z2z2[LIMBS], u1[LIMBS], u2[LIMBS], s1[LIMBS], s2[LIMBS], h[LIMBS], i[LIMBS], j[LIMBS];

	gsqr(P->z,z1z1);
	gsqr(Q->z,z2z2);
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

// P+=Q, Mixed Point addition where P is Jacobian and Q is Affine
void add_mixed(ECp_Aff *Q, ECp_Jac *P){
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

// Normalise all of P[i] using one inversion - Montgomery's trick
// Assume P[0] is already Affine and table parameter V=3

void multi_norm(ECp_Jac P[][V], ECp_Aff Q[][V]){
	int i, j, sz=M*V;
	type64 t1[LIMBS], t2[LIMBS], w[sz][LIMBS];
	gone(w[0]);             //z=1
	gcopy(P[0][1].z, w[1]);
	TMV_product(w[1], P[0][2].z, w[2]);

	for(i=1; i<M; i++)
        for(j=0; j<V; j++)
            TMV_product(w[V*i+j-1], P[i][j].z, w[V*i+j]);

	gcopy(w[i*j-1], t1);
	ginv(t1);       //a^{-1}_{n}

	sz-=1;
	for(i=M-1; i>=0; i--){
        TMV_product(w[sz-1], t1, w[sz]);        //z^{-1}_{n}
        TMV_product(t1, P[i][2].z, t2);        //a^{-1}_{n-1}
        TMV_product(w[sz-2], t2, w[sz-1]);      //z^{-1}_{n-1}
        TMV_product(t2, P[i][1].z, t1);        //a^{-1}_{n-2}
        if(i!=0){
            TMV_product(w[sz-3], t1, w[sz-2]);      //z^{-1}_{n-2}
            TMV_product(t1, P[i][0].z, t2);        //a^{-1}_{n-3}
            gcopy(t2, t1);
            sz-=3;
        }
	}

	gcopy(P[0][0].x, Q[0][0].x);
    gcopy(P[0][0].y, Q[0][0].y);
    gsqr(w[1],t1);
    TMV_product(P[0][1].x, t1, Q[0][1].x);
	TMV_product(t1, w[1], t2);
	TMV_product(P[0][1].y, t2, Q[0][1].y);
	gsqr(w[2],t1);
    TMV_product(P[0][2].x, t1, Q[0][2].x);
	TMV_product(t1, w[2], t2);
	TMV_product(P[0][2].y, t2, Q[0][2].y);

    for(i=1; i<M; i++){
        for(j=0; j<V; j++){
            gsqr(w[3*i+j],t1);
            TMV_product(P[i][j].x, t1, Q[i][j].x);
            TMV_product(t1, w[3*i+j], t2);
            TMV_product(P[i][j].y, t2, Q[i][j].y);
        }
    }

}

// Off-line Pre-computation
void Offline_precomp(ECp_Jac *P, ECp_Aff W[][V]){
    ECp_Jac Q, Jac_PreT[M][V];
	ECp_Jac exp_d[WINDOW-1], exp_e[V];
	Jac_copy(P, &Q);
	Jac_copy(P, &exp_e[0]);
	int32_t i, j, ti;

	for(i=1; i<V; i++){
        for(j=0; j<E; j++)
            doubl_Jac(&Q);
        Jac_copy(&Q, &exp_e[i]);
    }
    Jac_copy(P, &Q);
    for(i=1; i<WINDOW; i++){
        for(j=0; j<D; j++)
            doubl_Jac(&Q);
        Jac_copy(&Q, &exp_d[i-1]);
    }

    for(i=0; i<M; i++){
        for(j=0; j<V; j++){
            ti=i;                           //computation to be preformed on the value of i
            Jac_copy(&exp_e[j], &Jac_PreT[i][j]);
            for(uint32_t jj=0; jj<WINDOW-1; jj++){      //For the bits of i
                if(ti&1){                       //checking Least Significant Bit (LSB)
                    Jac_copy(&exp_d[jj], &Q);
                    if(j){                      //for j>0, extra doubling required
                        for(uint32_t k=0; k<j*E; k++)
                            doubl_Jac(&Q);
                    }
                    add_Jac(&Q, &Jac_PreT[i][j]);
                }
                ti>>=1;
            }
        }
    }
    //printf("PreT[0][0].z = %"PRIu64"/n", Jac_PreT[0][0].z[0]);
    multi_norm(Jac_PreT, W);

}


//Constant-time table look-up - borrowed from ed25519

void fe_cmov(type64 f[], type64 g[], int ib){
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

static void cmov(ECp_Aff *w, ECp_Aff *u, int b){
  fe_cmov(w->x,u->x,b);
  fe_cmov(w->y,u->y,b);
}

// return 1 if b==c, no branching
static int equal(int b, int c){
	int x=b^c;
	x-=1;  // if x=0, x now -1
	return (x>>31);
}

static void ct_select(ECp_Aff *T, ECp_Aff W[][V], const uint32_t col, const int32_t s, const uint32_t b){
  ECp_Aff MT;

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
  neg_Aff(T,&MT);       // minus t
  cmov(T,&MT,s>>1);     //because s={-1, 1}, therefore, we shift only by two bits
}

// Scalar Multiplication - scalar is 521-bit
void scalar_mul(int32_t S[][E], uint32_t T[][E], ECp_Aff W[][V], ECp_Jac *JP){
    ECp_Aff AQ;

    ct_select(&AQ, W, 0, S[0][E-1], T[0][E-1]);
    Aff_copy_Jac(&AQ, JP);

	for(uint32_t i=1; i<V; i++){
        ct_select(&AQ, W, i, S[i][E-1], T[i][E-1]);
        add_mixed(&AQ, JP);
	}

    for(int32_t i=E-2; i>=0; i--){
        doubl_Jac(JP);
        for(uint32_t j=0; j<V; j++){
            ct_select(&AQ, W, j, S[j][i], T[j][i]);
            add_mixed(&AQ, JP);
        }
    }
    norm(JP);

}

/////////////////////////////////////////////////////////////////////////////////////////
int main(void){
	uint64_t bef, aft, min_ccycle=1000000000000000;
	int32_t i, j, lpz=10000;
	ECp_Aff PreT[M][V];
	ECp_Jac JP;

/* Base point on NIST Curve */

    JP.x[0]=0x17E7E31C2E5BD66LL;              	JP.x[1]=0x22CF0615A90A6FELL;                  	JP.x[2]=0x127A2FFA8DE334LL;
	JP.x[3]=0x1DFBF9D64A3F877LL;              	JP.x[4]=0x6B4D3DBAA14B5ELL;                   	JP.x[5]=0x14FED487E0A2BD8LL;
	JP.x[6]=0x15B4429C6481390LL;              	JP.x[7]=0x3A73678FB2D988ELL;                  	JP.x[8]=0xC6858E06B70404LL;

	JP.y[0]=0xBE94769FD16650LL;               	JP.y[1]=0x31C21A89CB09022LL;                  	JP.y[2]=0x39013FAD0761353LL;
	JP.y[3]=0x2657BD099031542LL;              	JP.y[4]=0x3273E662C97EE72LL;                  	JP.y[5]=0x1E6D11A05EBEF45LL;
	JP.y[6]=0x3D1BD998F544495LL;              	JP.y[7]=0x3001172297ED0B1LL;                  	JP.y[8]=0x11839296A789A3BLL;

	JP.z[0]=0x1LL;                          	JP.z[1]=0x0LL;                                  JP.z[2]=0x0LL;
	JP.z[3]=0x0LL;                          	JP.z[4]=0x0LL;                                  JP.z[5]=0x0LL;
	JP.z[6]=0x0LL;                          	JP.z[7]=0x0LL;                                  JP.z[8]=0x0LL;

    /*We have chosen the following 521-bit odd integer at random for testing
    Private key =Scalar=4513347784704396292717128854612512358876187488483123398413374978717743734485614150005012626109238630014295998715254757289500358152056299961832918676652788179
    The mLSB-set of the above private key for window width-4 is:
    mLSB-set={1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, 1, 1, -1,
    -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1,
    1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, 1, 1, 0, -1, -1, 0, 0, 0, 1, 0, -1, 0, 1, -1, -1, 0, 0, 1, 0,
    -1, 0, 1, 1, 1, 1, 1, 0, 0, -1, 0, 1, 1, 0, 0, -1, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, 0, -1, 1, 0, 0, 0, 0, 1, 0, 0, 0, -1, 1, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 1,
    0, 0, -1, 0, -1, 0, 0, 1, 0, 1, 1, -1, -1, 0, 1, 1, 0, -1, 1, -1, 0, -1, 1, 0, 0, 1, -1, 0, 1, 1, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1, 0, -1, 0,
    1, 0, 1, 0, 0, 0, -1, 0, 0, 0, 1, 1, 0, -1, 1, 0, -1, -1, 1, 1, -1, 0, -1, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 1, -1, -1, -1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, -1, 0,
    -1, 0, -1, 0, 1, 0, 1, 1, -1, -1, 1, -1, -1, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1, -1, -1, -1, 1, -1, 0, 1, -1, 0, -1, 0, 0, -1, -1, 1, 0, 0, -1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0,
    0, 1, 1, 0, 0, 0, 1, -1, -1, 0, -1, 1, 0, -1, 1, 1, 1, 0, 1, 0, 1, 0, -1, 1, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 1, -1, -1, 1, -1, 0, 0, 0, 1, -1, -1, -1, 1, 1, 1, 0, 1, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 1, -1, 1, 0, 0, 0, -1, -1, 0, 1, 0, -1, 1, 0, -1, 0, -1, 0, 0, -1, 0, 0, 1, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 1,
    0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, -1, 0, 0, 0, -1, -1, 0, 1, 0, 0, 0, 1, 0, -1, 1, -1, 1, 0, 0, 0, 0, 0, 0}
    The mLSB-set of the above private key for window width-5 is:
    mLSB-set={1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, 1, 1, -1, -1,
    -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1,
    -1, 1, 1, 1, -1, 1, 1, -1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, -1, -1, 1, 1, -1, 0, -1, 0, 1, 1, 0, 0, 1, -1, 0, 1, 1, 0, -1, 0, 0, 1, 0, 1, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 1, -1, 0, 0,
    0, -1, -1, 0, 0, 1, -1, 1, 0, -1, 0, 0, -1, 0, -1, 0, 0, 0, -1, -1, 1, 1, 0, 0, 0, 0, 0, -1, -1, -1, 1, -1, 1, 0, -1, -1, -1, 1, 0, 0, 0, 1, -1, 0, -1, 0, 0, 0, 0, -1, 1, 0, 1, -1,
    1, 1, -1, -1, 1, -1, 0, 0, 1, -1, -1, 0, -1, 0, -1, 0, 0, -1, -1, -1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, 1, -1, -1, 0, 1, 0, -1, -1,
    -1, 1, -1, 0, 0, 0, -1, 0, -1, 0, -1, -1, 0, -1, -1, 0, 0, 0, 1, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 1, -1, -1, -1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0, 1, -1, 0, 0, 1, 0, 1, 0, -1, -1, 0,
    -1, 1, 1, 1, -1, -1, 1, -1, -1, 0, 1, 1, -1, 0, 0, 0, 1, 1, 0, 0, 0, -1, -1, 1, 0, 0, -1, 0, -1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, -1, 1, 1, -1, -1, -1, 0, -1, -1, -1, -1, 0, 0, 1, 0,
    0, -1, 0, -1, 0, -1, -1, 1, -1, -1, 0, 0, 1, 1, 0, 1, 0, -1, -1, -1, -1, 1, -1, 0, 1, 0, 0, 0, 0, 1, -1, -1, 1, -1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1,
    0, -1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, 0, 1, 1, 0, -1, 0, 1, 1, -1, -1, 0, 0, 0, 1, 1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, -1, 0, 1, 0, 0, 1, -1, -1, -1, -1,
    1, 0, 0, -1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, -1, -1, 0, -1, -1, -1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0}
    The mLSB-set of the above private key for window width-5 is:
    mLSB-set={1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, 1, 1, -1, -1, -1,
    1, -1, -1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 0, -1, 0, 1, 0, 1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0,
    0, -1, 0, 0, 0, 0, 0, 1, 0, -1, -1, 1, 0, 0, -1, -1, -1, 0, 0, 0, 1, 1, 0, -1, -1, 1, 1, 1, -1, 1, 1, 0, -1, 0, 0, -1, -1, 0, 0, 1, -1, 1, 1, -1, -1, 1, -1, 0, -1, 0, 0, 0, 0, -1, 1, 1,
    1, 0, 1, 0, 0, -1, 0, 0, 1, -1, 1, 1, -1, 0, 1, 0, 0, -1, 1, 0, 0, 0, 0, -1, -1, 1, 0, -1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, -1, -1, -1, 0, 0, 1, 1, 1, 0, 0, -1, 0, 0,
    1, -1, 1, 1, 0, -1, -1, 0, 0, -1, -1, -1, 1, -1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, -1, 0, 1, 0, -1, -1, 1, 0, 1, 1, 1, -1, -1, 0, -1,
    0, -1, 1, 1, 0, 0, -1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, -1, 0, -1, 0, -1, 1, 0, 0, -1, 0, -1, 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, -1, 0, -1, 0, -1, 0, 1, 0, -1, 0, 0, -1, -1, 0, 1, -1, -1,
    0, 0, 0, 0, 0, 1, -1, -1, -1, -1, 0, 1, 0, 0, 0, -1, 0, 1, 1, -1, -1, 0, 0, 1, 0, 0, 0, 0, 1, -1, -1, 0, 0, 1, 0, -1, -1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, -1, -1, -1, 1, -1, 1, 1, 1, 0, -1,
    -1, 0, 1, 1, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 1, -1, 1, 0, -1, 0, 0, 0, -1, 0, -1, 1, 0, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 1, 0, 0, 0, 1, 0, -1, 1, 0, 1, 0, 0, -1, 0, 1, -1,
    0, -1, 0, 0, -1, 0, 0, 0, 1, 1, 1, 0, 1, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 1, 0, 0, 1, -1, -1, -1, 1, 0, 0, 0, -1, 0, 0, 1, 1, -1, -1, 1, 0, 0, 0, -1, 0, 0, 0, -1, 0,
    0, 0, 0, 1, 0, -1, -1, -1, -1, 1, -1, 1, 1, 0, -1, 1}
    */

    uint32_t TT[V][E];         //TT is the table, T, of the absolute values
    int32_t  SS[V][E];         //SS is the table, s, of (first row) signs

#ifndef TEST
#if WINDOW==6
    TT[0][0]=24; TT[0][1]=13; TT[0][2]=30; TT[0][3]=23; TT[0][4]=0; TT[0][5]=29; TT[0][6]=4; TT[0][7]=4; TT[0][8]=22; TT[0][9]=7; TT[0][10]=26; TT[0][11]=28; TT[0][12]=10; TT[0][13]=23; TT[0][14]=6; TT[0][15]=14; TT[0][16]=16; TT[0][17]=9; TT[0][18]=12; TT[0][19]=8; TT[0][20]=16; TT[0][21]=20; TT[0][22]=30; TT[0][23]=11; TT[0][24]=18; TT[0][25]=17; TT[0][26]=19; TT[0][27]=29; TT[0][28]=8;
    TT[1][0]=4; TT[1][1]=15; TT[1][2]=11; TT[1][3]=15; TT[1][4]=8; TT[1][5]=12; TT[1][6]=14; TT[1][7]=11; TT[1][8]=11; TT[1][9]=4; TT[1][10]=25; TT[1][11]=31; TT[1][12]=5; TT[1][13]=9; TT[1][14]=27; TT[1][15]=3; TT[1][16]=7; TT[1][17]=19; TT[1][18]=20; TT[1][19]=27; TT[1][20]=18; TT[1][21]=24; TT[1][22]=5; TT[1][23]=7; TT[1][24]=2; TT[1][25]=22; TT[1][26]=11; TT[1][27]=15; TT[1][28]=25;
    TT[2][0]=21; TT[2][1]=27; TT[2][2]=21; TT[2][3]=19; TT[2][4]=1; TT[2][5]=14; TT[2][6]=5; TT[2][7]=24; TT[2][8]=14; TT[2][9]=6; TT[2][10]=12; TT[2][11]=17; TT[2][12]=9; TT[2][13]=9; TT[2][14]=1; TT[2][15]=0; TT[2][16]=21; TT[2][17]=4; TT[2][18]=20; TT[2][19]=21; TT[2][20]=28; TT[2][21]=18; TT[2][22]=21; TT[2][23]=27; TT[2][24]=17; TT[2][25]=25; TT[2][26]=7; TT[2][27]=16; TT[2][28]=23;

    SS[0][0]=1; SS[0][1]=-1; SS[0][2]=-1; SS[0][3]=1; SS[0][4]=-1; SS[0][5]=1; SS[0][6]=1; SS[0][7]=1; SS[0][8]=-1; SS[0][9]=-1; SS[0][10]=1; SS[0][11]=-1; SS[0][12]=-1; SS[0][13]=-1; SS[0][14]=1; SS[0][15]=1; SS[0][16]=-1; SS[0][17]=-1; SS[0][18]=-1; SS[0][19]=1; SS[0][20]=1; SS[0][21]=1; SS[0][22]=1; SS[0][23]=1; SS[0][24]=1; SS[0][25]=-1; SS[0][26]=-1; SS[0][27]=1; SS[0][28]=1;
    SS[1][0]=1; SS[1][1]=-1; SS[1][2]=-1; SS[1][3]=-1; SS[1][4]=1; SS[1][5]=-1; SS[1][6]=1; SS[1][7]=1; SS[1][8]=1; SS[1][9]=-1; SS[1][10]=-1; SS[1][11]=-1; SS[1][12]=1; SS[1][13]=1; SS[1][14]=1; SS[1][15]=-1; SS[1][16]=1; SS[1][17]=1; SS[1][18]=-1; SS[1][19]=-1; SS[1][20]=-1; SS[1][21]=1; SS[1][22]=-1; SS[1][23]=-1; SS[1][24]=-1; SS[1][25]=-1; SS[1][26]=1; SS[1][27]=-1; SS[1][28]=1;
    SS[2][0]=1; SS[2][1]=-1; SS[2][2]=-1; SS[2][3]=1; SS[2][4]=-1; SS[2][5]=-1; SS[2][6]=-1; SS[2][7]=-1; SS[2][8]=1; SS[2][9]=-1; SS[2][10]=-1; SS[2][11]=-1; SS[2][12]=1; SS[2][13]=1; SS[2][14]=1; SS[2][15]=1; SS[2][16]=1; SS[2][17]=-1; SS[2][18]=-1; SS[2][19]=-1; SS[2][20]=-1; SS[2][21]=-1; SS[2][22]=1; SS[2][23]=-1; SS[2][24]=1; SS[2][25]=1; SS[2][26]=-1; SS[2][27]=-1; SS[2][28]=1;
#endif
#if WINDOW==5
    TT[0][0]=3; TT[0][1]=7; TT[0][2]=6; TT[0][3]=2; TT[0][4]=6; TT[0][5]=5; TT[0][6]=12; TT[0][7]=15; TT[0][8]=14; TT[0][9]=14; TT[0][10]=4; TT[0][11]=14; TT[0][12]=13; TT[0][13]=11; TT[0][14]=13; TT[0][15]=13; TT[0][16]=7; TT[0][17]=10; TT[0][18]=11; TT[0][19]=10; TT[0][20]=15; TT[0][21]=7; TT[0][22]=2; TT[0][23]=8; TT[0][24]=11; TT[0][25]=5; TT[0][26]=12; TT[0][27]=7; TT[0][28]=11; TT[0][29]=10; TT[0][30]=13; TT[0][31]=10; TT[0][32]=4; TT[0][33]=1; TT[0][34]=0;
    TT[1][0]=13; TT[1][1]=13; TT[1][2]=12; TT[1][3]=1; TT[1][4]=9; TT[1][5]=0; TT[1][6]=8; TT[1][7]=1; TT[1][8]=6; TT[1][9]=7; TT[1][10]=6; TT[1][11]=7; TT[1][12]=7; TT[1][13]=6; TT[1][14]=4; TT[1][15]=2; TT[1][16]=5; TT[1][17]=7; TT[1][18]=6; TT[1][19]=14; TT[1][20]=11; TT[1][21]=11; TT[1][22]=5; TT[1][23]=8; TT[1][24]=1; TT[1][25]=6; TT[1][26]=8; TT[1][27]=15; TT[1][28]=8; TT[1][29]=15; TT[1][30]=14; TT[1][31]=12; TT[1][32]=6; TT[1][33]=7; TT[1][34]=9;
    TT[2][0]=1; TT[2][1]=5; TT[2][2]=14; TT[2][3]=10; TT[2][4]=12; TT[2][5]=0; TT[2][6]=4; TT[2][7]=5; TT[2][8]=7; TT[2][9]=7; TT[2][10]=13; TT[2][11]=13; TT[2][12]=1; TT[2][13]=6; TT[2][14]=11; TT[2][15]=3; TT[2][16]=3; TT[2][17]=3; TT[2][18]=4; TT[2][19]=14; TT[2][20]=12; TT[2][21]=7; TT[2][22]=13; TT[2][23]=8; TT[2][24]=9; TT[2][25]=8; TT[2][26]=14; TT[2][27]=4; TT[2][28]=10; TT[2][29]=3; TT[2][30]=9; TT[2][31]=4; TT[2][32]=7; TT[2][33]=1; TT[2][34]=7;

    SS[0][0]=1; SS[0][1]=-1; SS[0][2]=-1; SS[0][3]=1; SS[0][4]=-1; SS[0][5]=1; SS[0][6]=1; SS[0][7]=1; SS[0][8]=-1; SS[0][9]=-1; SS[0][10]=1; SS[0][11]=-1; SS[0][12]=-1; SS[0][13]=-1; SS[0][14]=1; SS[0][15]=1; SS[0][16]=-1; SS[0][17]=-1; SS[0][18]=-1; SS[0][19]=1; SS[0][20]=1; SS[0][21]=1; SS[0][22]=1; SS[0][23]=1; SS[0][24]=1; SS[0][25]=-1; SS[0][26]=-1; SS[0][27]=1; SS[0][28]=1; SS[0][29]=1; SS[0][30]=-1; SS[0][31]=-1; SS[0][32]=-1; SS[0][33]=1; SS[0][34]=-1;
    SS[1][0]=1; SS[1][1]=1; SS[1][2]=1; SS[1][3]=-1; SS[1][4]=-1; SS[1][5]=-1; SS[1][6]=1; SS[1][7]=1; SS[1][8]=1; SS[1][9]=-1; SS[1][10]=1; SS[1][11]=1; SS[1][12]=-1; SS[1][13]=-1; SS[1][14]=-1; SS[1][15]=1; SS[1][16]=-1; SS[1][17]=-1; SS[1][18]=-1; SS[1][19]=-1; SS[1][20]=1; SS[1][21]=-1; SS[1][22]=1; SS[1][23]=1; SS[1][24]=-1; SS[1][25]=-1; SS[1][26]=1; SS[1][27]=-1; SS[1][28]=-1; SS[1][29]=-1; SS[1][30]=-1; SS[1][31]=1; SS[1][32]=-1; SS[1][33]=-1; SS[1][34]=-1;
    SS[2][0]=1; SS[2][1]=1; SS[2][2]=1; SS[2][3]=1; SS[2][4]=1; SS[2][5]=-1; SS[2][6]=-1; SS[2][7]=-1; SS[2][8]=-1; SS[2][9]=-1; SS[2][10]=1; SS[2][11]=-1; SS[2][12]=1; SS[2][13]=1; SS[2][14]=-1; SS[2][15]=-1; SS[2][16]=-1; SS[2][17]=1; SS[2][18]=1; SS[2][19]=-1; SS[2][20]=-1; SS[2][21]=1; SS[2][22]=-1; SS[2][23]=-1; SS[2][24]=-1; SS[2][25]=1; SS[2][26]=1; SS[2][27]=1; SS[2][28]=1; SS[2][29]=-1; SS[2][30]=1; SS[2][31]=1; SS[2][32]=1; SS[2][33]=-1; SS[2][34]=1;
    SS[2][0]=1; SS[2][1]=1; SS[2][2]=1; SS[2][3]=1; SS[2][4]=1; SS[2][5]=-1; SS[2][6]=-1; SS[2][7]=-1; SS[2][8]=-1; SS[2][9]=-1; SS[2][10]=1; SS[2][11]=-1; SS[2][12]=1; SS[2][13]=1; SS[2][14]=-1; SS[2][15]=-1; SS[2][16]=-1; SS[2][17]=1; SS[2][18]=1; SS[2][19]=-1; SS[2][20]=-1; SS[2][21]=1; SS[2][22]=-1; SS[2][23]=-1; SS[2][24]=-1; SS[2][25]=1; SS[2][26]=1; SS[2][27]=1; SS[2][28]=1; SS[2][29]=-1; SS[2][30]=1; SS[2][31]=1; SS[2][32]=1; SS[2][33]=-1; SS[2][34]=1;
#endif
#if WINDOW==4
    //Use these values of TT[] and SS[] to compare the result to that of the variable-base
    /*TT[0][0]=5; TT[0][1]=6; TT[0][2]=7; TT[0][3]=3; TT[0][4]=3; TT[0][5]=3; TT[0][6]=2; TT[0][7]=0; TT[0][8]=3; TT[0][9]=3; TT[0][10]=3; TT[0][11]=3; TT[0][12]=6; TT[0][13]=4; TT[0][14]=6; TT[0][15]=3; TT[0][16]=2; TT[0][17]=3; TT[0][18]=7; TT[0][19]=4; TT[0][20]=4; TT[0][21]=7; TT[0][22]=2; TT[0][23]=1; TT[0][24]=5; TT[0][25]=7; TT[0][26]=3; TT[0][27]=7; TT[0][28]=7; TT[0][29]=4; TT[0][30]=5; TT[0][31]=5; TT[0][32]=2; TT[0][33]=2; TT[0][34]=6; TT[0][35]=2; TT[0][36]=4; TT[0][37]=6; TT[0][38]=5; TT[0][39]=2; TT[0][40]=2; TT[0][41]=0; TT[0][42]=0; TT[0][43]=4;
    TT[1][0]=0; TT[1][1]=7; TT[1][2]=7; TT[1][3]=4; TT[1][4]=4; TT[1][5]=2; TT[1][6]=1; TT[1][7]=5; TT[1][8]=5; TT[1][9]=2; TT[1][10]=0; TT[1][11]=2; TT[1][12]=5; TT[1][13]=2; TT[1][14]=0; TT[1][15]=2; TT[1][16]=1; TT[1][17]=4; TT[1][18]=5; TT[1][19]=7; TT[1][20]=3; TT[1][21]=4; TT[1][22]=7; TT[1][23]=1; TT[1][24]=2; TT[1][25]=7; TT[1][26]=5; TT[1][27]=1; TT[1][28]=3; TT[1][29]=0; TT[1][30]=2; TT[1][31]=7; TT[1][32]=2; TT[1][33]=2; TT[1][34]=7; TT[1][35]=2; TT[1][36]=2; TT[1][37]=6; TT[1][38]=6; TT[1][39]=6; TT[1][40]=2; TT[1][41]=6; TT[1][42]=1; TT[1][43]=3;
    TT[2][0]=3; TT[2][1]=5; TT[2][2]=4; TT[2][3]=2; TT[2][4]=5; TT[2][5]=6; TT[2][6]=6; TT[2][7]=6; TT[2][8]=5; TT[2][9]=4; TT[2][10]=6; TT[2][11]=3; TT[2][12]=0; TT[2][13]=5; TT[2][14]=6; TT[2][15]=6; TT[2][16]=7; TT[2][17]=4; TT[2][18]=2; TT[2][19]=4; TT[2][20]=2; TT[2][21]=5; TT[2][22]=4; TT[2][23]=4; TT[2][24]=4; TT[2][25]=7; TT[2][26]=1; TT[2][27]=0; TT[2][28]=4; TT[2][29]=4; TT[2][30]=5; TT[2][31]=5; TT[2][32]=1; TT[2][33]=3; TT[2][34]=3; TT[2][35]=1; TT[2][36]=6; TT[2][37]=1; TT[2][38]=2; TT[2][39]=1; TT[2][40]=2; TT[2][41]=3; TT[2][42]=3; TT[2][43]=1;

    SS[0][0]=-1; SS[0][1]=1; SS[0][2]=1; SS[0][3]=-1; SS[0][4]=-1; SS[0][5]=1; SS[0][6]=-1; SS[0][7]=1; SS[0][8]=1; SS[0][9]=1; SS[0][10]=-1; SS[0][11]=1; SS[0][12]=1; SS[0][13]=1; SS[0][14]=-1; SS[0][15]=-1; SS[0][16]=1; SS[0][17]=-1; SS[0][18]=-1; SS[0][19]=-1; SS[0][20]=1; SS[0][21]=1; SS[0][22]=-1; SS[0][23]=-1; SS[0][24]=-1; SS[0][25]=-1; SS[0][26]=1; SS[0][27]=-1; SS[0][28]=-1; SS[0][29]=-1; SS[0][30]=-1; SS[0][31]=1; SS[0][32]=-1; SS[0][33]=1; SS[0][34]=1; SS[0][35]=1; SS[0][36]=-1; SS[0][37]=1; SS[0][38]=1; SS[0][39]=1; SS[0][40]=-1; SS[0][41]=-1; SS[0][42]=1; SS[0][43]=1;
    SS[1][0]=1; SS[1][1]=1; SS[1][2]=1; SS[1][3]=1; SS[1][4]=1; SS[1][5]=-1; SS[1][6]=1; SS[1][7]=-1; SS[1][8]=-1; SS[1][9]=-1; SS[1][10]=1; SS[1][11]=1; SS[1][12]=-1; SS[1][13]=1; SS[1][14]=1; SS[1][15]=-1; SS[1][16]=1; SS[1][17]=-1; SS[1][18]=1; SS[1][19]=-1; SS[1][20]=1; SS[1][21]=1; SS[1][22]=-1; SS[1][23]=-1; SS[1][24]=-1; SS[1][25]=1; SS[1][26]=1; SS[1][27]=-1; SS[1][28]=-1; SS[1][29]=1; SS[1][30]=1; SS[1][31]=1; SS[1][32]=1; SS[1][33]=1; SS[1][34]=-1; SS[1][35]=-1; SS[1][36]=-1; SS[1][37]=-1; SS[1][38]=-1; SS[1][39]=-1; SS[1][40]=1; SS[1][41]=1; SS[1][42]=1; SS[1][43]=1;
    SS[2][0]=-1; SS[2][1]=-1; SS[2][2]=1; SS[2][3]=-1; SS[2][4]=-1; SS[2][5]=1; SS[2][6]=1; SS[2][7]=-1; SS[2][8]=-1; SS[2][9]=-1; SS[2][10]=-1; SS[2][11]=-1; SS[2][12]=1; SS[2][13]=1; SS[2][14]=-1; SS[2][15]=-1; SS[2][16]=-1; SS[2][17]=1; SS[2][18]=-1; SS[2][19]=-1; SS[2][20]=1; SS[2][21]=1; SS[2][22]=1; SS[2][23]=1; SS[2][24]=-1; SS[2][25]=-1; SS[2][26]=-1; SS[2][27]=-1; SS[2][28]=1; SS[2][29]=1; SS[2][30]=1; SS[2][31]=-1; SS[2][32]=1; SS[2][33]=1; SS[2][34]=-1; SS[2][35]=-1; SS[2][36]=1; SS[2][37]=-1; SS[2][38]=-1; SS[2][39]=1; SS[2][40]=-1; SS[2][41]=1; SS[2][42]=-1; SS[2][43]=1;*/


    TT[0][0]=4; TT[0][1]=5; TT[0][2]=7; TT[0][3]=4; TT[0][4]=0; TT[0][5]=0; TT[0][6]=3; TT[0][7]=6; TT[0][8]=5; TT[0][9]=6; TT[0][10]=7; TT[0][11]=5; TT[0][12]=3; TT[0][13]=2; TT[0][14]=2; TT[0][15]=7; TT[0][16]=6; TT[0][17]=5; TT[0][18]=6; TT[0][19]=5; TT[0][20]=7; TT[0][21]=5; TT[0][22]=1; TT[0][23]=5; TT[0][24]=2; TT[0][25]=0; TT[0][26]=3; TT[0][27]=0; TT[0][28]=1; TT[0][29]=3; TT[0][30]=6; TT[0][31]=2; TT[0][32]=3; TT[0][33]=6; TT[0][34]=4; TT[0][35]=7; TT[0][36]=3; TT[0][37]=0; TT[0][38]=1; TT[0][39]=4; TT[0][40]=4; TT[0][41]=1; TT[0][42]=6; TT[0][43]=0;
    TT[1][0]=7; TT[1][1]=5; TT[1][2]=0; TT[1][3]=4; TT[1][4]=2; TT[1][5]=6; TT[1][6]=1; TT[1][7]=2; TT[1][8]=4; TT[1][9]=2; TT[1][10]=1; TT[1][11]=7; TT[1][12]=4; TT[1][13]=3; TT[1][14]=2; TT[1][15]=6; TT[1][16]=2; TT[1][17]=2; TT[1][18]=6; TT[1][19]=2; TT[1][20]=1; TT[1][21]=6; TT[1][22]=4; TT[1][23]=1; TT[1][24]=0; TT[1][25]=0; TT[1][26]=0; TT[1][27]=2; TT[1][28]=3; TT[1][29]=0; TT[1][30]=0; TT[1][31]=1; TT[1][32]=2; TT[1][33]=3; TT[1][34]=6; TT[1][35]=6; TT[1][36]=3; TT[1][37]=2; TT[1][38]=5; TT[1][39]=3; TT[1][40]=3; TT[1][41]=5; TT[1][42]=2; TT[1][43]=5;
    TT[2][0]=1; TT[2][1]=2; TT[2][2]=3; TT[2][3]=3; TT[2][4]=1; TT[2][5]=4; TT[2][6]=3; TT[2][7]=1; TT[2][8]=6; TT[2][9]=6; TT[2][10]=7; TT[2][11]=1; TT[2][12]=0; TT[2][13]=3; TT[2][14]=1; TT[2][15]=0; TT[2][16]=4; TT[2][17]=5; TT[2][18]=7; TT[2][19]=7; TT[2][20]=0; TT[2][21]=5; TT[2][22]=0; TT[2][23]=2; TT[2][24]=2; TT[2][25]=6; TT[2][26]=4; TT[2][27]=2; TT[2][28]=6; TT[2][29]=0; TT[2][30]=2; TT[2][31]=2; TT[2][32]=6; TT[2][33]=3; TT[2][34]=5; TT[2][35]=7; TT[2][36]=5; TT[2][37]=6; TT[2][38]=1; TT[2][39]=2; TT[2][40]=3; TT[2][41]=0; TT[2][42]=1; TT[2][43]=0;

    SS[0][0]=1; SS[0][1]=-1; SS[0][2]=-1; SS[0][3]=1; SS[0][4]=-1; SS[0][5]=1; SS[0][6]=1; SS[0][7]=1; SS[0][8]=-1; SS[0][9]=-1; SS[0][10]=1; SS[0][11]=-1; SS[0][12]=-1; SS[0][13]=-1; SS[0][14]=1; SS[0][15]=1; SS[0][16]=-1; SS[0][17]=-1; SS[0][18]=-1; SS[0][19]=1; SS[0][20]=1; SS[0][21]=1; SS[0][22]=1; SS[0][23]=1; SS[0][24]=1; SS[0][25]=-1; SS[0][26]=-1; SS[0][27]=1; SS[0][28]=1; SS[0][29]=1; SS[0][30]=-1; SS[0][31]=-1; SS[0][32]=-1; SS[0][33]=1; SS[0][34]=-1; SS[0][35]=1; SS[0][36]=1; SS[0][37]=1; SS[0][38]=-1; SS[0][39]=-1; SS[0][40]=-1; SS[0][41]=1; SS[0][42]=1; SS[0][43]=1;
    SS[1][0]=-1; SS[1][1]=1; SS[1][2]=1; SS[1][3]=-1; SS[1][4]=-1; SS[1][5]=-1; SS[1][6]=1; SS[1][7]=-1; SS[1][8]=-1; SS[1][9]=-1; SS[1][10]=-1; SS[1][11]=1; SS[1][12]=-1; SS[1][13]=1; SS[1][14]=1; SS[1][15]=-1; SS[1][16]=-1; SS[1][17]=1; SS[1][18]=-1; SS[1][19]=-1; SS[1][20]=-1; SS[1][21]=-1; SS[1][22]=1; SS[1][23]=-1; SS[1][24]=-1; SS[1][25]=-1; SS[1][26]=1; SS[1][27]=1; SS[1][28]=1; SS[1][29]=1; SS[1][30]=1; SS[1][31]=-1; SS[1][32]=-1; SS[1][33]=-1; SS[1][34]=-1; SS[1][35]=-1; SS[1][36]=1; SS[1][37]=-1; SS[1][38]=1; SS[1][39]=1; SS[1][40]=-1; SS[1][41]=-1; SS[1][42]=-1; SS[1][43]=1;
    SS[2][0]=1; SS[2][1]=-1; SS[2][2]=-1; SS[2][3]=1; SS[2][4]=-1; SS[2][5]=-1; SS[2][6]=-1; SS[2][7]=1; SS[2][8]=1; SS[2][9]=1; SS[2][10]=1; SS[2][11]=-1; SS[2][12]=1; SS[2][13]=1; SS[2][14]=1; SS[2][15]=-1; SS[2][16]=1; SS[2][17]=1; SS[2][18]=1; SS[2][19]=1; SS[2][20]=1; SS[2][21]=-1; SS[2][22]=-1; SS[2][23]=1; SS[2][24]=-1; SS[2][25]=-1; SS[2][26]=-1; SS[2][27]=-1; SS[2][28]=1; SS[2][29]=-1; SS[2][30]=-1; SS[2][31]=1; SS[2][32]=1; SS[2][33]=1; SS[2][34]=-1; SS[2][35]=1; SS[2][36]=-1; SS[2][37]=1; SS[2][38]=-1; SS[2][39]=-1; SS[2][40]=1; SS[2][41]=1; SS[2][42]=1; SS[2][43]=1;
#endif
#else
#if WINDOW==6
    TT[0][0]=1; TT[0][1]=0; TT[0][2]=1; TT[0][3]=0; TT[0][4]=3; TT[0][5]=1; TT[0][6]=3; TT[0][7]=3; TT[0][8]=0; TT[0][9]=1; TT[0][10]=0; TT[0][11]=0; TT[0][12]=1; TT[0][13]=1; TT[0][14]=1; TT[0][15]=1; TT[0][16]=0; TT[0][17]=2; TT[0][18]=1; TT[0][19]=0; TT[0][20]=3; TT[0][21]=2; TT[0][22]=2; TT[0][23]=1; TT[0][24]=1; TT[0][25]=2; TT[0][26]=3; TT[0][27]=0; TT[0][28]=1;
    TT[1][0]=2; TT[1][1]=0; TT[1][2]=3; TT[1][3]=2; TT[1][4]=2; TT[1][5]=2; TT[1][6]=3; TT[1][7]=3; TT[1][8]=0; TT[1][9]=2; TT[1][10]=1; TT[1][11]=2; TT[1][12]=0; TT[1][13]=2; TT[1][14]=0; TT[1][15]=0; TT[1][16]=1; TT[1][17]=1; TT[1][18]=1; TT[1][19]=3; TT[1][20]=3; TT[1][21]=2; TT[1][22]=3; TT[1][23]=3; TT[1][24]=0; TT[1][25]=1; TT[1][26]=1; TT[1][27]=1; TT[1][28]=3;
    TT[2][0]=1; TT[2][1]=0; TT[2][2]=1; TT[2][3]=2; TT[2][4]=2; TT[2][5]=0; TT[2][6]=0; TT[2][7]=3; TT[2][8]=1; TT[2][9]=3; TT[2][10]=1; TT[2][11]=3; TT[2][12]=2; TT[2][13]=0; TT[2][14]=0; TT[2][15]=3; TT[2][16]=2; TT[2][17]=2; TT[2][18]=3; TT[2][19]=3; TT[2][20]=2; TT[2][21]=3; TT[2][22]=1; TT[2][23]=2; TT[2][24]=1; TT[2][25]=2; TT[2][26]=0; TT[2][27]=2; TT[2][28]=16;

    SS[0][0]=-1; SS[0][1]=-1; SS[0][2]=1; SS[0][3]=-1; SS[0][4]=-1; SS[0][5]=-1; SS[0][6]=-1; SS[0][7]=-1; SS[0][8]=-1; SS[0][9]=1; SS[0][10]=-1; SS[0][11]=-1; SS[0][12]=1; SS[0][13]=1; SS[0][14]=-1; SS[0][15]=-1; SS[0][16]=-1; SS[0][17]=-1; SS[0][18]=1; SS[0][19]=1; SS[0][20]=1; SS[0][21]=-1; SS[0][22]=-1; SS[0][23]=1; SS[0][24]=-1; SS[0][25]=-1; SS[0][26]=-1; SS[0][27]=1; SS[0][28]=-1;
    SS[1][0]=-1; SS[1][1]=1; SS[1][2]=-1; SS[1][3]=1; SS[1][4]=1; SS[1][5]=1; SS[1][6]=1; SS[1][7]=-1; SS[1][8]=-1; SS[1][9]=-1; SS[1][10]=1; SS[1][11]=1; SS[1][12]=1; SS[1][13]=-1; SS[1][14]=1; SS[1][15]=1; SS[1][16]=-1; SS[1][17]=1; SS[1][18]=1; SS[1][19]=1; SS[1][20]=1; SS[1][21]=1; SS[1][22]=-1; SS[1][23]=1; SS[1][24]=1; SS[1][25]=-1; SS[1][26]=1; SS[1][27]=1; SS[1][28]=-1;
    SS[2][0]=1; SS[2][1]=1; SS[2][2]=1; SS[2][3]=-1; SS[2][4]=1; SS[2][5]=-1; SS[2][6]=1; SS[2][7]=1; SS[2][8]=1; SS[2][9]=-1; SS[2][10]=1; SS[2][11]=-1; SS[2][12]=1; SS[2][13]=1; SS[2][14]=1; SS[2][15]=1; SS[2][16]=-1; SS[2][17]=-1; SS[2][18]=-1; SS[2][19]=1; SS[2][20]=-1; SS[2][21]=-1; SS[2][22]=-1; SS[2][23]=1; SS[2][24]=1; SS[2][25]=1; SS[2][26]=-1; SS[2][27]=-1; SS[2][28]=1;
#endif
#if WINDOW==5
    TT[0][0]=2; TT[0][1]=0; TT[0][2]=3; TT[0][3]=2; TT[0][4]=2; TT[0][5]=3; TT[0][6]=0; TT[0][7]=0; TT[0][8]=1; TT[0][9]=0; TT[0][10]=1; TT[0][11]=0; TT[0][12]=2; TT[0][13]=3; TT[0][14]=3; TT[0][15]=0; TT[0][16]=2; TT[0][17]=3; TT[0][18]=2; TT[0][19]=0; TT[0][20]=0; TT[0][21]=3; TT[0][22]=1; TT[0][23]=1; TT[0][24]=0; TT[0][25]=2; TT[0][26]=2; TT[0][27]=3; TT[0][28]=0; TT[0][29]=3; TT[0][30]=2; TT[0][31]=2; TT[0][32]=1; TT[0][33]=3; TT[0][34]=0;
    TT[1][0]=0; TT[1][1]=1; TT[1][2]=3; TT[1][3]=0; TT[1][4]=2; TT[1][5]=1; TT[1][6]=0; TT[1][7]=3; TT[1][8]=3; TT[1][9]=2; TT[1][10]=0; TT[1][11]=0; TT[1][12]=3; TT[1][13]=1; TT[1][14]=3; TT[1][15]=2; TT[1][16]=3; TT[1][17]=0; TT[1][18]=0; TT[1][19]=0; TT[1][20]=1; TT[1][21]=0; TT[1][22]=0; TT[1][23]=1; TT[1][24]=0; TT[1][25]=0; TT[1][26]=1; TT[1][27]=1; TT[1][28]=1; TT[1][29]=1; TT[1][30]=0; TT[1][31]=0; TT[1][32]=0; TT[1][33]=0; TT[1][34]=0;
    TT[2][0]=0; TT[2][1]=0; TT[2][2]=0; TT[2][3]=1; TT[2][4]=1; TT[2][5]=1; TT[2][6]=1; TT[2][7]=0; TT[2][8]=0; TT[2][9]=0; TT[2][10]=0; TT[2][11]=0; TT[2][12]=0; TT[2][13]=0; TT[2][14]=0; TT[2][15]=0; TT[2][16]=1; TT[2][17]=1; TT[2][18]=1; TT[2][19]=1; TT[2][20]=0; TT[2][21]=1; TT[2][22]=0; TT[2][23]=0; TT[2][24]=1; TT[2][25]=0; TT[2][26]=1; TT[2][27]=0; TT[2][28]=1; TT[2][29]=1; TT[2][30]=0; TT[2][31]=8; TT[2][32]=9; TT[2][33]=1; TT[2][34]=1;

    SS[0][0]=-1; SS[0][1]=-1; SS[0][2]=1; SS[0][3]=-1; SS[0][4]=-1; SS[0][5]=-1; SS[0][6]=-1; SS[0][7]=-1; SS[0][8]=-1; SS[0][9]=1; SS[0][10]=-1; SS[0][11]=-1; SS[0][12]=1; SS[0][13]=1; SS[0][14]=-1; SS[0][15]=-1; SS[0][16]=-1; SS[0][17]=-1; SS[0][18]=1; SS[0][19]=1; SS[0][20]=1; SS[0][21]=-1; SS[0][22]=-1; SS[0][23]=1; SS[0][24]=-1; SS[0][25]=-1; SS[0][26]=-1; SS[0][27]=1; SS[0][28]=-1; SS[0][29]=-1; SS[0][30]=1; SS[0][31]=-1; SS[0][32]=1; SS[0][33]=1; SS[0][34]=1;
    SS[1][0]=1; SS[1][1]=-1; SS[1][2]=-1; SS[1][3]=-1; SS[1][4]=1; SS[1][5]=1; SS[1][6]=1; SS[1][7]=-1; SS[1][8]=1; SS[1][9]=1; SS[1][10]=-1; SS[1][11]=1; SS[1][12]=1; SS[1][13]=1; SS[1][14]=1; SS[1][15]=1; SS[1][16]=-1; SS[1][17]=1; SS[1][18]=1; SS[1][19]=-1; SS[1][20]=1; SS[1][21]=1; SS[1][22]=-1; SS[1][23]=1; SS[1][24]=1; SS[1][25]=1; SS[1][26]=-1; SS[1][27]=1; SS[1][28]=-1; SS[1][29]=1; SS[1][30]=1; SS[1][31]=1; SS[1][32]=-1; SS[1][33]=1; SS[1][34]=-1;
    SS[2][0]=1; SS[2][1]=1; SS[2][2]=1; SS[2][3]=1; SS[2][4]=-1; SS[2][5]=-1; SS[2][6]=-1; SS[2][7]=1; SS[2][8]=-1; SS[2][9]=-1; SS[2][10]=-1; SS[2][11]=1; SS[2][12]=1; SS[2][13]=1; SS[2][14]=-1; SS[2][15]=-1; SS[2][16]=1; SS[2][17]=1; SS[2][18]=-1; SS[2][19]=-1; SS[2][20]=1; SS[2][21]=-1; SS[2][22]=-1; SS[2][23]=-1; SS[2][24]=1; SS[2][25]=-1; SS[2][26]=-1; SS[2][27]=-1; SS[2][28]=1; SS[2][29]=1; SS[2][30]=1; SS[2][31]=-1; SS[2][32]=1; SS[2][33]=1; SS[2][34]=1;
#endif
#if WINDOW==4
    TT[0][0]=1; TT[0][1]=1; TT[0][2]=0; TT[0][3]=0; TT[0][4]=0; TT[0][5]=1; TT[0][6]=0; TT[0][7]=1; TT[0][8]=1; TT[0][9]=0; TT[0][10]=1; TT[0][11]=0; TT[0][12]=0; TT[0][13]=1; TT[0][14]=0; TT[0][15]=1; TT[0][16]=1; TT[0][17]=1; TT[0][18]=1; TT[0][19]=0; TT[0][20]=1; TT[0][21]=1; TT[0][22]=0; TT[0][23]=1; TT[0][24]=1; TT[0][25]=0; TT[0][26]=0; TT[0][27]=0; TT[0][28]=1; TT[0][29]=1; TT[0][30]=1; TT[0][31]=1; TT[0][32]=1; TT[0][33]=0; TT[0][34]=1; TT[0][35]=0; TT[0][36]=1; TT[0][37]=1; TT[0][38]=1; TT[0][39]=1; TT[0][40]=0; TT[0][41]=0; TT[0][42]=0; TT[0][43]=0;
    TT[1][0]=0; TT[1][1]=0; TT[1][2]=1; TT[1][3]=1; TT[1][4]=0; TT[1][5]=0; TT[1][6]=1; TT[1][7]=1; TT[1][8]=0; TT[1][9]=0; TT[1][10]=0; TT[1][11]=0; TT[1][12]=0; TT[1][13]=0; TT[1][14]=0; TT[1][15]=1; TT[1][16]=1; TT[1][17]=1; TT[1][18]=1; TT[1][19]=1; TT[1][20]=1; TT[1][21]=1; TT[1][22]=1; TT[1][23]=0; TT[1][24]=0; TT[1][25]=1; TT[1][26]=0; TT[1][27]=1; TT[1][28]=1; TT[1][29]=0; TT[1][30]=0; TT[1][31]=1; TT[1][32]=0; TT[1][33]=0; TT[1][34]=0; TT[1][35]=0; TT[1][36]=1; TT[1][37]=0; TT[1][38]=1; TT[1][39]=0; TT[1][40]=1; TT[1][41]=0; TT[1][42]=0; TT[1][43]=0;
    TT[2][0]=0; TT[2][1]=0; TT[2][2]=1; TT[2][3]=1; TT[2][4]=0; TT[2][5]=0; TT[2][6]=1; TT[2][7]=0; TT[2][8]=0; TT[2][9]=0; TT[2][10]=0; TT[2][11]=1; TT[2][12]=1; TT[2][13]=1; TT[2][14]=0; TT[2][15]=1; TT[2][16]=0; TT[2][17]=0; TT[2][18]=0; TT[2][19]=1; TT[2][20]=1; TT[2][21]=0; TT[2][22]=0; TT[2][23]=1; TT[2][24]=0; TT[2][25]=0; TT[2][26]=0; TT[2][27]=1; TT[2][28]=1; TT[2][29]=0; TT[2][30]=0; TT[2][31]=0; TT[2][32]=1; TT[2][33]=0; TT[2][34]=1; TT[2][35]=0; TT[2][36]=0; TT[2][37]=5; TT[2][38]=5; TT[2][39]=4; TT[2][40]=4; TT[2][41]=4; TT[2][42]=4; TT[2][43]=4;

    SS[0][0]=-1; SS[0][1]=-1; SS[0][2]=1; SS[0][3]=-1; SS[0][4]=-1; SS[0][5]=-1; SS[0][6]=-1; SS[0][7]=-1; SS[0][8]=-1; SS[0][9]=1; SS[0][10]=-1; SS[0][11]=-1; SS[0][12]=1; SS[0][13]=1; SS[0][14]=-1; SS[0][15]=-1; SS[0][16]=-1; SS[0][17]=-1; SS[0][18]=1; SS[0][19]=1; SS[0][20]=1; SS[0][21]=-1; SS[0][22]=-1; SS[0][23]=1; SS[0][24]=-1; SS[0][25]=-1; SS[0][26]=-1; SS[0][27]=1; SS[0][28]=-1; SS[0][29]=-1; SS[0][30]=1; SS[0][31]=-1; SS[0][32]=1; SS[0][33]=1; SS[0][34]=1; SS[0][35]=1; SS[0][36]=-1; SS[0][37]=-1; SS[0][38]=-1; SS[0][39]=1; SS[0][40]=1; SS[0][41]=1; SS[0][42]=-1; SS[0][43]=1;
    SS[1][0]=1; SS[1][1]=-1; SS[1][2]=1; SS[1][3]=1; SS[1][4]=1; SS[1][5]=1; SS[1][6]=1; SS[1][7]=-1; SS[1][8]=1; SS[1][9]=1; SS[1][10]=-1; SS[1][11]=1; SS[1][12]=1; SS[1][13]=-1; SS[1][14]=1; SS[1][15]=1; SS[1][16]=1; SS[1][17]=-1; SS[1][18]=1; SS[1][19]=-1; SS[1][20]=1; SS[1][21]=1; SS[1][22]=1; SS[1][23]=-1; SS[1][24]=1; SS[1][25]=-1; SS[1][26]=1; SS[1][27]=1; SS[1][28]=1; SS[1][29]=1; SS[1][30]=-1; SS[1][31]=-1; SS[1][32]=-1; SS[1][33]=1; SS[1][34]=-1; SS[1][35]=-1; SS[1][36]=-1; SS[1][37]=1; SS[1][38]=1; SS[1][39]=1; SS[1][40]=-1; SS[1][41]=-1; SS[1][42]=1; SS[1][43]=1;
    SS[2][0]=-1; SS[2][1]=-1; SS[2][2]=1; SS[2][3]=-1; SS[2][4]=-1; SS[2][5]=-1; SS[2][6]=1; SS[2][7]=-1; SS[2][8]=-1; SS[2][9]=-1; SS[2][10]=1; SS[2][11]=1; SS[2][12]=1; SS[2][13]=-1; SS[2][14]=1; SS[2][15]=1; SS[2][16]=-1; SS[2][17]=-1; SS[2][18]=1; SS[2][19]=-1; SS[2][20]=-1; SS[2][21]=1; SS[2][22]=1; SS[2][23]=1; SS[2][24]=-1; SS[2][25]=1; SS[2][26]=-1; SS[2][27]=1; SS[2][28]=1; SS[2][29]=-1; SS[2][30]=1; SS[2][31]=1; SS[2][32]=1; SS[2][33]=-1; SS[2][34]=1; SS[2][35]=1; SS[2][36]=1; SS[2][37]=-1; SS[2][38]=-1; SS[2][39]=-1; SS[2][40]=-1; SS[2][41]=-1; SS[2][42]=-1; SS[2][43]=1;
#endif
#endif

    Offline_precomp(&JP, PreT);

	for(i=0;i<10;i++){
		bef=rdtsc();
		for (j=0;j<lpz;j++){
			init(&JP);
			scalar_mul(SS, TT, PreT, &JP);
		}
		aft=rdtscp();
		if((aft-bef)/(lpz) < min_ccycle)
				min_ccycle= (aft-bef)/(lpz);
	}
	printf("Window width :: %u\n%s%lu\n", WINDOW, "The minimum clock cycles count is :: ", min_ccycle);
	output(&JP);

}
