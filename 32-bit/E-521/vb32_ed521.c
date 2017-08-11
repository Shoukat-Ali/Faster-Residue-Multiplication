/* Test program for ed521.c constant-time variable-base scalar multiplication
 This curve is recommended by Daniel J. Bernstein and Tanja Lange at http://safecurves.cr.yp.to
 Curve Equation :: x^2+y^2 = 1-376014x^2y^2 where d=-376014
 Uses Bernstein et al.s point multiplication method from Curve41417 paper
 Cache safety thanks to ed25519
 Thanks to M. Scott for making his code public and is available at http://indigo.ie/~mscott/ed521.cpp
 We have replaced gmul() by TMV_product()
 We have added rdtscp() to measure the clock cycles count as suggested by Paoloni; http://www.intel.com.tr/content/dam/www/public/us/en/documents/white-papers/ia-32-ia-64-benchmark-code-execution-paper.pdf
 We have amended scr(), gmuli(), gsqr(), and gsqr2() according to modulus p
 gcc -Wall -m32 -O3 vb32_ed521.c -o vb32_ed521.exe */


#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#define WINDOW 4 //5 //6

#if WINDOW==4
#define PANES 131
#endif

#if WINDOW==5
#define PANES 105
#endif

#if WINDOW==6
#define PANES 87
#endif


#define LIMBS 18
#define M (1<<(WINDOW-1))
#define CURVE_PARAMETER 376014

//#define TEST  //define to multiply by group order


static const int32_t lower29bits = 0x1fffffff;
static const int32_t lower28bits = 0xfffffff;


__inline__ uint64_t rdtsc(){
   uint32_t lo, hi;
   __asm__ __volatile__ ("xorl %%eax,%%eax \n        cpuid" ::: "%eax", "%ebx", "%ecx", "%edx");
   __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));

   return (uint64_t)hi << 32 | lo;
}

__inline__ uint64_t rdtscp(){
   uint32_t lo, hi;
   __asm__ __volatile__ ("rdtscp" : "=a" (lo), "=d" (hi) :: "%ecx");
   __asm__ __volatile__ ("cpuid" ::: "%eax", "%ebx","%ecx", "%edx");
   return (uint64_t)hi << 32 | lo;
}

///////////////////////// Field Operations ////////////////////////////////////
/*Because of the limbs sizes, Point Arithmetic formulas and efficiency purpose we are using the strategy of in-place reduction
because both from the worst-case analysis of the size and the testing result, the inputs to multiplication and squaring routines should
be in the range of scr otherwise, overflow is resulted and the result will be incorrect.*/

//w=x+y and the result is reduced in-place
void gadd(int32_t x[], int32_t y[], int32_t w[]){
	int32_t t;
	t=x[0]+y[0];                 //Could be at most 31-bit
	w[0]=t&lower29bits;

	t=(x[1]+y[1])+(t>>29);              w[1]=t&lower29bits;
	t=(x[2]+y[2])+(t>>29);              w[2]=t&lower29bits;
	t=(x[3]+y[3])+(t>>29);              w[3]=t&lower29bits;
	t=(x[4]+y[4])+(t>>29);              w[4]=t&lower29bits;
	t=(x[5]+y[5])+(t>>29);              w[5]=t&lower29bits;
	t=(x[6]+y[6])+(t>>29);              w[6]=t&lower29bits;
	t=(x[7]+y[7])+(t>>29);              w[7]=t&lower29bits;
	t=(x[8]+y[8])+(t>>29);              w[8]=t&lower29bits;
	t=(x[9]+y[9])+(t>>29);              w[9]=t&lower29bits;
	t=(x[10]+y[10])+(t>>29);            w[10]=t&lower29bits;
	t=(x[11]+y[11])+(t>>29);            w[11]=t&lower29bits;
	t=(x[12]+y[12])+(t>>29);            w[12]=t&lower29bits;
	t=(x[13]+y[13])+(t>>29);            w[13]=t&lower29bits;
	t=(x[14]+y[14])+(t>>29);            w[14]=t&lower29bits;
	t=(x[15]+y[15])+(t>>29);            w[15]=t&lower29bits;
	t=(x[16]+y[16])+(t>>29);            w[16]=t&lower29bits;
	t=(x[17]+y[17])+(t>>29);            w[17]=t&lower28bits;               //Most significant limb is at most 28-bit

	w[0]+=(t>>28);            //Least significant limb could be at most 30-bit

}

// w=x-y and the result is reduced in-place
void gsub(int32_t x[], int32_t y[], int32_t w[]){
    int32_t t;
	t=x[0]-y[0];                 //Could be at most 31-bit
	w[0]=t&lower29bits;

	t=(x[1]-y[1])+(t>>29);              w[1]=t&lower29bits;
	t=(x[2]-y[2])+(t>>29);              w[2]=t&lower29bits;
	t=(x[3]-y[3])+(t>>29);              w[3]=t&lower29bits;
	t=(x[4]-y[4])+(t>>29);              w[4]=t&lower29bits;
	t=(x[5]-y[5])+(t>>29);              w[5]=t&lower29bits;
	t=(x[6]-y[6])+(t>>29);              w[6]=t&lower29bits;
	t=(x[7]-y[7])+(t>>29);              w[7]=t&lower29bits;
	t=(x[8]-y[8])+(t>>29);              w[8]=t&lower29bits;
	t=(x[9]-y[9])+(t>>29);              w[9]=t&lower29bits;
	t=(x[10]-y[10])+(t>>29);            w[10]=t&lower29bits;
	t=(x[11]-y[11])+(t>>29);            w[11]=t&lower29bits;
	t=(x[12]-y[12])+(t>>29);            w[12]=t&lower29bits;
	t=(x[13]-y[13])+(t>>29);            w[13]=t&lower29bits;
	t=(x[14]-y[14])+(t>>29);            w[14]=t&lower29bits;
	t=(x[15]-y[15])+(t>>29);            w[15]=t&lower29bits;
	t=(x[16]-y[16])+(t>>29);            w[16]=t&lower29bits;
	t=(x[17]-y[17])+(t>>29);            w[17]=t&lower28bits;               //Most significant limb is at most 28-bit

	w[0]+=(t>>28);            //Least significant limb could be at most 30-bit

}

// w=x-2y and the result is reduced in-place
void gsub2(int32_t x[], int32_t y[], int32_t w[]){
    int32_t t;
	t=x[0]-2*y[0];                 //Could be at most 31-bit
	w[0]=t&lower29bits;

	t=(x[1]-2*y[1])+(t>>29);              w[1]=t&lower29bits;
	t=(x[2]-2*y[2])+(t>>29);              w[2]=t&lower29bits;
	t=(x[3]-2*y[3])+(t>>29);              w[3]=t&lower29bits;
	t=(x[4]-2*y[4])+(t>>29);              w[4]=t&lower29bits;
	t=(x[5]-2*y[5])+(t>>29);              w[5]=t&lower29bits;
	t=(x[6]-2*y[6])+(t>>29);              w[6]=t&lower29bits;
	t=(x[7]-2*y[7])+(t>>29);              w[7]=t&lower29bits;
	t=(x[8]-2*y[8])+(t>>29);              w[8]=t&lower29bits;
	t=(x[9]-2*y[9])+(t>>29);              w[9]=t&lower29bits;
	t=(x[10]-2*y[10])+(t>>29);            w[10]=t&lower29bits;
	t=(x[11]-2*y[11])+(t>>29);            w[11]=t&lower29bits;
	t=(x[12]-2*y[12])+(t>>29);            w[12]=t&lower29bits;
	t=(x[13]-2*y[13])+(t>>29);            w[13]=t&lower29bits;
	t=(x[14]-2*y[14])+(t>>29);            w[14]=t&lower29bits;
	t=(x[15]-2*y[15])+(t>>29);            w[15]=t&lower29bits;
	t=(x[16]-2*y[16])+(t>>29);            w[16]=t&lower29bits;
	t=(x[17]-2*y[17])+(t>>29);            w[17]=t&lower28bits;               //Most significant limb is at most 28-bit

	w[0]+=(t>>28);            //Least significant limb could be at most 30-bit

}

// w=x
void gcopy(int32_t x[], int32_t w[]){
	w[0]=x[0];      	w[1]=x[1];          	w[2]=x[2];          	w[3]=x[3];
	w[4]=x[4];      	w[5]=x[5];          	w[6]=x[6];          	w[7]=x[7];
	w[8]=x[8];      	w[9]=x[9];          	w[10]=x[10];        	w[11]=x[11];
	w[12]=x[12];    	w[13]=x[13];        	w[14]=x[14];        	w[15]=x[15];
	w[16]=x[16];        w[17]=x[17];
}


// w*=2 and the result is reduced in-place
void gmul2(int32_t w[]){
    int32_t t;
	t=w[0]*2;		            //Could be at most 31-bit
	w[0]=t&lower29bits;

	t=(w[1]*2)+(t>>29);         w[1]=t&lower29bits;
	t=(w[2]*2)+(t>>29);         w[2]=t&lower29bits;
	t=(w[3]*2)+(t>>29);         w[3]=t&lower29bits;
	t=(w[4]*2)+(t>>29);       	w[4]=t&lower29bits;
	t=(w[5]*2)+(t>>29);       	w[5]=t&lower29bits;
	t=(w[6]*2)+(t>>29);       	w[6]=t&lower29bits;
	t=(w[7]*2)+(t>>29);       	w[7]=t&lower29bits;
	t=(w[8]*2)+(t>>29);       	w[8]=t&lower29bits;
    t=(w[9]*2)+(t>>29);       	w[9]=t&lower29bits;
	t=(w[10]*2)+(t>>29);        w[10]=t&lower29bits;
	t=(w[11]*2)+(t>>29);      	w[11]=t&lower29bits;
	t=(w[12]*2)+(t>>29);      	w[12]=t&lower29bits;
	t=(w[13]*2)+(t>>29);      	w[13]=t&lower29bits;
	t=(w[14]*2)+(t>>29);      	w[14]=t&lower29bits;
	t=(w[15]*2)+(t>>29);      	w[15]=t&lower29bits;
	t=(w[16]*2)+(t>>29);        w[16]=t&lower29bits;
	t=(w[17]*2)+(t>>29);        w[17]=t&lower28bits;           //most significant limb is at most 28-bit

	w[0]+=(t>>28);				//least significant limb could be at most 30-bit

}

// w-=2*x and the result is reduced in-place
/*void gsb2(int32_t x[], int32_t w[]){
    int32_t t;
	t=w[0]-2*x[0];                 //Could be at most 32-bit
	w[0]=t&lower29bits;

	t=(w[1]-2*x[1])+(t>>29);              w[1]=t&lower29bits;
	t=(w[2]-2*x[2])+(t>>29);              w[2]=t&lower29bits;
	t=(w[3]-2*x[3])+(t>>29);              w[3]=t&lower29bits;
	t=(w[4]-2*x[4])+(t>>29);              w[4]=t&lower29bits;
	t=(w[5]-2*x[5])+(t>>29);              w[5]=t&lower29bits;
	t=(w[6]-2*x[6])+(t>>29);              w[6]=t&lower29bits;
	t=(w[7]-2*x[7])+(t>>29);              w[7]=t&lower29bits;
	t=(w[8]-2*x[8])+(t>>29);              w[8]=t&lower29bits;
	t=(w[9]-2*x[9])+(t>>29);              w[9]=t&lower29bits;
	t=(w[10]-2*x[10])+(t>>29);            w[10]=t&lower29bits;
	t=(w[11]-2*x[11])+(t>>29);            w[11]=t&lower29bits;
	t=(w[12]-2*x[12])+(t>>29);            w[12]=t&lower29bits;
	t=(w[13]-2*x[13])+(t>>29);            w[13]=t&lower29bits;
	t=(w[14]-2*x[14])+(t>>29);            w[14]=t&lower29bits;
	t=(w[15]-2*x[15])+(t>>29);            w[15]=t&lower29bits;
	t=(w[16]-2*x[16])+(t>>29);            w[16]=t&lower29bits;
	t=(w[17]-2*x[17])+(t>>29);            w[17]=t&lower28bits;               //Most significant limb is at most 28-bit

	w[0]+=(t>>28);            //Least significant limb could be at most 30-bit

}*/

//w=w-x-y and the result is reduced in-place
void gtsb(int32_t x[], int32_t y[], int32_t w[]){
    int32_t t;
	t=w[0]-x[0]-y[0];                 //Could be at most 32-bit
	w[0]=t&lower29bits;

	t=(w[1]-x[1]-y[1])+(t>>29);              w[1]=t&lower29bits;
	t=(w[2]-x[2]-y[2])+(t>>29);              w[2]=t&lower29bits;
	t=(w[3]-x[3]-y[3])+(t>>29);              w[3]=t&lower29bits;
	t=(w[4]-x[4]-y[4])+(t>>29);              w[4]=t&lower29bits;
	t=(w[5]-x[5]-y[5])+(t>>29);              w[5]=t&lower29bits;
	t=(w[6]-x[6]-y[6])+(t>>29);              w[6]=t&lower29bits;
	t=(w[7]-x[7]-y[7])+(t>>29);              w[7]=t&lower29bits;
	t=(w[8]-x[8]-y[8])+(t>>29);              w[8]=t&lower29bits;
	t=(w[9]-x[9]-y[9])+(t>>29);              w[9]=t&lower29bits;
	t=(w[10]-x[10]-y[10])+(t>>29);            w[10]=t&lower29bits;
	t=(w[11]-x[11]-y[11])+(t>>29);            w[11]=t&lower29bits;
	t=(w[12]-x[12]-y[12])+(t>>29);            w[12]=t&lower29bits;
	t=(w[13]-x[13]-y[13])+(t>>29);            w[13]=t&lower29bits;
	t=(w[14]-x[14]-y[14])+(t>>29);            w[14]=t&lower29bits;
	t=(w[15]-x[15]-y[15])+(t>>29);            w[15]=t&lower29bits;
	t=(w[16]-x[16]-y[16])+(t>>29);            w[16]=t&lower29bits;
	t=(w[17]-x[17]-y[17])+(t>>29);            w[17]=t&lower28bits;               //Most significant limb is at most 28-bit

	w[0]+=(t>>28);            //Least significant limb could be at most 30-bit

}

//Short Coefficient reduction form = scr
void scr(int32_t w[], int32_t d[]){
	int32_t t1;
	d[0]=w[0]&lower29bits;		//In spite of the Reduced limb form w[0] is expected to be 29-bit

	t1=w[1]+(w[0]>>29);         //Arithmetic/Singed right shift operation is performed
	d[1]=t1&lower29bits;

	t1=w[2]+(t1>>29);           d[2]=t1&lower29bits;
	t1=w[3]+(t1>>29);       	d[3]=t1&lower29bits;
	t1=w[4]+(t1>>29);       	d[4]=t1&lower29bits;
	t1=w[5]+(t1>>29);       	d[5]=t1&lower29bits;
	t1=w[6]+(t1>>29);       	d[6]=t1&lower29bits;
	t1=w[7]+(t1>>29);       	d[7]=t1&lower29bits;
	t1=w[8]+(t1>>29);       	d[8]=t1&lower29bits;
	t1=w[9]+(t1>>29);           d[9]=t1&lower29bits;
	t1=w[10]+(t1>>29);      	d[10]=t1&lower29bits;
	t1=w[11]+(t1>>29);      	d[11]=t1&lower29bits;
	t1=w[12]+(t1>>29);      	d[12]=t1&lower29bits;
	t1=w[13]+(t1>>29);      	d[13]=t1&lower29bits;
	t1=w[14]+(t1>>29);         	d[14]=t1&lower29bits;
	t1=w[15]+(t1>>29);      	d[15]=t1&lower29bits;
	t1=w[16]+(t1>>29);      	d[16]=t1&lower29bits;
	t1=w[17]+(t1>>29);          d[17]=t1&lower28bits;           //most significant limb is at most 28-bit

	d[0]=d[0]+(t1>>28);				//least significant limb could be at most 30-bit

}

// multiply w by a constant, w*=i where the product is double-word i.e. 64-bit and the final result result is in scr form
void gmuli(int32_t w[], int32_t i){
	int64_t t;

	t=(int64_t)w[0]*i;                 //Empirical result shows that t could be at most 48-bit
	w[0]=((int32_t)t)&lower29bits;

	t=(int64_t)w[1]*i+(t>>29);         w[1]=((int32_t)t)&lower29bits;
	t=(int64_t)w[2]*i+(t>>29);         w[2]=((int32_t)t)&lower29bits;
	t=(int64_t)w[3]*i+(t>>29);         w[3]=((int32_t)t)&lower29bits;
	t=(int64_t)w[4]*i+(t>>29);         w[4]=((int32_t)t)&lower29bits;
	t=(int64_t)w[5]*i+(t>>29);         w[5]=((int32_t)t)&lower29bits;
	t=(int64_t)w[6]*i+(t>>29);         w[6]=((int32_t)t)&lower29bits;
	t=(int64_t)w[7]*i+(t>>29);         w[7]=((int32_t)t)&lower29bits;
	t=(int64_t)w[8]*i+(t>>29);         w[8]=((int32_t)t)&lower29bits;
	t=(int64_t)w[9]*i+(t>>29);         w[9]=((int32_t)t)&lower29bits;
	t=(int64_t)w[10]*i+(t>>29);         w[10]=((int32_t)t)&lower29bits;
	t=(int64_t)w[11]*i+(t>>29);         w[11]=((int32_t)t)&lower29bits;
	t=(int64_t)w[12]*i+(t>>29);         w[12]=((int32_t)t)&lower29bits;
	t=(int64_t)w[13]*i+(t>>29);         w[13]=((int32_t)t)&lower29bits;
	t=(int64_t)w[14]*i+(t>>29);         w[14]=((int32_t)t)&lower29bits;
	t=(int64_t)w[15]*i+(t>>29);         w[15]=((int32_t)t)&lower29bits;
	t=(int64_t)w[16]*i+(t>>29);         w[16]=((int32_t)t)&lower29bits;
	t=(int64_t)w[17]*i+(t>>29);         w[17]=((int32_t)t)&lower28bits;               //Most significant limb is at most 28-bit

	w[0]+=(int32_t)(t>>28);            //Least significant limb could be at most 30-bit
}

//z=x^2 (mod p)
/*The function computes the Z = X*X mod (2**521) - 1] where X and Z are 18-limb.
X is an array of integers of size 18 where the first element is at most 30-bit, the next 16 elements are 29-bit at most and the last element is 28-bit integer at most.
The output Z_res is one-dimensional array of integers such that Z_res[0] in [0,2**30-1], Z_res[i] in [0,2**29-1] and Z_res[17] in [0,2**28-1]*/
void gsqr(int32_t X[], int32_t Z_res[]){
    uint64_t C = 2*((int64_t)X[0]*X[17]+(int64_t)X[1]*X[16]+(int64_t)X[2]*X[15]+(int64_t)X[3]*X[14]+(int64_t)X[4]*X[13]+(int64_t)X[5]*X[12]+(int64_t)X[6]*X[11]+(int64_t)X[7]*X[10]+(int64_t)X[8]*X[9]);
    uint64_t CC = C&lower28bits;
    C = (int64_t)X[0]*X[0] + 4*((int64_t)X[1]*X[17]+(int64_t)X[2]*X[16]+(int64_t)X[3]*X[15]+(int64_t)X[4]*X[14]+(int64_t)X[5]*X[13]+(int64_t)X[6]*X[12]+(int64_t)X[7]*X[11]+(int64_t)X[8]*X[10])+ 2*((int64_t)X[9]*X[9]) + (C>>28);
    Z_res[0] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[1])+ 4*((int64_t)X[2]*X[17]+(int64_t)X[3]*X[16]+(int64_t)X[4]*X[15]+(int64_t)X[5]*X[14]+(int64_t)X[6]*X[13]+(int64_t)X[7]*X[12]+(int64_t)X[8]*X[11]+(int64_t)X[9]*X[10]) + (C>>29);
    Z_res[1] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[2]+(int64_t)X[10]*X[10]) + (int64_t)X[1]*X[1]+ 4*((int64_t)X[3]*X[17]+(int64_t)X[4]*X[16]+(int64_t)X[5]*X[15]+(int64_t)X[6]*X[14]+(int64_t)X[7]*X[13]+(int64_t)X[8]*X[12]+(int64_t)X[9]*X[11]) + (C>>29);
    Z_res[2] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[3]+(int64_t)X[1]*X[2]) + 4*((int64_t)X[4]*X[17]+(int64_t)X[5]*X[16]+(int64_t)X[6]*X[15]+(int64_t)X[7]*X[14]+(int64_t)X[8]*X[13]+(int64_t)X[9]*X[12]+(int64_t)X[10]*X[11]) + (C>>29);
    Z_res[3] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[4]+(int64_t)X[1]*X[3]+(int64_t)X[11]*X[11]) + (int64_t)X[2]*X[2] + 4*((int64_t)X[5]*X[17]+(int64_t)X[6]*X[16]+(int64_t)X[7]*X[15]+(int64_t)X[8]*X[14]+(int64_t)X[9]*X[13]+(int64_t)X[10]*X[12]) + (C>>29);
    Z_res[4] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[5]+(int64_t)X[1]*X[4]+(int64_t)X[2]*X[3]) + 4*((int64_t)X[6]*X[17]+(int64_t)X[7]*X[16]+(int64_t)X[8]*X[15]+(int64_t)X[9]*X[14]+(int64_t)X[10]*X[13]+(int64_t)X[11]*X[12]) + (C>>29);
    Z_res[5] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[6]+(int64_t)X[1]*X[5]+(int64_t)X[2]*X[4]+(int64_t)X[12]*X[12]) + (int64_t)X[3]*X[3] + 4*((int64_t)X[7]*X[17]+(int64_t)X[8]*X[16]+(int64_t)X[9]*X[15]+(int64_t)X[10]*X[14]+(int64_t)X[11]*X[13]) + (C>>29);
    Z_res[6] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[7]+(int64_t)X[1]*X[6]+(int64_t)X[2]*X[5]+(int64_t)X[3]*X[4]) +4*((int64_t)X[8]*X[17]+(int64_t)X[9]*X[16]+(int64_t)X[10]*X[15]+(int64_t)X[11]*X[14]+(int64_t)X[12]*X[13]) + (C>>29);
    Z_res[7] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[8]+(int64_t)X[1]*X[7]+(int64_t)X[2]*X[6]+(int64_t)X[3]*X[5]+(int64_t)X[13]*X[13])+ (int64_t)X[4]*X[4] +4*((int64_t)X[9]*X[17]+(int64_t)X[10]*X[16]+(int64_t)X[11]*X[15]+(int64_t)X[12]*X[14]) + (C>>29);
    Z_res[8] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[9]+(int64_t)X[1]*X[8]+(int64_t)X[2]*X[7]+(int64_t)X[3]*X[6]+(int64_t)X[4]*X[5]) +4*((int64_t)X[10]*X[17]+(int64_t)X[11]*X[16]+(int64_t)X[12]*X[15]+(int64_t)X[13]*X[14]) + (C>>29);
    Z_res[9] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[10]+(int64_t)X[1]*X[9]+(int64_t)X[2]*X[8]+(int64_t)X[3]*X[7]+(int64_t)X[4]*X[6]+(int64_t)X[14]*X[14])+ (int64_t)X[5]*X[5] +4*((int64_t)X[11]*X[17]+(int64_t)X[12]*X[16]+(int64_t)X[13]*X[15]) + (C>>29);
    Z_res[10] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[11]+(int64_t)X[1]*X[10]+(int64_t)X[2]*X[9]+(int64_t)X[3]*X[8]+(int64_t)X[4]*X[7]+(int64_t)X[5]*X[6]) +4*((int64_t)X[12]*X[17]+(int64_t)X[13]*X[16]+(int64_t)X[14]*X[15]) + (C>>29);
    Z_res[11] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[12]+(int64_t)X[1]*X[11]+(int64_t)X[2]*X[10]+(int64_t)X[3]*X[9]+(int64_t)X[4]*X[8]+(int64_t)X[5]*X[7]+(int64_t)X[15]*X[15]) + (int64_t)X[6]*X[6] +4*((int64_t)X[13]*X[17]+(int64_t)X[14]*X[16]) + (C>>29);
    Z_res[12] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[13]+(int64_t)X[1]*X[12]+(int64_t)X[2]*X[11]+(int64_t)X[3]*X[10]+(int64_t)X[4]*X[9]+(int64_t)X[5]*X[8]+(int64_t)X[6]*X[7]) +4*((int64_t)X[14]*X[17]+(int64_t)X[15]*X[16]) + (C>>29);
    Z_res[13] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[14]+(int64_t)X[1]*X[13]+(int64_t)X[2]*X[12]+(int64_t)X[3]*X[11]+(int64_t)X[4]*X[10]+(int64_t)X[5]*X[9]+(int64_t)X[6]*X[8]+(int64_t)X[16]*X[16]) + (int64_t)X[7]*X[7] +4*((int64_t)X[15]*X[17]) + (C>>29);
    Z_res[14] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[15]+(int64_t)X[1]*X[14]+(int64_t)X[2]*X[13]+(int64_t)X[3]*X[12]+(int64_t)X[4]*X[11]+(int64_t)X[5]*X[10]+(int64_t)X[6]*X[9]+(int64_t)X[7]*X[8]) +4*((int64_t)X[16]*X[17]) + (C>>29);
    Z_res[15] = ((int32_t)C)&lower29bits;
    C = 2*((int64_t)X[0]*X[16]+(int64_t)X[1]*X[15]+(int64_t)X[2]*X[14]+(int64_t)X[3]*X[13]+(int64_t)X[4]*X[12]+(int64_t)X[5]*X[11]+(int64_t)X[6]*X[10]+(int64_t)X[7]*X[9]+(int64_t)X[17]*X[17]) +(int64_t)X[8]*X[8] + (C>>29);
    Z_res[16] = ((int32_t)C)&lower29bits;
    CC+= C>>29;
    Z_res[17]= ((int32_t)CC)&lower28bits;
    Z_res[0]+= (int32_t)(CC>>28);

}


//Modular Multiplication i.e. z=x*y mod p
/*The function computes TMVP of size 6 by using the 2-way decomposition formula in our paper.
mat = is one-dimensional array of 11 elements because it's Toeplitz
vec = is one-dimensional array of size 6
res = is one-dimensional array of size 6 such that it represents the resultant product/vector */
inline void TMVP_level2(const int32_t mat[], const int32_t vec[], int64_t res[]){
    int32_t temp_vec[3];
    int32_t temp_mat[4];

    temp_vec[0]= vec[0]-vec[3];
    temp_vec[1]= vec[1]-vec[4];
    temp_vec[2]= vec[2]-vec[5];

    res[3]= res[0]= (int64_t)mat[0]*temp_vec[0]+(int64_t)mat[1]*temp_vec[1]+(int64_t)mat[2]*temp_vec[2];
    res[4]= res[1]= (int64_t)mat[6]*temp_vec[0]+(int64_t)mat[0]*temp_vec[1]+(int64_t)mat[1]*temp_vec[2];
    res[5]= res[2]= (int64_t)mat[7]*temp_vec[0]+(int64_t)mat[6]*temp_vec[1]+(int64_t)mat[0]*temp_vec[2];

    temp_mat[0]= mat[3]+mat[0];
    temp_mat[1]= mat[4]+mat[1];
    temp_mat[2]= mat[2]+mat[6];
    temp_mat[3]= mat[1]+mat[7];

    res[0]+= (int64_t)temp_mat[0]*vec[3]+(int64_t)temp_mat[1]*vec[4]+(int64_t)(mat[5]+mat[2])*vec[5];
    res[1]+= (int64_t)temp_mat[2]*vec[3]+(int64_t)temp_mat[0]*vec[4]+(int64_t)temp_mat[1]*vec[5];
    res[2]+= (int64_t)temp_mat[3]*vec[3]+(int64_t)temp_mat[2]*vec[4]+(int64_t)temp_mat[0]*vec[5];

    temp_mat[0]= mat[8]+mat[0];
    temp_mat[1]= temp_mat[3];
    temp_mat[3]= mat[9]+mat[6];

    res[3]= (int64_t)temp_mat[0]*vec[0]+(int64_t)temp_mat[1]*vec[1]+(int64_t)temp_mat[2]*vec[2] - res[3];
    res[4]= (int64_t)temp_mat[3]*vec[0]+(int64_t)temp_mat[0]*vec[1]+(int64_t)temp_mat[1]*vec[2] - res[4];
    res[5]= (int64_t)(mat[10]+mat[7])*vec[0]+(int64_t)temp_mat[3]*vec[1]+(int64_t)temp_mat[0]*vec[2] - res[5];

}

/*The function computes the Z = X*Y mod (2**521) - 1] using TMVP where X, Y, and Z are 18-limb.
X_mat and Y_vec are one-dimensional arrays of integers where the first element is at most 30-bit, the next 16 elements are 29-bit at most and the last element is 28-bit integer at most.
The output Z_res is one-dimensional array of integers such that Z_res[0] in [0,2**30-1], Z_res[i] in [0,2**29-1] and Z_res[17] in [0,2**28-1]*/
void TMV_product(const int32_t X_mat[], const int32_t Y_vec[], int32_t Z_res[]){
    int32_t A3_2[]={2*X_mat[17], 2*X_mat[16], 2*X_mat[15], 2*X_mat[14], 2*X_mat[13]};           //A bit sloppy but it's okay because of the range of X
    int32_t temp_vec[6], temp_mat[11];                  //For second level call and taken signed to make compiler generate efficient code
    uint32_t temp1_mat[11];              //Very important to be taken unsigned 32-bit otherwise wrong result
    int64_t P2[6], P3[6], P4[6];

    //Computing P2
    temp_vec[0] = Y_vec[0] - Y_vec[6];          //No explicit type casting required because the result is within the range of signed integer
    temp_vec[1] = Y_vec[1] - Y_vec[7];          temp_vec[2] = Y_vec[2] - Y_vec[8];              temp_vec[3] = Y_vec[3] - Y_vec[9];
    temp_vec[4] = Y_vec[4] - Y_vec[10];         temp_vec[5] = Y_vec[5] - Y_vec[11];

    temp_mat[0] =  X_mat[6];        //No explicit type casting required because X is within the range of signed integer
    temp_mat[1] =  X_mat[5];                    temp_mat[2] =  X_mat[4];                        temp_mat[3] =  X_mat[3];
    temp_mat[4] =  X_mat[2];                    temp_mat[5] =  X_mat[1];                        temp_mat[6] =  X_mat[7];
    temp_mat[7] =  X_mat[8];                    temp_mat[8] =  X_mat[9];                        temp_mat[9] =  X_mat[10];
    temp_mat[10] = X_mat[11];

    TMVP_level2(temp_mat, temp_vec, P2);

    //Computing P4
    temp_vec[0] = Y_vec[6] - Y_vec[12];             //No explicit type casting required because the result is within the range of signed integer
    temp_vec[1] = Y_vec[7] - Y_vec[13];         temp_vec[2] = Y_vec[8] - Y_vec[14];             temp_vec[3] = Y_vec[9] - Y_vec[15];
    temp_vec[4] = Y_vec[10] - Y_vec[16];        temp_vec[5] = Y_vec[11] - Y_vec[17];

    temp_mat[0] =  2*X_mat[12];                 //No explicit type casting required because 2*X is within the range of signed integer
    temp_mat[1] =  2*X_mat[11];                 temp_mat[2] =  2*X_mat[10];                     temp_mat[3] =  2*X_mat[9];
    temp_mat[4] =  2*X_mat[8];                  temp_mat[5] =  2*X_mat[7];                      temp_mat[6] =  A3_2[4];
    temp_mat[7] =  A3_2[3];                     temp_mat[8] =  A3_2[2];                         temp_mat[9] =  A3_2[1];
    temp_mat[10] = A3_2[0];

    TMVP_level2(temp_mat, temp_vec, P4);

    //Computing P3
    temp_vec[0] = Y_vec[0] - Y_vec[12];             //No explicit type casting required because the result is within the range of signed integer
    temp_vec[1] = Y_vec[1] - Y_vec[13];         temp_vec[2] = Y_vec[2] - Y_vec[14];             temp_vec[3] = Y_vec[3] - Y_vec[15];
    temp_vec[4] = Y_vec[4] - Y_vec[16];         temp_vec[5] = Y_vec[5] - Y_vec[17];

    temp_mat[0] =  X_mat[0];            //No explicit type casting required because 2*X is within the range of signed integer
    temp_mat[1] =  A3_2[0];                     temp_mat[2] =  A3_2[1];                         temp_mat[3] =  A3_2[2];
    temp_mat[4] =  A3_2[3];                     temp_mat[5] =  A3_2[4];                         temp_mat[6] =  X_mat[1];
    temp_mat[7] =  X_mat[2];                    temp_mat[8] =  X_mat[3];                        temp_mat[9] =  X_mat[4];
    temp_mat[10] = X_mat[5];

    TMVP_level2(temp_mat, temp_vec, P3);


    //Computing P1
    temp1_mat[0]=X_mat[6] + X_mat[12];      //The result is 30-bit unsigned integer at the worst-case
    temp1_mat[1]=X_mat[5] + X_mat[11];          temp1_mat[2]=X_mat[4] + X_mat[10];              temp1_mat[3]=X_mat[3] + X_mat[9];
    temp1_mat[4]=X_mat[2] + X_mat[8];           temp1_mat[5]=X_mat[1] + X_mat[7];               temp1_mat[6]=X_mat[7] + X_mat[13];
    temp1_mat[7]=X_mat[8] + X_mat[14];          temp1_mat[8]=X_mat[9] + X_mat[15];              temp1_mat[9]=X_mat[10] + X_mat[16];
    temp1_mat[10]=X_mat[11] + X_mat[17];

    temp_mat[0]=temp1_mat[0] + X_mat[0];       //The result is 31-bit unsigned integer at the worst-case
    temp_mat[1]=temp1_mat[1] + A3_2[0];        temp_mat[2]=temp1_mat[2] + A3_2[1];            temp_mat[3]=temp1_mat[3] + A3_2[2];
    temp_mat[4]=temp1_mat[4] + A3_2[3];        temp_mat[5]=temp1_mat[5] + A3_2[4];            temp_mat[6]=temp1_mat[6] + X_mat[1];
    temp_mat[7]=temp1_mat[7] + X_mat[2];       temp_mat[8]=temp1_mat[8] + X_mat[3];           temp_mat[9]=temp1_mat[9] + X_mat[4];
    temp_mat[10]=temp1_mat[10] + X_mat[5];


    //Computing P6, remember that matrix size is 32-bit unsigned and same is true for P5
    temp1_mat[0]+=temp_mat[0];                 temp1_mat[1]+=temp_mat[1];                     temp1_mat[2]+=temp_mat[2];
    temp1_mat[3]+=temp_mat[3];                 temp1_mat[4]+=temp_mat[4];                     temp1_mat[5]+=temp_mat[5];
    temp1_mat[6]+=temp_mat[6];                 temp1_mat[7]+=temp_mat[7];                     temp1_mat[8]+=temp_mat[8];
    temp1_mat[9]+=temp_mat[9];                 temp1_mat[10]+=temp_mat[10];


    //Carry Propagation
    uint64_t temp= (int64_t)temp_mat[10]*Y_vec[0] + (int64_t)temp_mat[9]*Y_vec[1] + (int64_t)temp_mat[8]*Y_vec[2] + (int64_t)temp_mat[7]*Y_vec[3] + (int64_t)temp_mat[6]*Y_vec[4]  + (int64_t)temp_mat[0]*Y_vec[5] - P2[5]-P3[5];
    uint64_t C = temp&lower28bits;            //From testing we have found that it should be double-word especially at the end
    temp= (P3[0]+P4[0]+ (uint64_t)temp1_mat[0]*(uint32_t)Y_vec[12] + (uint64_t)temp1_mat[1]*(uint32_t)Y_vec[13]  + (uint64_t)temp1_mat[2]*(uint32_t)Y_vec[14] + (uint64_t)temp1_mat[3]*(uint32_t)Y_vec[15] + (uint64_t)temp1_mat[4]*(uint32_t)Y_vec[16]  + (uint64_t)temp1_mat[5]*(uint32_t)Y_vec[17]) + (temp>>28);        //Testing shows this 64-bit
    Z_res[0]= ((int32_t)temp)&lower29bits;
    temp= (P3[1]+P4[1]+ (uint64_t)temp1_mat[6]*(uint32_t)Y_vec[12] + (uint64_t)temp1_mat[0]*(uint32_t)Y_vec[13]  + (uint64_t)temp1_mat[1]*(uint32_t)Y_vec[14] + (uint64_t)temp1_mat[2]*(uint32_t)Y_vec[15] + (uint64_t)temp1_mat[3]*(uint32_t)Y_vec[16]  + (uint64_t)temp1_mat[4]*(uint32_t)Y_vec[17]) + (uint64_t)(temp>>29);    //To ensure the carry is +ve
    Z_res[1]= ((int32_t)temp)&lower29bits;
    temp= (P3[2]+P4[2]+ (uint64_t)temp1_mat[7]*(uint32_t)Y_vec[12] + (uint64_t)temp1_mat[6]*(uint32_t)Y_vec[13]  + (uint64_t)temp1_mat[0]*(uint32_t)Y_vec[14] + (uint64_t)temp1_mat[1]*(uint32_t)Y_vec[15] + (uint64_t)temp1_mat[2]*(uint32_t)Y_vec[16]  + (uint64_t)temp1_mat[3]*(uint32_t)Y_vec[17]) + (temp>>29);
    Z_res[2]= ((int32_t)temp)&lower29bits;
    temp= (P3[3]+P4[3]+ (uint64_t)temp1_mat[8]*(uint32_t)Y_vec[12] + (uint64_t)temp1_mat[7]*(uint32_t)Y_vec[13]  + (uint64_t)temp1_mat[6]*(uint32_t)Y_vec[14] + (uint64_t)temp1_mat[0]*(uint32_t)Y_vec[15] + (uint64_t)temp1_mat[1]*(uint32_t)Y_vec[16]  + (uint64_t)temp1_mat[2]*(uint32_t)Y_vec[17]) + (temp>>29);
    Z_res[3]= ((int32_t)temp)&lower29bits;
    temp= (P3[4]+P4[4]+ (uint64_t)temp1_mat[9]*(uint32_t)Y_vec[12] + (uint64_t)temp1_mat[8]*(uint32_t)Y_vec[13]  + (uint64_t)temp1_mat[7]*(uint32_t)Y_vec[14] + (uint64_t)temp1_mat[6]*(uint32_t)Y_vec[15] + (uint64_t)temp1_mat[0]*(uint32_t)Y_vec[16]  + (uint64_t)temp1_mat[1]*(uint32_t)Y_vec[17]) + (temp>>29);
    Z_res[4]= ((int32_t)temp)&lower29bits;
    temp= (P3[5]+P4[5]+ (uint64_t)temp1_mat[10]*(uint32_t)Y_vec[12] + (uint64_t)temp1_mat[9]*(uint32_t)Y_vec[13] + (uint64_t)temp1_mat[8]*(uint32_t)Y_vec[14] + (uint64_t)temp1_mat[7]*(uint32_t)Y_vec[15] + (uint64_t)temp1_mat[6]*(uint32_t)Y_vec[16]  + (uint64_t)temp1_mat[0]*(uint32_t)Y_vec[17]) + (temp>>29);
    Z_res[5]= ((int32_t)temp)&lower29bits;
    //For P5
    temp1_mat[0]=temp_mat[0]+X_mat[12];        //At most 32-bit unsigned

    temp= (P2[0]-P4[0]+ (uint64_t)temp1_mat[0]*(uint32_t)Y_vec[6]  + (uint64_t)temp1_mat[10]*(uint32_t)Y_vec[7] + (uint64_t)temp1_mat[9]*(uint32_t)Y_vec[8]  + (uint64_t)temp1_mat[8]*(uint32_t)Y_vec[9]  + (uint64_t)temp1_mat[7]*(uint32_t)Y_vec[10]  + (uint64_t)temp1_mat[6]*(uint32_t)Y_vec[11])+ (temp>>29);
    Z_res[6]= ((int32_t)temp)&lower29bits;
    temp= (P2[1]-P4[1]+ (int64_t)temp_mat[5]*Y_vec[6]  + (uint64_t)temp1_mat[0]*(uint32_t)Y_vec[7]  + (uint64_t)temp1_mat[10]*(uint32_t)Y_vec[8] + (uint64_t)temp1_mat[9]*(uint32_t)Y_vec[9]  + (uint64_t)temp1_mat[8]*(uint32_t)Y_vec[10]  + (uint64_t)temp1_mat[7]*(uint32_t)Y_vec[11]) + (temp>>29);
    Z_res[7]= ((int32_t)temp)&lower29bits;
    temp= (P2[2]-P4[2]+ (int64_t)temp_mat[4]*Y_vec[6]  + (int64_t)temp_mat[5]*Y_vec[7]  + (uint64_t)temp1_mat[0]*(uint32_t)Y_vec[8]  + (uint64_t)temp1_mat[10]*(uint32_t)Y_vec[9] + (uint64_t)temp1_mat[9]*(uint32_t)Y_vec[10]  + (uint64_t)temp1_mat[8]*(uint32_t)Y_vec[11])+ (temp>>29);
    Z_res[8]= ((int32_t)temp)&lower29bits;
    temp= (P2[3]-P4[3]+ (int64_t)temp_mat[3]*Y_vec[6]  + (int64_t)temp_mat[4]*Y_vec[7]  + (int64_t)temp_mat[5]*Y_vec[8]  + (uint64_t)temp1_mat[0]*(uint32_t)Y_vec[9]  + (uint64_t)temp1_mat[10]*(uint32_t)Y_vec[10] + (uint64_t)temp1_mat[9]*(uint32_t)Y_vec[11]) + (temp>>29);
    Z_res[9]= ((int32_t)temp)&lower29bits;
    temp= (P2[4]-P4[4]+ (int64_t)temp_mat[2]*Y_vec[6]  + (int64_t)temp_mat[3]*Y_vec[7]  + (int64_t)temp_mat[4]*Y_vec[8]  + (int64_t)temp_mat[5]*Y_vec[9]  + (uint64_t)temp1_mat[0]*(uint32_t)Y_vec[10]  + (uint64_t)temp1_mat[10]*(uint32_t)Y_vec[11]) + (temp>>29);
    Z_res[10]= ((int32_t)temp)&lower29bits;
    temp= (P2[5]-P4[5]+ (int64_t)temp_mat[1]*Y_vec[6]  + (int64_t)temp_mat[2]*Y_vec[7]  + (int64_t)temp_mat[3]*Y_vec[8]  + (int64_t)temp_mat[4]*Y_vec[9]  + (int64_t)temp_mat[5]*Y_vec[10]  + (int64_t)temp1_mat[0]*(uint32_t)Y_vec[11]) + (temp>>29);
    Z_res[11]= ((int32_t)temp)&lower29bits;
    temp= ((int64_t)temp_mat[0]*Y_vec[0] + (int64_t)temp_mat[1]*Y_vec[1]  + (int64_t)temp_mat[2]*Y_vec[2] + (int64_t)temp_mat[3]*Y_vec[3] + (int64_t)temp_mat[4]*Y_vec[4]  + (int64_t)temp_mat[5]*Y_vec[5] - P2[0]-P3[0]) + (temp>>29);
    Z_res[12]= ((int32_t)temp)&lower29bits;
    temp= ((int64_t)temp_mat[6]*Y_vec[0] + (int64_t)temp_mat[0]*Y_vec[1]  + (int64_t)temp_mat[1]*Y_vec[2] + (int64_t)temp_mat[2]*Y_vec[3] + (int64_t)temp_mat[3]*Y_vec[4]  + (int64_t)temp_mat[4]*Y_vec[5] - P2[1]-P3[1]) + (temp>>29);
    Z_res[13]= ((int32_t)temp)&lower29bits;
    temp= ((int64_t)temp_mat[7]*Y_vec[0] + (int64_t)temp_mat[6]*Y_vec[1]  + (int64_t)temp_mat[0]*Y_vec[2] + (int64_t)temp_mat[1]*Y_vec[3] + (int64_t)temp_mat[2]*Y_vec[4]  + (int64_t)temp_mat[3]*Y_vec[5] - P2[2]-P3[2]) + (temp>>29);
    Z_res[14]= ((int32_t)temp)&lower29bits;
    temp= ((int64_t)temp_mat[8]*Y_vec[0] + (int64_t)temp_mat[7]*Y_vec[1]  + (int64_t)temp_mat[6]*Y_vec[2] + (int64_t)temp_mat[0]*Y_vec[3] + (int64_t)temp_mat[1]*Y_vec[4]  + (int64_t)temp_mat[2]*Y_vec[5] - P2[3]-P3[3]) + (temp>>29);
    Z_res[15]= ((int32_t)temp)&lower29bits;
    temp= ((int64_t)temp_mat[9]*Y_vec[0] + (int64_t)temp_mat[8]*Y_vec[1]  + (int64_t)temp_mat[7]*Y_vec[2] + (int64_t)temp_mat[6]*Y_vec[3] + (int64_t)temp_mat[0]*Y_vec[4]  + (int64_t)temp_mat[1]*Y_vec[5] -P2[4]-P3[4]) + (temp>>29);
    Z_res[16]= ((int32_t)temp)&lower29bits;
    C+= (temp >> 29);
    Z_res[17]= ((int32_t)C)&lower28bits;
    Z_res[0]+= (int32_t)(C>>28);

}


// Inverse x = 1/x = x^(p-2) mod p
// 13 muls, 520 sqrs
void ginv(int32_t x[]){
	int32_t x127[LIMBS], w[LIMBS], t[LIMBS], z[LIMBS];
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

///////////////// Curve Operations ////////////////////////////////////////////////
// Point Structure

typedef struct {
    int32_t x[LIMBS];
    int32_t y[LIMBS];
    int32_t z[LIMBS];
    int32_t t[LIMBS];
} ECp;

//P=0 Neutral element
void inf(ECp *P){
	for (int32_t i=0;i<LIMBS;i++)
		P->x[i]=P->y[i]=P->z[i]=P->t[i]=0;
	P->y[0]=P->z[0]=1;
}

// Initialize P
void init(int32_t *x, int32_t *y, ECp *P){
	for (int32_t i=0;i<LIMBS;i++){
		P->x[i]=x[i];
		P->y[i]=y[i];
		P->z[i]=0;
	}
	P->z[0]=1;
	TMV_product(x,y,P->t);
}

//P=Q
void copy(ECp *Q,ECp *P){
	for (int32_t i=0;i<LIMBS;i++){
		P->x[i]=Q->x[i];
		P->y[i]=Q->y[i];
		P->z[i]=Q->z[i];
		P->t[i]=Q->t[i];
	}
}

// P=-Q
void neg(ECp *Q,ECp *P){
    int32_t xt[LIMBS], tt[LIMBS];
	for (int32_t i=0;i<LIMBS;i++){
		xt[i]=-Q->x[i];             //Could be at most 31-bit signed
		P->y[i]=Q->y[i];
		P->z[i]=Q->z[i];
		tt[i]=-Q->t[i];             //Could be at most 31-bit signed
	}
	scr(xt, P->x);          //To avoid overflow
	scr(tt, P->t);          //To avoid overflow
}

//From Extended projective to Affine coordinate
void norm(ECp *P){
	int32_t w[LIMBS], t[LIMBS];
	gcopy(P->z,w);
	ginv(w);
	TMV_product(P->x, w, t);            //t is in intermediate form
	scr(t, P->x);                       //P->x is a unique residue in modulo p
	TMV_product(P->y,w,t);              //t is in intermediate form
	scr(t, P->y);                       //P->y is a unique residue in modulo p
	/*TMV_product(P->z,w,t);            //t is in intermediate form
	scr(t, P->z);                       //P->z is a unique residue in modulo p
	TMV_product(P->t,w,t);              //t is in intermediate form
	scr(t, P->t);                       //P->t is a unique residue in modulo p*/
}


// P+=P where P->z=1
void double_1(ECp *P){
	int32_t a[LIMBS],b[LIMBS],e[LIMBS],f[LIMBS],g[LIMBS], h[LIMBS];
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
	gsub2(P->z,g,P->z);
	TMV_product(e,h,P->t);
}

// P+=P where the extended coordinate t is not required
void double_2(ECp *P){
	int32_t a[LIMBS],b[LIMBS],c[LIMBS],e[LIMBS],f[LIMBS],g[LIMBS],h[LIMBS];
	gsqr(P->x,a);
	gsqr(P->y,b);
	gsqr(P->z,c);
	gadd(P->x,P->y,g); gsqr(g,e); gtsb(a, b, e);
	gadd(a,b,g);
	gsub2(g,c,f);
	gsub(a,b,h);
	TMV_product(e,f,P->x);
	TMV_product(g,h,P->y);
	TMV_product(f,g,P->z);
}

// P+=P where the extended coordinate t is required
void double_3(ECp *P){
	int32_t a[LIMBS],b[LIMBS],c[LIMBS],e[LIMBS],f[LIMBS],g[LIMBS],h[LIMBS];
	gsqr(P->x,a);
	gsqr(P->y,b);
	gsqr(P->z,c);
	gadd(P->x,P->y,g); gsqr(g,e); gtsb(a, b, e);
	gadd(a,b,g);
	gsub2(g,c,f);
	gsub(a,b,h);
	TMV_product(e,f,P->x);
	TMV_product(g,h,P->y);
	TMV_product(f,g,P->z);
	TMV_product(e,h,P->t);
}

//P+=Q, when P->t is required and Q->z is 1

void add_1(ECp *Q,ECp *P){
	int32_t a[LIMBS],b[LIMBS],c[LIMBS],d[LIMBS],e[LIMBS],f[LIMBS],g[LIMBS],h[LIMBS];
	TMV_product(P->x,Q->x,a);
	TMV_product(P->y,Q->y,b);
	TMV_product(P->t,Q->t,c);
	gadd(P->z,c,f);     //reversed sign as d is negative
	gsub(P->z,c,g);     //reversed sign as d is negative
	gsub(b,a,h);
	gadd(P->x,P->y,c); gadd(Q->x,Q->y,d); TMV_product(c,d,e); gtsb(a, b, e);
	TMV_product(e,f,P->x);
	TMV_product(g,h,P->y);
	TMV_product(f,g,P->z);
	TMV_product(e,h,P->t);
}

//P+=Q where the extended coordinate t is not required
void add_2(ECp *Q,ECp *P){
	int32_t a[LIMBS],b[LIMBS],c[LIMBS],d[LIMBS],e[LIMBS],f[LIMBS],g[LIMBS],h[LIMBS];
	TMV_product(P->x,Q->x,a);
	TMV_product(P->y,Q->y,b);
	TMV_product(P->t,Q->t,c);
	TMV_product(P->z,Q->z,d);
	gadd(d,c,f);        //reversed sign as d is negative
	gsub(d,c,g);        //reversed sign as d is negative
	gsub(b,a,h);
	gadd(P->x,P->y,c); gadd(Q->x,Q->y,d); TMV_product(c,d,e); gtsb(a, b, e);
	TMV_product(e,f,P->x);
	TMV_product(g,h,P->y);
	TMV_product(f,g,P->z);
}


// Pre-computation
void precomp(ECp *P,ECp W[]){
	inf(&W[0]);                                             //W[0]= Neutral element/point
	copy(P,&W[1]);      gmuli(W[1].t,CURVE_PARAMETER);      //W[1]=P
	copy(P,&W[2]);      double_1(&W[2]);            //W[2]=2P
	copy(&W[2],&W[3]);  add_1(&W[1],&W[3]);         //W[3]=3P
	copy(&W[2],&W[4]);  double_3(&W[4]);            //W[4]=4P
	copy(&W[4],&W[5]);  add_1(&W[1],&W[5]);         //W[5]=5P
	copy(&W[3],&W[6]);  double_3(&W[6]);            //W[6]=6P
	copy(&W[6],&W[7]);  add_1(&W[1],&W[7]);         //W[7]=7P
	copy(&W[4],&W[8]);  double_3(&W[8]);            //W[8]=8P

#if WINDOW>4
	copy(&W[8],&W[9]); add_1(&W[1],&W[9]);          //W[9]=9P
	copy(&W[5],&W[10]); double_3(&W[10]);           //W[10]=10P
	copy(&W[10],&W[11]); add_1(&W[1],&W[11]);       //W[11]=11P
	copy(&W[6],&W[12]); double_3(&W[12]);           //W[12]=12P
	copy(&W[12],&W[13]); add_1(&W[1],&W[13]);       //W[13]=13P
	copy(&W[7],&W[14]); double_3(&W[14]);           //W[14]=14P
	copy(&W[14],&W[15]); add_1(&W[1],&W[15]);       //W[15]=15P
	copy(&W[8],&W[16]); double_3(&W[16]);           //W[16]=16P
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

/* premultiply t parameter by curve constant */

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
        printf("x[%d]= %x\n", i, P->x[i]);

	puts("");
	for(i=0; i<LIMBS; i++)
        printf("y[%d]= %x\n", i, P->y[i]);

    puts("");
}

/*Constant time table look-up - borrowed from ed25519 */

void fe_cmov(int32_t f[], int32_t g[], int32_t ib){
  int32_t b=-ib;
  f[0]^=(f[0]^g[0])&b;
  f[1]^=(f[1]^g[1])&b;
  f[2]^=(f[2]^g[2])&b;
  f[3]^=(f[3]^g[3])&b;
  f[4]^=(f[4]^g[4])&b;
  f[5]^=(f[5]^g[5])&b;
  f[6]^=(f[6]^g[6])&b;
  f[7]^=(f[7]^g[7])&b;
  f[8]^=(f[8]^g[8])&b;
  f[9]^=(f[9]^g[9])&b;
  f[10]^=(f[10]^g[10])&b;
  f[11]^=(f[11]^g[11])&b;
  f[12]^=(f[12]^g[12])&b;
  f[13]^=(f[13]^g[13])&b;
  f[14]^=(f[14]^g[14])&b;
  f[15]^=(f[15]^g[15])&b;
  f[16]^=(f[16]^g[16])&b;
  f[17]^=(f[17]^g[17])&b;
}

static void cmov(ECp *w,ECp *u,int32_t b){
  fe_cmov(w->x,u->x,b);
  fe_cmov(w->y,u->y,b);
  fe_cmov(w->z,u->z,b);
  fe_cmov(w->t,u->t,b);
}

// return 1 if b==c, no branching
static int32_t equal(int32_t b,int32_t c){
	int32_t x=b^c;
	x-=1;  // if x=0, x now -1
	return ((x>>31)&1);
}

static void ct_select(ECp *T,ECp W[],int32_t b){
  ECp MT;
  int32_t m=b>>31;
  int32_t babs=(b^m)-m;

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


/* Point Multiplication - exponent is 519 bits */
void scalar_mul(int32_t *w,ECp *P){
	ECp W[M+1], Q;

	precomp(P,W);
	copy(&W[w[PANES-1]], P);
	for (int32_t i=PANES-2; i>=0; i--){
		ct_select(&Q, W, w[i]);
		window(&Q,P);
    }
	norm(P);
}


/////////////////////////////////////////////////////////////////////////////
int main(){
	uint64_t bef, aft, min_ccycle=100000000;
	int32_t w[PANES], i, j, lpz=10000;
	int32_t xs[LIMBS],ys[LIMBS];
	ECp P;

/* Base point on Edwards Curve (from SafeCurves Website) */

    xs[0]=0xf19ba6c;        xs[1]=0x154a051;        xs[2]=0x120e2a8c;       xs[3]=0x1f6266c;        xs[4]=0x19c6059d;       xs[5]=0xeab47e4;
	xs[6]=0x12c6ba52;       xs[7]=0x1998e486;       xs[8]=0x13f6ecc5;       xs[9]=0x63101c8;        xs[10]=0x2fcf270;       xs[11]=0xd9031d9;
	xs[12]=0x1d9f42fc;      xs[13]=0x143c51df;      xs[14]=0x12c8a5ac;      xs[15]=0x313bf21;       xs[16]=0x1c48648b;      xs[17]=0x3a965a2;

	ys[0]=0xc;              ys[1]=0x0;              ys[2]=0x0;              ys[3]=0x0;              ys[4]=0x0;              ys[5]=0x0;
	ys[6]=0x0;              ys[7]=0x0;              ys[8]=0x0;              ys[9]=0x0;              ys[10]=0x0;             ys[11]=0x0;
	ys[12]=0x0;             ys[13]=0x0;             ys[14]=0x0;             ys[15]=0x0;             ys[16]=0x0;             ys[17]=0x0;


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
	printf("Window width :: %u \nThe minimum clock cycles count is :: %"PRIu64"\n", WINDOW, min_ccycle);
	output(&P);

}
