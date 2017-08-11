/*The following information are taken from http://safecurves.cr.yp.to/base.html
Curve Equation :: x^2+y^2 = 1-376014x^2y^2 where d=-376014
bitlength(ECsGord)= 519
Test program for ed521 scalar point multiplication
Uses Daniel J. Bernstein et al.s point multiplication method from Curve41417 paper
Cache safety thanks to ed25519
Fully Tested and debugged
We have replaced gmul() by TMV_product()
We have added rdtscp() to measure the clock cycles count as suggested by Paoloni; http://www.intel.com.tr/content/dam/www/public/us/en/documents/white-papers/ia-32-ia-64-benchmark-code-execution-paper.pdf
We have amended scr(), gmuli(), gsqr(), and gsqr2() according to modulus p
gcc -Wall -m32 -O3 fb32_ed521.c -o fb32_ed521.exe */


#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>
#include<math.h>


#define WINDOW 4 //6 //5

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
#define LIMBS 18
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
   __asm__ __volatile__ ("rdtscp": "=a" (lo), "=d" (hi) :: "%rcx");
   __asm__ __volatile__ ("cpuid"::: "%rax", "%rbx", "%rcx", "%rdx");
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
void TMV_product(int32_t X_mat[], int32_t Y_vec[], int32_t Z_res[]){
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


///////////////// Elliptic Curve sub-Group Operations //////////////////////////////

// Point Structure
typedef struct{
    int32_t x[LIMBS];
    int32_t y[LIMBS];
    int32_t z[LIMBS];
    int32_t t[LIMBS];
} ECp;


//P=0, Neutral element
void inf(ECp *P){
	for (int i=0;i<LIMBS;i++)
		P->x[i]=P->y[i]=P->z[i]=P->t[i]=0;
	P->y[0]=P->z[0]=1;
}

// Initialize P
void init(int32_t x[], int32_t y[], ECp *P){
	for (int32_t i=0;i<LIMBS;i++){
		P->x[i]=x[i];
		P->y[i]=y[i];
		P->z[i]=0;
	}
	P->z[0]=1;
	TMV_product(x,y,P->t);
}

//P=Q
void copy(ECp *Q, ECp *P){
	for (int32_t i=0;i<LIMBS;i++){
		P->x[i]=Q->x[i];
		P->y[i]=Q->y[i];
		P->z[i]=Q->z[i];
		P->t[i]=Q->t[i];
	}
}

// P=-Q
void neg(ECp *Q, ECp *P){
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


// P+=P where the extended coordinate t is not required
void doubl_2(ECp *P){
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
void doubl_3(ECp *P){
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

//P+=Q where the extended coordinate t is not required
void add_2(ECp *Q, ECp *P){
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

//P+=Q where the extended coordinate t is required
void add_3(ECp *Q, ECp *P){
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
	TMV_product(e,h,P->t);
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


// Printing the Affine coordinates of the computed new point on the screen
void output(ECp *P){
    int32_t i;
    for(i=0; i<LIMBS; i++)
        printf("x[%d]= %X\n", i, P->x[i]);

	puts("");
	for(i=0; i<LIMBS; i++)
        printf("y[%d]= %X\n", i, P->y[i]);

    puts("");

}


/*Constant time table look-up - borrowed from ed25519*/
void fe_cmov(int32_t f[], int32_t g[], int32_t ib){
  f[0]^=(f[0]^g[0])&ib;
  f[1]^=(f[1]^g[1])&ib;
  f[2]^=(f[2]^g[2])&ib;
  f[3]^=(f[3]^g[3])&ib;
  f[4]^=(f[4]^g[4])&ib;
  f[5]^=(f[5]^g[5])&ib;
  f[6]^=(f[6]^g[6])&ib;
  f[7]^=(f[7]^g[7])&ib;
  f[8]^=(f[8]^g[8])&ib;
  f[9]^=(f[9]^g[9])&ib;
  f[10]^=(f[10]^g[10])&ib;
  f[11]^=(f[11]^g[11])&ib;
  f[12]^=(f[12]^g[12])&ib;
  f[13]^=(f[13]^g[13])&ib;
  f[14]^=(f[14]^g[14])&ib;
  f[15]^=(f[15]^g[15])&ib;
  f[16]^=(f[16]^g[16])&ib;
  f[17]^=(f[17]^g[17])&ib;
}

static void cmov(ECp *w, ECp *u, int32_t b){
  fe_cmov(w->x,u->x,b);
  fe_cmov(w->y,u->y,b);
  fe_cmov(w->z,u->z,b);
  fe_cmov(w->t,u->t,b);
}

// Due the sign table we have changed the equality function as follows
// return -1 if b==c, else 0, no branching
static int32_t equal(const uint32_t b, const uint32_t c){
	int32_t x=b^c;
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

    //P+=Q where P is neutral element
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

///////////////// Main Function //////////////////////////////
int main(void){
	uint64_t bef, aft, min_ccycle=10000000000;
	int32_t i, j, lpz=10000;
	ECp P, PreT[M][V];
	int32_t xs[LIMBS], ys[LIMBS];

/* Base point on Edwards Curve (from SafeCurves Website) */

    xs[0]=0xf19ba6c;        xs[1]=0x154a051;        xs[2]=0x120e2a8c;       xs[3]=0x1f6266c;        xs[4]=0x19c6059d;       xs[5]=0xeab47e4;
	xs[6]=0x12c6ba52;       xs[7]=0x1998e486;       xs[8]=0x13f6ecc5;       xs[9]=0x63101c8;        xs[10]=0x2fcf270;       xs[11]=0xd9031d9;
	xs[12]=0x1d9f42fc;      xs[13]=0x143c51df;      xs[14]=0x12c8a5ac;  	xs[15]=0x313bf21;       xs[16]=0x1c48648b;      xs[17]=0x3a965a2;

	ys[0]=0xc;              ys[1]=0x0;              ys[2]=0x0;          	ys[3]=0x0;              ys[4]=0x0;              ys[5]=0x0;
	ys[6]=0x0;              ys[7]=0x0;              ys[8]=0x0;          	ys[9]=0x0;              ys[10]=0x0;             ys[11]=0x0;
	ys[12]=0x0;             ys[13]=0x0;             ys[14]=0x0;             ys[15]=0x0;             ys[16]=0x0;             ys[17]=0x0;

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
	printf("Window width :: %u\n%s%"PRIu64"\n", WINDOW, "The minimum clock cycles count is :: ", min_ccycle);
	output(&P);

}
