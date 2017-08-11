/* Test program for ws521 scalar point multiplication
 This is NIST standard (short) Weierstrass curve p-521
 Fully Tested and debugged
 Uses constant time method described by Bos et al. - http://eprint.iacr.org/2014/130
 Cache safety thanks to ed25519
 M.Scott 27/02/2015
 We have replaced gmul() by TMV_product()
 We have added rdtscp() to measure the clock cycles count as suggested by Paoloni; http://www.intel.com.tr/content/dam/www/public/us/en/documents/white-papers/ia-32-ia-64-benchmark-code-execution-paper.pdf
 We have amended scr() and gsqr() according to modulus p
 gcc -Wall -m32 -O3 vb32_ws521.c -o vb32_ws521.exe */

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
#define LIMBS 18

//#define TEST  /* define to multiply by the group order */

static const int32_t lower29bits = 0x1fffffff;
static const int32_t lower28bits = 0xfffffff;

__inline__ uint64_t rdtsc(){
   uint32_t lo, hi;
   __asm__ __volatile__ ("xorl %%eax,%%eax \n        cpuid" ::: "%rax", "%rbx", "%rcx", "%rdx");
   __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));

   return (uint64_t)hi << 32 | lo;
}

__inline__ uint64_t rdtscp(){
   uint32_t lo, hi;
   __asm__ __volatile__ ("rdtscp" : "=a" (lo), "=d" (hi) :: "%rcx");
   __asm__ __volatile__ ("cpuid" ::: "%rax", "%rbx","%rcx", "%rdx");
   return (uint64_t)hi << 32 | lo;
}

///////////////////////// Field Operations ////////////////////////////////////
// w=1
void gone(int32_t w[]){
	w[0]=1;
	w[1]= w[2]= w[3]= w[4]= w[5]= w[6]= w[7]= w[8]= w[9]= 0;
	w[10]= w[11]= w[12]= w[13]= w[14]= w[15]= w[16]= w[17]=0;
}

// w=x+y and the result is in reduced form
void gadd(int32_t x[], int32_t y[], int32_t w[]){

	w[0]=x[0]+y[0];
	w[1]=x[1]+y[1];
	w[2]=x[2]+y[2];
	w[3]=x[3]+y[3];
	w[4]=x[4]+y[4];
	w[5]=x[5]+y[5];
	w[6]=x[6]+y[6];
	w[7]=x[7]+y[7];
	w[8]=x[8]+y[8];
	w[9]=x[9]+y[9];
	w[10]=x[10]+y[10];
	w[11]=x[11]+y[11];
	w[12]=x[12]+y[12];
	w[13]=x[13]+y[13];
	w[14]=x[14]+y[14];
	w[15]=x[15]+y[15];
	w[16]=x[16]+y[16];
	w[17]=x[17]+y[17];
}

// w=x-y and the result is in reduced form
void gsub(int32_t x[], int32_t y[], int32_t w[]){
	w[0]=x[0]-y[0];
	w[1]=x[1]-y[1];
	w[2]=x[2]-y[2];
	w[3]=x[3]-y[3];
	w[4]=x[4]-y[4];
	w[5]=x[5]-y[5];
	w[6]=x[6]-y[6];
	w[7]=x[7]-y[7];
	w[8]=x[8]-y[8];
	w[9]=x[9]-y[9];
	w[10]=x[10]-y[10];
	w[11]=x[11]-y[11];
	w[12]=x[12]-y[12];
	w[13]=x[13]-y[13];
	w[14]=x[14]-y[14];
	w[15]=x[15]-y[15];
	w[16]=x[16]-y[16];
	w[17]=x[17]-y[17];
}


//w=w-x-y and the result in scr
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

//w=w-x-2y and the result in reduced limb form
void gtsb2(int32_t x[], int32_t y[], int32_t w[]){
    int32_t t;
	t=w[0]-x[0]-2*y[0];                 //Could be at most 32-bit
	w[0]=t&lower29bits;

	t=(w[1]-x[1]-2*y[1])+(t>>29);              w[1]=t&lower29bits;
	t=(w[2]-x[2]-2*y[2])+(t>>29);              w[2]=t&lower29bits;
	t=(w[3]-x[3]-2*y[3])+(t>>29);              w[3]=t&lower29bits;
	t=(w[4]-x[4]-2*y[4])+(t>>29);              w[4]=t&lower29bits;
	t=(w[5]-x[5]-2*y[5])+(t>>29);              w[5]=t&lower29bits;
	t=(w[6]-x[6]-2*y[6])+(t>>29);              w[6]=t&lower29bits;
	t=(w[7]-x[7]-2*y[7])+(t>>29);              w[7]=t&lower29bits;
	t=(w[8]-x[8]-2*y[8])+(t>>29);              w[8]=t&lower29bits;
	t=(w[9]-x[9]-2*y[9])+(t>>29);              w[9]=t&lower29bits;
	t=(w[10]-x[10]-2*y[10])+(t>>29);            w[10]=t&lower29bits;
	t=(w[11]-x[11]-2*y[11])+(t>>29);            w[11]=t&lower29bits;
	t=(w[12]-x[12]-2*y[12])+(t>>29);            w[12]=t&lower29bits;
	t=(w[13]-x[13]-2*y[13])+(t>>29);            w[13]=t&lower29bits;
	t=(w[14]-x[14]-2*y[14])+(t>>29);            w[14]=t&lower29bits;
	t=(w[15]-x[15]-2*y[15])+(t>>29);            w[15]=t&lower29bits;
	t=(w[16]-x[16]-2*y[16])+(t>>29);            w[16]=t&lower29bits;
	t=(w[17]-x[17]-2*y[17])+(t>>29);            w[17]=t&lower28bits;               //Most significant limb is at most 28-bit

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

// w*=2 and the result is in reduced form
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


// w*=3 and the result is in reduced form
void gmul3(int32_t w[]){
	int32_t t;
	t=w[0]*3;		            //Could be at most 31-bit
	w[0]=t&lower29bits;

	t=(w[1]*3)+(t>>29);         w[1]=t&lower29bits;
	t=(w[2]*3)+(t>>29);         w[2]=t&lower29bits;
	t=(w[3]*3)+(t>>29);         w[3]=t&lower29bits;
	t=(w[4]*3)+(t>>29);       	w[4]=t&lower29bits;
	t=(w[5]*3)+(t>>29);       	w[5]=t&lower29bits;
	t=(w[6]*3)+(t>>29);       	w[6]=t&lower29bits;
	t=(w[7]*3)+(t>>29);       	w[7]=t&lower29bits;
	t=(w[8]*3)+(t>>29);       	w[8]=t&lower29bits;
    t=(w[9]*3)+(t>>29);       	w[9]=t&lower29bits;
	t=(w[10]*3)+(t>>29);        w[10]=t&lower29bits;
	t=(w[11]*3)+(t>>29);      	w[11]=t&lower29bits;
	t=(w[12]*3)+(t>>29);      	w[12]=t&lower29bits;
	t=(w[13]*3)+(t>>29);      	w[13]=t&lower29bits;
	t=(w[14]*3)+(t>>29);      	w[14]=t&lower29bits;
	t=(w[15]*3)+(t>>29);      	w[15]=t&lower29bits;
	t=(w[16]*3)+(t>>29);        w[16]=t&lower29bits;
	t=(w[17]*3)+(t>>29);        w[17]=t&lower28bits;           //most significant limb is at most 28-bit

	w[0]+=(t>>28);				//least significant limb could be at most 30-bit
}

// Before calling make sure that all the limbs are at most 29-bit
// w*=4 and the result is in reduced form
void gmul4(int32_t w[]){
	int32_t t;
	t=w[0]*4;
	w[0]=t&lower29bits;

	t=(w[1]*4)+(t>>29);         w[1]=t&lower29bits;
	t=(w[2]*4)+(t>>29);         w[2]=t&lower29bits;
	t=(w[3]*4)+(t>>29);         w[3]=t&lower29bits;
	t=(w[4]*4)+(t>>29);       	w[4]=t&lower29bits;
	t=(w[5]*4)+(t>>29);       	w[5]=t&lower29bits;
	t=(w[6]*4)+(t>>29);       	w[6]=t&lower29bits;
	t=(w[7]*4)+(t>>29);       	w[7]=t&lower29bits;
	t=(w[8]*4)+(t>>29);       	w[8]=t&lower29bits;
    t=(w[9]*4)+(t>>29);       	w[9]=t&lower29bits;
	t=(w[10]*4)+(t>>29);        w[10]=t&lower29bits;
	t=(w[11]*4)+(t>>29);      	w[11]=t&lower29bits;
	t=(w[12]*4)+(t>>29);      	w[12]=t&lower29bits;
	t=(w[13]*4)+(t>>29);      	w[13]=t&lower29bits;
	t=(w[14]*4)+(t>>29);      	w[14]=t&lower29bits;
	t=(w[15]*4)+(t>>29);      	w[15]=t&lower29bits;
	t=(w[16]*4)+(t>>29);        w[16]=t&lower29bits;
	t=(w[17]*4)+(t>>29);        w[17]=t&lower28bits;           //most significant limb is at most 28-bit

	w[0]+=(t>>28);				//least significant limb could be at most 30-bit
}

// w=2(x-y) and the result in reduced form
void g2sb(int32_t x[], int32_t y[], int32_t w[]){
    int32_t t;
	t=2*(x[0]-y[0]);                 //Could be at most 32-bit
	w[0]=t&lower29bits;

	t=2*(x[1]-y[1])+(t>>29);              w[1]=t&lower29bits;
	t=2*(x[2]-y[2])+(t>>29);              w[2]=t&lower29bits;
	t=2*(x[3]-y[3])+(t>>29);              w[3]=t&lower29bits;
	t=2*(x[4]-y[4])+(t>>29);              w[4]=t&lower29bits;
	t=2*(x[5]-y[5])+(t>>29);              w[5]=t&lower29bits;
	t=2*(x[6]-y[6])+(t>>29);              w[6]=t&lower29bits;
	t=2*(x[7]-y[7])+(t>>29);              w[7]=t&lower29bits;
	t=2*(x[8]-y[8])+(t>>29);              w[8]=t&lower29bits;
	t=2*(x[9]-y[9])+(t>>29);              w[9]=t&lower29bits;
	t=2*(x[10]-y[10])+(t>>29);            w[10]=t&lower29bits;
	t=2*(x[11]-y[11])+(t>>29);            w[11]=t&lower29bits;
	t=2*(x[12]-y[12])+(t>>29);            w[12]=t&lower29bits;
	t=2*(x[13]-y[13])+(t>>29);            w[13]=t&lower29bits;
	t=2*(x[14]-y[14])+(t>>29);            w[14]=t&lower29bits;
	t=2*(x[15]-y[15])+(t>>29);            w[15]=t&lower29bits;
	t=2*(x[16]-y[16])+(t>>29);            w[16]=t&lower29bits;
	t=2*(x[17]-y[17])+(t>>29);            w[17]=t&lower28bits;               //Most significant limb is at most 28-bit

	w[0]+=(t>>28);            //Least significant limb could be at most 30-bit

}


// w=3(w-x) and the result in reduced form
void g3sb(int32_t x[], int32_t w[]){
    int32_t t;
	t=3*(w[0]-x[0]);                 //Could be at most 32-bit
	w[0]=t&lower29bits;

	t=3*(w[1]-x[1])+(t>>29);              w[1]=t&lower29bits;
	t=3*(w[2]-x[2])+(t>>29);              w[2]=t&lower29bits;
	t=3*(w[3]-x[3])+(t>>29);              w[3]=t&lower29bits;
	t=3*(w[4]-x[4])+(t>>29);              w[4]=t&lower29bits;
	t=3*(w[5]-x[5])+(t>>29);              w[5]=t&lower29bits;
	t=3*(w[6]-x[6])+(t>>29);              w[6]=t&lower29bits;
	t=3*(w[7]-x[7])+(t>>29);              w[7]=t&lower29bits;
	t=3*(w[8]-x[8])+(t>>29);              w[8]=t&lower29bits;
	t=3*(w[9]-x[9])+(t>>29);              w[9]=t&lower29bits;
	t=3*(w[10]-x[10])+(t>>29);            w[10]=t&lower29bits;
	t=3*(w[11]-x[11])+(t>>29);            w[11]=t&lower29bits;
	t=3*(w[12]-x[12])+(t>>29);            w[12]=t&lower29bits;
	t=3*(w[13]-x[13])+(t>>29);            w[13]=t&lower29bits;
	t=3*(w[14]-x[14])+(t>>29);            w[14]=t&lower29bits;
	t=3*(w[15]-x[15])+(t>>29);            w[15]=t&lower29bits;
	t=3*(w[16]-x[16])+(t>>29);            w[16]=t&lower29bits;
	t=3*(w[17]-x[17])+(t>>29);            w[17]=t&lower28bits;               //Most significant limb is at most 28-bit

	w[0]+=(t>>28);            //Least significant limb could be at most 30-bit

}

// w-=2*x and the result in reduced form
void gsb2(int32_t x[], int32_t w[]){
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
}


// reduce w - Reduced Limb Form (RLF)
/* By testing we have found [0,2^{30}-1]x[0,2^{29}-1]^{16}x[0,2^{28}-1]
 where least significant limb is in [0,2^{30}-1] and the most significant limb is in [0,2^{28}-1]*/
void scr(int32_t w[]){
	int32_t t1, t2;
	t1=w[0]&lower29bits;

	t2=w[1]+(w[0]>>29);         w[1]=t2&lower29bits;
	t2=w[2]+(t2>>29);           w[2]=t2&lower29bits;
	t2=w[3]+(t2>>29);       	w[3]=t2&lower29bits;
	t2=w[4]+(t2>>29);       	w[4]=t2&lower29bits;
	t2=w[5]+(t2>>29);       	w[5]=t2&lower29bits;
	t2=w[6]+(t2>>29);       	w[6]=t2&lower29bits;
	t2=w[7]+(t2>>29);       	w[7]=t2&lower29bits;
	t2=w[8]+(t2>>29);       	w[8]=t2&lower29bits;
	t2=w[9]+(t2>>29);           w[9]=t2&lower29bits;
	t2=w[10]+(t2>>29);      	w[10]=t2&lower29bits;
	t2=w[11]+(t2>>29);      	w[11]=t2&lower29bits;
	t2=w[12]+(t2>>29);      	w[12]=t2&lower29bits;
	t2=w[13]+(t2>>29);      	w[13]=t2&lower29bits;
	t2=w[14]+(t2>>29);         	w[14]=t2&lower29bits;
	t2=w[15]+(t2>>29);      	w[15]=t2&lower29bits;
	t2=w[16]+(t2>>29);      	w[16]=t2&lower29bits;
	t2=w[17]+(t2>>29);          w[17]=t2&lower28bits;           //most significant limb is at most 28-bit
	w[0]=t1+(t2>>28);				//least significant limb could be at most 30-bit

}



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


//Multiplication i.e. z=x*y
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


///////////////// Elliptic Curve sub-Group Operations //////////////////////////////

// Point Structure
typedef struct{
    int32_t x[LIMBS];
    int32_t y[LIMBS];
    int32_t z[LIMBS];
    int32_t inf;
} ECp;

// P=0
void inf(ECp *P){
	P->inf=1;
}

// Initialise P
void init(int32_t *x, int32_t *y, ECp *P){
	for (int i=0;i<LIMBS;i++){
		P->x[i]=x[i];
		P->y[i]=y[i];
		P->z[i]=0;
	}
	P->z[0]=1;
	P->inf=0;
}

// P=Q
void copy(ECp *Q,ECp *P){
	for (int i=0;i<LIMBS;i++){
		P->x[i]=Q->x[i];
		P->y[i]=Q->y[i];
		P->z[i]=Q->z[i];
	}
	P->inf=Q->inf;
}

// P=-Q, The points to be negated are affine, therefore, we ignore the z-coordinate
void neg(ECp *Q,ECp *P){
    for (int32_t i=0;i<LIMBS;i++){
		P->x[i]=Q->x[i];
		P->y[i]=-Q->y[i];
		//P->z[i]=0;
	}
	scr(P->y);
	//P->z[0]=1;
	//P->inf=Q->inf;
}

// From Jacobian to Affine coordinate
void norm(ECp *P){
	int32_t iz2[LIMBS], iz3[LIMBS], w[LIMBS], t[LIMBS];
	if (P->inf) return;
	gcopy(P->z,w);
	ginv(w);
	gsqr(w,iz2);
	TMV_product(iz2,w,iz3);
	TMV_product(P->x,iz2,t);        //t is in intermediate form
	gcopy(t, P->x);
	scr(P->x);                    //P->x is a unique residue in modulo p
	TMV_product(P->y,iz3,t);        //t is in intermediate form
	gcopy(t, P->y);
	scr(P->y);                    //P->y is a unique residue in modulo p
	gone(P->z);
}

// P=2P, Point doubling where P is Affine
void doubl_1(ECp *P){
    int32_t r0[LIMBS],r1[LIMBS],r2[LIMBS],r3[LIMBS];

	gcopy(P->y, P->z);
	gmul2(P->z);

	gsqr(P->x,r0);
	gsqr(P->y,r1);
	gsqr(r1,r2);
	gadd(P->x,r1,r3);
	scr(r3);
	gsqr(r3,r1);
	gtsb(r0,r2,r1);
	gmul2(r1);
	r0[0]-=1;
	gmul3(r0);
	gsqr(r0,P->x);
	gsb2(r1,P->x);

	gsub(r1,P->x,r3);
	scr(r3);
	TMV_product(r3,r0,P->y);
	scr(r2);                  //To make one bit more space for multiplication by 4
	gmul4(r2);
	gsb2(r2,P->y);

}

//Point doubling operation i.e. 2P, where P is Jacobian
void doubl_2(ECp *P){
	int32_t r0[LIMBS],r1[LIMBS],r2[LIMBS],r3[LIMBS], r4[LIMBS];

	gsqr(P->z,r0);
	gsqr(P->y,r1);
	TMV_product(P->x,r1,r2);
	gadd(P->y, P->z, r3);
	gadd(P->x, r0, P->y);
	scr(P->y);

	g3sb(r0,P->x);             //P->x is in reduced form
	TMV_product(P->x,P->y,P->z);
	gsqr(P->z,P->x);
	scr(r2);                  //To make one bit more space for multiplication by 4
	gmul4(r2);
    gcopy(r2, r4);          //4*r2 will be used later
    gsb2(r4, P->x);         //P->x in reduced form

	gsub(r2, P->x, r2);
	scr(r2);
    TMV_product(P->z,r2,P->y);
	gsqr(r1,r2);
	scr(r2);                  //To make one bit more space for multiplication by 4
	gmul4(r2);
    gsb2(r2, P->y);         //P->y is in reduced form

	gsqr(r3,P->z);
	gtsb(r0,r1,P->z);       //P->z is in reduced form
}

// P+=Q, Mixed Point addition where Q is Jacobian and P is Affine. Then the resultant P will be Jacobian
void mixed_add_1(ECp *Q, ECp *P){
    int32_t r0[LIMBS],r1[LIMBS],r2[LIMBS],r3[LIMBS];

	gsqr(Q->z,r0);
	TMV_product(P->x,r0,r1);
	TMV_product(P->y,r0,r2);
	TMV_product(r2,Q->z,P->x);
	gsub(r1,Q->x,r3);
	scr(r3);
	gsqr(r3,r2);
	gadd(Q->z,r3,P->y);
	scr(P->y);
	gsqr(P->y,P->z);
	gtsb(r0,r2,P->z);

	scr(r2);                //To make one bit more space for multiplication by 4
	gmul4(r2);
	TMV_product(r2,r3,r0);
	g2sb(P->x,Q->y,r1);
	gsqr(r1,P->x);
	TMV_product(Q->x,r2,P->y);
	gtsb2(r0,P->y,P->x);

	gcopy(P->y,r2);
	gsub(r2,P->x,r3);
	scr(r3);
	TMV_product(r1,r3,P->y);
	TMV_product(r0,Q->y,r2);
	gsb2(r2,P->y);

}

// Mixed Point Addition i.e. P+=Q, where Q is Affine and P is Jacobian
void mixed_add_2(ECp *Q, ECp *P){
	int32_t r0[LIMBS], r1[LIMBS], r2[LIMBS], r3[LIMBS], r4[LIMBS];

	gsqr(P->z,r3);
	TMV_product(Q->x,r3,r4);
	TMV_product(Q->y,r3,r2);
	TMV_product(r2,P->z,r1);

	gsub(r4,P->x,r4); //r4=h
	g2sb(r1,P->y,r2); //r2=r
	gadd(P->z,r4,r1); //r1=z1+h
	scr(r4);
	gsqr(r4,r0); //r0=hh
	scr(r1);
	gsqr(r1,P->z);
	gtsb(r3,r0,P->z);     //P->z in reduced form

	scr(r0);              //To make one bit more space for multiplication by 4
	gmul4(r0);            //r0 in reduced form, r0=i
	TMV_product(r4,r0,r1); //j
	TMV_product(r1,P->y,r3); //r3=j*y1
 	TMV_product(r0,P->x,P->y); //P->y=v
 	gsqr(r2,P->x); //P->x=r^2
	gcopy(P->y, r4); //r4=v
    gtsb2(r1, r4, P->x);         //P->x is in reduced form

	gsub(P->y,P->x,r0);
	scr(r0);
	TMV_product(r0,r2,P->y);
	gsb2(r3,P->y);              //P->y is in reduced form
}

//Point Addition i.e. P+=Q, where both are Jacobian
void Jac_add(ECp *Q,ECp *P){
	int32_t z1z1[LIMBS], z2z2[LIMBS], u1[LIMBS], u2[LIMBS], s1[LIMBS], s2[LIMBS], h[LIMBS], i[LIMBS], j[LIMBS];

	gsqr(P->z,z1z1);
	gsqr(Q->z,z2z2);
	TMV_product(P->x,z2z2,u1);
	TMV_product(Q->x,z1z1,u2);
	TMV_product(P->y,Q->z,h);
	TMV_product(h,z2z2,s1);
	TMV_product(Q->y,P->z,h);
	TMV_product(h,z1z1,s2);

	gsub(u2,u1,h);
	scr(h);
	gcopy(h,j);
	gmul2(j);
	gsqr(j,i);
	TMV_product(h,i,j);
	TMV_product(u1,i,u2); //v
	g2sb(s2,s1,u1); //r
	gsqr(u1,P->x);
	gtsb2(j, u2, P->x);         //P->x is in reduced form

	TMV_product(s1,j,i);
	gsub(u2,P->x,u2);
	scr(u2);
	TMV_product(u2,u1,P->y);
	gsb2(i,P->y);           //P->y is in reduced limb form

	gadd(P->z,Q->z,i);
	scr(i);
	gsqr(i,j);
	gtsb(z1z1,z2z2,j);
	TMV_product(h,j,P->z);         //P->z is in reduced limb form
}


// Printing the Affine coordinates of the computed new point on the screen
void output(ECp *P){
    int32_t i;
    for(i=0; i<LIMBS; i++)
        printf("x[%d]= %X\n", i, P->x[i]);

	puts("");
	for(i=0; i<LIMBS; i++)
        printf("y[%d]= %X\n", i, P->y[i]);
}


/* Normalise all of P[i] using one inversion - Montgomery's trick
 Assume P[0] is already in Affine form */
void multi_norm(ECp P[]){
	int32_t i;
	int32_t t1[LIMBS], t2[LIMBS], t3[LIMBS], w[M][LIMBS];
	gone(w[1]);   // 0->1
	gcopy(P[1].z,w[2]); // 0-1, 1-2
	for (i=3;i<M;i++)   // 2-3
		TMV_product(w[i-1],P[i-1].z,w[i]);

	TMV_product(w[M-1],P[M-1].z,t1);
	ginv(t1);

	gcopy(P[M-1].z,t2);
	TMV_product(w[M-1],t1,t3);
	gcopy(t3,w[M-1]);

	for (i=M-2;;i--){
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

    for (i=1;i<M;i++)  // 0-1
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
void precomp(ECp *P,ECp W[]){
	ECp Q;
	copy(P,&Q);
	doubl_1(&Q);            //because Q is affine
	copy(P,&W[0]);
	copy(&W[0],&W[1]);
	mixed_add_1(&Q,&W[1]);  //Q is Jacobian and W[1] is Affine

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


/* Windows of width 4-6 */
void window(ECp *Q, ECp *P){
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
	mixed_add_2(Q,P);
}

/*Constant time table look-up - borrowed from ed25519*/
void fe_cmov(int32_t f[],int32_t g[],int32_t ib){
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

static void cmov(ECp *w, ECp *u, int32_t b){
  fe_cmov(w->x,u->x,b);
  fe_cmov(w->y,u->y,b);
}

// return 1 if b==c, no branching
static int32_t equal(const int32_t b, const int32_t c){
	int32_t x=b^c;
	x-=1;  // if x=0, x now -1
	return ((x>>31)&1);
}

static void ct_select(ECp *T, ECp W[], int32_t b){
  ECp MT;
  int32_t m=b>>31;
  int32_t babs=(b^m)-m;

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


// Point Multiplication - scalar is 521 bits
void scalar_mul(int32_t *w, ECp *P){
	ECp W[M], Q;
	precomp(P,W);
	multi_norm(W);
	copy(&W[(w[PANES-1]-1)/2], P);

	for(int32_t i=PANES-2; i>=0; i--){
		ct_select(&Q, W, w[i]);
		window(&Q, P);
	}
	norm(P);

}


///////////////// Main Function //////////////////////////////
int main(void){
	uint64_t bef, aft, min_ccycle=10000000000;
	int32_t i, j, w[PANES], lpz=10000;
	ECp P;
	int32_t xs[LIMBS], ys[LIMBS];

/* Base point on NIST Curve */

	xs[0]=0x2e5bd66;        xs[1]=0xbf3f18e;        xs[2]=0x1a90a6fe;       xs[3]=0x1167830a;       xs[4]=0x1a8de334;       xs[5]=0x93d17f;
	xs[6]=0x4a3f877;        xs[7]=0xefdfceb;        xs[8]=0x1aa14b5e;   	xs[9]=0x35a69ed;        xs[10]=0x1e0a2bd8;      xs[11]=0xa7f6a43;
	xs[12]=0x6481390;       xs[13]=0xada214e;       xs[14]=0x1b2d988e;  	xs[15]=0x1d39b3c7;      xs[16]=0x6b70404;       xs[17]=0x6342c70;

	ys[0]=0x1fd16650;       ys[1]=0x5f4a3b4;        ys[2]=0x1cb09022;   	ys[3]=0x18e10d44;       ys[4]=0x10761353;       ys[5]=0x1c809fd6;
	ys[6]=0x19031542;       ys[7]=0x132bde84;       ys[8]=0xc97ee72;    	ys[9]=0x1939f331;       ys[10]=0x5ebef45;       ys[11]=0xf3688d0;
	ys[12]=0xf544495;       ys[13]=0x1e8deccc;      ys[14]=0x97ed0b1;   	ys[15]=0x18008b91;      ys[16]=0xa789a3b;       ys[17]=0x8c1c94b;

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

	for(i=0;i<10;i++){
		bef=rdtsc();
		for (j=0; j<lpz; j++){
			init(xs,ys,&P);
			scalar_mul(w,&P);
		}
		aft=rdtscp();
		if((aft-bef)/(lpz) < min_ccycle)
				min_ccycle= (aft-bef)/(lpz);
        }
	printf("Window width :: %u\n%s%"PRIu64"\n", WINDOW, "The minimum clock cycles count is :: ", min_ccycle);
	output(&P);

}
