/*The prime order of the group, ECsGord, is taken from FIPS PUB 186-4
ECsGord=6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449
bitlength(ECsGord)= t=521
We have implemented the Algorithm 7, Fixed-base scalar multiplication, proposed by J. W. Bos et al.
The Point operation formulas are selected from Explicit Formulas Database https://www.hyperelliptic.org/EFD/index.html
We are counting the clock cycles of the Evaluation Stage of Algorithm 7 by using the fastest explicit formulas for (Jacobian) doubling and (mixed) addition.
The program is tested 10 times for 10000 iteration in order to obtain consistent cycle count
gcc -Wall -m32 -O3 fb32_ws521.c -o fb32_ws521.exe */


#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>
#include<math.h>


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
   __asm__ __volatile__ ("rdtscp": "=a" (lo), "=d" (hi) :: "%rcx");
   __asm__ __volatile__ ("cpuid"::: "%rax", "%rbx", "%rcx", "%rdx");
   return (uint64_t)hi << 32 | lo;
}

///////////////////////// Field Operations ////////////////////////////////////
// w=1
void gone(int32_t w[]){
	w[0]=1;
	w[1]= w[2]= w[3]= w[4]= w[5]= w[6]= w[7]= w[8]= w[9]= 0;
	w[10]= w[11]= w[12]= w[13]= w[14]= w[15]= w[16]= w[17]=0;
}

// w=x+y and the result in reduced form
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

// w=x-y and the result in reduced form
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


//w=w-x-y and the result in reduce limb form
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



// w*=4 and the result in reduced form
void gmul4(int32_t w[]){
	int32_t t;
	t=w[0]*4;		            //Could be at most 31-bit
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

//Squaring i.e. z=x^2

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

//
// Inverse x = 1/x = x^(p-2) mod p
// 13 muls, 520 sqrs
//
void ginv(int32_t *x){
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


///////////////// Curve Operations //////////////////////////////

// Point Structure
typedef struct{
    int32_t x[LIMBS];
    int32_t y[LIMBS];
    int32_t z[LIMBS];
    int32_t inf;
} ECp_Jac;

typedef struct{
    int32_t x[LIMBS];
    int32_t y[LIMBS];
} ECp_Aff;

// P=0
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


// P=Q where both P and Q are Affine
void Aff_copy_Jac(ECp_Aff *Q, ECp_Jac *P){
	for (int i=0;i<LIMBS;i++){
		P->x[i]=Q->x[i];
		P->y[i]=Q->y[i];
		P->z[i]=0;
	}
	P->z[0]=1;
	P->inf=0;
}

// P=Q where both P and Q are Affine
void Jac_copy(ECp_Jac *Q, ECp_Jac *P){
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
	scr(P->y);
}

// From Jacobian to Affine coordinate
void norm(ECp_Jac *P){
	int32_t iz2[LIMBS], iz3[LIMBS], w[LIMBS], t[LIMBS];
	if (P->inf) return;
	gcopy(P->z,w);
	ginv(w);
	gsqr(w,iz2);
	TMV_product(iz2,w,iz3);
	TMV_product(P->x,iz2,t);        //t is in intermediate form
	gcopy(t, P->x);
	scr(P->x);                   //P->x is a unique residue in modulo p
	TMV_product(P->y,iz3,t);        //t is in intermediate form
	gcopy(t, P->y);
	scr(P->y);                    //P->y is a unique residue in modulo p
	gone(P->z);
}


//Point doubling operation i.e. 2P, where P is Jacobian
void doubl_Jac(ECp_Jac *P){
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

// P+=Q, where P and Q are Jacobian
void add_Jac(ECp_Jac *Q, ECp_Jac *P){
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
	gcopy(h,j);
	gmul2(j);
	gsqr(j,i);
	scr(h);
	TMV_product(h,i,j);

	TMV_product(u1,i,u2);
	g2sb(s2,s1,u1);
	gsqr(u1,i);
	gcopy(i, P->x);
	gtsb2(j, u2, P->x);

	TMV_product(s1,j,i);
	gsub(u2,P->x,j);
	scr(j);
	TMV_product(j,u1,P->y);
	gsb2(i,P->y);

	gadd(P->z,Q->z,i);
	scr(i);
	gsqr(i,j);
	gtsb(z1z1,z2z2,j);
	TMV_product(h,j,P->z);
}

// Mixed Point Addition i.e. P+=Q, where Q is Affine and P is Jacobian
void add_mixed(ECp_Aff *Q, ECp_Jac *P){
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

// Normalise all of P[i] using one inversion - Montgomery's trick
// Assume P[0] is already Affine and table parameter V=3
void multi_norm(ECp_Jac P[][V], ECp_Aff Q[][V]){
	int i, j, sz=M*V;
	int32_t t1[LIMBS], t2[LIMBS], w[sz][LIMBS];
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


// Printing the Affine coordinates of the computed new point on the screen
void output(ECp_Jac *P){
    int32_t i;
    for(i=0; i<LIMBS; i++)
        printf("x[%d]= %X\n", i, P->x[i]);

	puts("");
	for(i=0; i<LIMBS; i++)
        printf("y[%d]= %X\n", i, P->y[i]);

    puts("");

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
    multi_norm(Jac_PreT, W);
}



/*Constant time table look-up - borrowed from ed25519*/
void fe_cmov(int32_t f[], int32_t g[],int32_t ib){
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

static void cmov(ECp_Aff *w, ECp_Aff *u,int32_t b){
  fe_cmov(w->x,u->x,b);
  fe_cmov(w->y,u->y,b);
}

// Due the sign table we have changed the equality function as follows
// return -1 if b==c, else 0, no branching
static int32_t equal(const uint32_t b, const uint32_t c){
	int32_t x=b^c;
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


///////////////// Main Function //////////////////////////////
int main(void){
	uint64_t bef, aft, min_ccycle=1000000000000;
	int32_t i, j, lpz=10000;
	ECp_Aff PreT[M][V];
	ECp_Jac JP;

/* Base point on NIST Curve */

	JP.x[0]=0x2e5bd66;        JP.x[1]=0xbf3f18e;         JP.x[2]=0x1a90a6fe;        JP.x[3]=0x1167830a;       JP.x[4]=0x1a8de334;        JP.x[5]=0x93d17f;
	JP.x[6]=0x4a3f877;        JP.x[7]=0xefdfceb;         JP.x[8]=0x1aa14b5e;    	JP.x[9]=0x35a69ed;        JP.x[10]=0x1e0a2bd8;       JP.x[11]=0xa7f6a43;
	JP.x[12]=0x6481390;       JP.x[13]=0xada214e;        JP.x[14]=0x1b2d988e;   	JP.x[15]=0x1d39b3c7;      JP.x[16]=0x6b70404;        JP.x[17]=0x6342c70;

	JP.y[0]=0x1fd16650;       JP.y[1]=0x5f4a3b4;         JP.y[2]=0x1cb09022;    	JP.y[3]=0x18e10d44;       JP.y[4]=0x10761353;        JP.y[5]=0x1c809fd6;
	JP.y[6]=0x19031542;       JP.y[7]=0x132bde84;        JP.y[8]=0xc97ee72;     	JP.y[9]=0x1939f331;       JP.y[10]=0x5ebef45;        JP.y[11]=0xf3688d0;
	JP.y[12]=0xf544495;       JP.y[13]=0x1e8deccc;       JP.y[14]=0x97ed0b1;    	JP.y[15]=0x18008b91;      JP.y[16]=0xa789a3b;        JP.y[17]=0x8c1c94b;

	JP.z[0]=0x1;              JP.z[1]=0x0;              JP.z[2]=0x0;    	        JP.z[3]=0x0;            JP.z[4]=0x0;                JP.z[5]=0x0;
	JP.z[6]=0x0;              JP.z[7]=0x0;              JP.z[8]=0x0;     	        JP.z[9]=0x0;            JP.z[10]=0x0;               JP.z[11]=0x0;
	JP.z[12]=0x0;             JP.z[13]=0x0;             JP.z[14]=0x0;    	        JP.z[15]=0x0;           JP.z[16]=0x0;               JP.z[17]=0x0;

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
	printf("Window width :: %u\n%s%"PRIu64"\n", WINDOW, "The minimum clock cycles count is :: ", min_ccycle);
	output(&JP);

}
