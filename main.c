// In this version we have computed P2, P3 and P4 using TMVP as signed integer both from input and output parameters.
// P1, P6 and P5 are computed using schoolbook such that both the multiplicands of P_1 are taken signed. While P_5 and P_6 are unsigned

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>
#include <gmp.h>

#define ITERATIONS 10000000
#define LIMBS 18


static const int32_t lower29bits = 0x1fffffff;
static const int32_t lower28bits = 0xfffffff;


__inline__ uint64_t rdtsc(){
   uint32_t lo, hi;
   __asm__ __volatile__ ("xorl %%eax,%%eax \n        cpuid"::: "%rax", "%rbx", "%rcx", "%rdx");
   __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
   return (uint64_t)hi << 32 | lo;
}

__inline__ uint64_t rdtscp(){
   uint32_t lo, hi;
   __asm__ __volatile__ ("rdtscp" : "=a" (lo), "=d" (hi) :: "%rcx");
   __asm__ __volatile__ ("cpuid" ::: "%rax", "%rbx","%rcx", "%rdx");
   return (uint64_t)hi << 32 | lo;
}

/*The function computes TMVP of size 6 by using the 2-way decomposition formula in our paper.
mat = is one-dimensional array of 11 elements because it's Toeplitz
vec = is one-dimensional array of size 6
res = is one-dimensional array of size 6 such that it represents the resultant product/vector */
inline void TMVP_level2(const int32_t mat[], const int32_t vec[], int64_t res[]){
    int32_t temp[4];

    temp[0]= vec[0]-vec[3];         temp[1]= vec[1]-vec[4];             temp[2]= vec[2]-vec[5];

    res[3]= res[0]= (int64_t)mat[0]*temp[0]+(int64_t)mat[1]*temp[1]+(int64_t)mat[2]*temp[2];
    res[4]= res[1]= (int64_t)mat[6]*temp[0]+(int64_t)mat[0]*temp[1]+(int64_t)mat[1]*temp[2];
    res[5]= res[2]= (int64_t)mat[7]*temp[0]+(int64_t)mat[6]*temp[1]+(int64_t)mat[0]*temp[2];

    temp[0]= mat[3]+mat[0];         temp[1]= mat[4]+mat[1];             temp[2]= mat[2]+mat[6];
    temp[3]= mat[1]+mat[7];

    res[0]+= (int64_t)temp[0]*vec[3]+(int64_t)temp[1]*vec[4]+(int64_t)(mat[5]+mat[2])*vec[5];
    res[1]+= (int64_t)temp[2]*vec[3]+(int64_t)temp[0]*vec[4]+(int64_t)temp[1]*vec[5];
    res[2]+= (int64_t)temp[3]*vec[3]+(int64_t)temp[2]*vec[4]+(int64_t)temp[0]*vec[5];

    temp[0]= mat[8]+mat[0];         temp[1]= temp[3];                   temp[3]= mat[9]+mat[6];

    res[3]= (int64_t)temp[0]*vec[0]+(int64_t)temp[1]*vec[1]+(int64_t)temp[2]*vec[2] - res[3];
    res[4]= (int64_t)temp[3]*vec[0]+(int64_t)temp[0]*vec[1]+(int64_t)temp[1]*vec[2] - res[4];
    res[5]= (int64_t)(mat[10]+mat[7])*vec[0]+(int64_t)temp[3]*vec[1]+(int64_t)temp[0]*vec[2] - res[5];

}


/*The function computes the Z = X*Y mod (2**521) - 1] using TMVP where X, Y, and Z are 18-limb.
X_mat and Y_vec are one-dimensional arrays of integers where the first element is at most 30-bit, the next 16 elements are 29-bit at most and the last element is 28-bit integer at most.
The output Z_res is one-dimensional array of integers such that Z_res[0] in [0,2**30-1], Z_res[i] in [0,2**29-1] and Z_res[17] in [0,2**28-1]*/
void TMVP_level1(const int32_t X_mat[], const int32_t Y_vec[], int32_t Z_res[]){
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


/*void scr(int32_t w[]){
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

}*/


//////////////////////////Timing Test /////////////////////////////////////////////////////////
/*This function computes the residue multiplication modulo 521-bit Mersenne prime
 and counts the clock cycles of the function by average the cycles count for 10000000 calls.
 Remember: size(X)=size(Y)=18, no matter what their binary length may be.
*/
int main(){
    int32_t X[LIMBS], Z[LIMBS], Y[LIMBS];
    uint64_t start, end, min_ccycle=100000;
    uint32_t i, j, values=100;

    srand(time(NULL));
	for(j=0; j<values; j++){
        for(i=0; i<LIMBS-1; i++){
            X[i] = (rand()<<15)%(1<<29);
            Y[i] = (rand()<<15)%(1<<29);
        }
        X[i] = (rand()<<13)%(1<<28);          //Most Significant Limb is 28-bit i.e X[17]
        Y[i] = (rand()<<13)%(1<<28);          //Most Significant Limb is 28-bit i.e Y[17]

        start=rdtsc();
        end=rdtscp();
        start=rdtsc();
        for (i=0; i<ITERATIONS; i++)
            TMVP_level1(X, Y, Z);
        end=rdtscp();
        if((end-start)/ITERATIONS < min_ccycle)
            min_ccycle = (end-start)/ITERATIONS;
        if(j%(values/10) == 0)
            printf("%u\n", j);

        }

	printf("The minimum number of clock cycles is :: %"PRIu64"\n", min_ccycle);

}

///////////////////////Checking for correctness /////////////////////////////////////////////
/*int main(){
    int32_t X[LIMBS], Z[LIMBS], Y[LIMBS];
    uint32_t i, j, values=1000000;
    mpz_t term, X_mpz, Y_mpz, Z_mpz, modulus, prod_mpz;

    mpz_init(term);
    mpz_init(X_mpz);
    mpz_init(Y_mpz);
    mpz_init(Z_mpz);
    mpz_init(prod_mpz);

    mpz_init(modulus);
    mpz_set_ui(modulus, 1);
    mpz_mul_2exp(modulus, modulus, 521);
    mpz_sub_ui(modulus, modulus, 1);
    //gmp_printf("Modulus :: %Zd\n", modulus);

    srand(time(NULL));
	for(j=0; j<values; j++){
        mpz_set_ui(X_mpz, 0);
        mpz_set_ui(Y_mpz, 0);
        mpz_set_ui(Z_mpz, 0);
        for(i=0; i<LIMBS-1; i++){
            X[i] = (rand()<<15)%(1<<29);
            Y[i] = (rand()<<15)%(1<<29);
            mpz_set_ui(term, X[i]);
			mpz_mul_2exp(term, term, 29*i);
			mpz_add(X_mpz, X_mpz, term);
			mpz_set_ui(term, Y[i]);
			mpz_mul_2exp(term, term, 29*i);
			mpz_add(Y_mpz, Y_mpz, term);
        }
        X[i] = (rand()<<13)%(1<<28);          //Most Significant Limb is 28-bit i.e X[17]
        //printf("%d\t", X[i]);
        Y[i] = (rand()<<13)%(1<<28);          //Most Significant Limb is 28-bit i.e Y[17]
        mpz_set_ui(term, X[i]);
		mpz_mul_2exp(term, term, 29*i);     //Remember the term (position/index) does not change and that's we have 29*17=493
		mpz_add(X_mpz, X_mpz, term);
		mpz_set_ui(term, Y[i]);
		mpz_mul_2exp(term, term, 29*i);     //Remember the term (position/index) does not change and that's we have 29*17=493
		mpz_add(Y_mpz, Y_mpz, term);

        TMVP_level1(X, Y, Z);
        if(Z[0] >= (1<<29))
            scr(Z);
        for(i=0; i<LIMBS; i++){
            mpz_set_ui(term, Z[i]);
			mpz_mul_2exp(term, term, 29*i);
			mpz_add(Z_mpz, Z_mpz, term);
        }
        mpz_mul(prod_mpz, X_mpz, Y_mpz);
        mpz_mod(prod_mpz, prod_mpz, modulus);     //Storing the remainder in prod_mpz
        if(mpz_cmp(Z_mpz, prod_mpz) != 0){
            printf("Not equal at %u\n", j);
            gmp_printf("X :: %Zd\n", X_mpz);
            gmp_printf("Y :: %Zd\n", Y_mpz);
        }
        if(j%(values/10) == 0)
            printf("%u\n", j);
    }

}*/
