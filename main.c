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

