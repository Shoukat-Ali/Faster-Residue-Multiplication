// Test program for ws521 scalar point multiplication
// This is NIST standard Weierstrass curve p-521
// Fully Tested and debugged
// Uses constant time method described by Bos et al. - http://eprint.iacr.org/2014/130
// Cache safety thanks to ed25519
// We have replaced gmul() by TMV_product()
// We have added rdtscp() to measure the clock cycles count as suggested by Paoloni; http://www.intel.com.tr/content/dam/www/public/us/en/documents/white-papers/ia-32-ia-64-benchmark-code-execution-paper.pdf
// We have amended scr() and gsqr() according to modulus p
// gcc -Wall -O3 vb64_ws521.c -o vb64_ws521.exe

