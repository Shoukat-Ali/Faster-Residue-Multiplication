/* Test program for ed521 scalar point multiplication
 Uses Bernstein et al.s point multiplication method from Curve41417 paper
 Cache safety thanks to ed25519
 Fully Tested and debugged
 We have replaced gmul() by TMV_product()
 We have added rdtscp() to measure the clock cycles count as suggested by Paoloni; http://www.intel.com.tr/content/dam/www/public/us/en/documents/white-papers/ia-32-ia-64-benchmark-code-execution-paper.pdf
 We have amended scr(), gmuli(), gsqr(), and gsqr2() according to modulus p
 gcc -Wall -O3 vb64_ed521.c -o vb64_ed521.exe */
