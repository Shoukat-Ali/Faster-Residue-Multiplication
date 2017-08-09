/* Test program for ed521.c constant-time variable-base scalar multiplication
 This curve is recommended by Daniel J. Bernstein and Tanja Lange at http://safecurves.cr.yp.to
 Curve Equation :: x^2+y^2 = 1-376014x^2y^2 where d=-376014
 Uses Bernstein et al.s point multiplication method from Curve41417 paper
 Cache safety thanks to ed25519
 Thanks to M. Scott for making his code public and is available at http://indigo.ie/~mscott/ed521.cpp
 We have replaced gmul() by TMV_product()
 We have added rdtscp() to measure the clock cycles count as suggested by Paoloni; http://www.intel.com.tr/content/dam/www/public/us/en/documents/white-papers/ia-32-ia-64-benchmark-code-execution-paper.pdf
 We have amended scr(), gmuli(), gsqr(), and gsqr2() according to modulus p
 gcc -Wall -O3 vb32_ed521.c -o vb32_ed521.exe */
