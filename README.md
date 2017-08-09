# Faster-Residue-Multiplication and Application to ECC
We present faster residue multiplication modulo 521-bit Mersenne prime for 32- and 64-bit platforms by using Toeplitz Matrix-Vector Product (TMVP). The total arithmetic cost of our proposed algorithms is less than the existing algorithms. We implemented constant-time variable- and fixed-base scalar multiplication for the standard NIST curve P-521 and Edwards curve E-521.

The constant-time variable-base scalar multiplication using fixed-window 'w' is implemented using the “Algorithm 1” of J. W. Bos et al. paper “Selecting elliptic curves for cryptography: an efficiency and security analysis”. Instead of 2^{w-2}, the size of our pre-computed table is 2^{w-1}. Similarly, the constant-time fixed-base scalar multiplication is implemented using the modified LSB-set comb method algorithm proposed by J. W. Bos et al. i.e. “Algorithm 7”. But the point addition and doubling formulas are selected from the Explicit Formulas Database (EFD), managed by Daniel J. Bernstein and Tanja Lange, and from the “Curve41417: Karatsuba revisited” paper of Bernstein et al.
