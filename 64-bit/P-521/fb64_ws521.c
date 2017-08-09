/*The prime order of the group, ECsGord, is taken from FIPS PUB 186-4
ECsGord=6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449
bitlength(ECsGord)= t=521
We have implemented the Algorithm 7, Fixed-base scalar multiplication, proposed by J. W. Bos et al.
But the curve operations are selected by us such the offline computation is performed using the slow Affine addition and doubling formula
While the Online computation is performed using the fastest explicit formulas for (Jacobian) doubling and (mixed) addition.
The program is tested 10 times for 10000 iteration in order to obtain consistent cycle count
gcc -Wall -O3 fb64_ws521.c -o fb64_ws521.exe */
