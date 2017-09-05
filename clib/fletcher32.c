/* 
fletcher32 check sum of an 16-byte integer array.

Taken verbatim from the optimized c version  as reported in
https://en.wikipedia.org/wiki/Fletcher%27s_checksum#Optimizations
except for the fact that the data size (words) is not an argument 
of the function but rather the address of its location is, as this
is what a fortran call passes as reference. 
the type of words and ndata variables is also defined as uint32_t 
instead of size_t because this is what fortran implicitely does.

SdG,  September 3rd 2017
*/

#include <unistd.h>
#include <stdio.h>
#include <stdint.h>

uint32_t fletcher32( uint16_t const *data, uint32_t *ndata )
{
        uint32_t sum1 = 0xffff, sum2 = 0xffff;
        size_t tlen; uint32_t words = *ndata ;

        while (words) {
                tlen = ((words >= 359) ? 359 : words);
                words -= tlen;
                do {
                        sum2 += sum1 += *data++;
                        tlen--;
                } while (tlen);
                sum1 = (sum1 & 0xffff) + (sum1 >> 16);
                sum2 = (sum2 & 0xffff) + (sum2 >> 16);
        }  
        /* Second reduction step to reduce sums to 16 bits */
        sum1 = (sum1 & 0xffff) + (sum1 >> 16);
        sum2 = (sum2 & 0xffff) + (sum2 >> 16);
        return (sum2 << 16) | sum1;
}
