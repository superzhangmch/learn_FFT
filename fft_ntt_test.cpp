#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdint.h>
#include "fft_ntt.h"

template <> long Zp<0>::P = 2113929217;
template <> long Zp<0>::G = 5;
template <> long Zp<1>::P = 2013265921;
template <> long Zp<1>::G = 31;
template <> long Zp<2>::P = 1811939329;
template <> long Zp<2>::G = 13;


int main()
{
    FntMul fft;
    int max_digit = 112000 * 2;
    int radix = 1000000000; // for radix=100, max_digit = 200w not ok

    //FftMul fft;
    //int max_digit = 200000*9;
    //int radix = 10; // for radix=100, max_digit = 200w not ok
    fft.init(max_digit * 4);
    srand(time(NULL));

    uint32_t *aa = new uint32_t[max_digit];
    uint32_t *bb = new uint32_t[max_digit];
    uint32_t *cc = new uint32_t[max_digit * 2];
    int aa_s, bb_s, cc_s = 0;

    struct timeval tpstart, tpend;
    aa_s = bb_s = max_digit/2;

    for (int k = 0; k < 100; ++k) {
        for (int i = 0; i< aa_s; ++i) {
            aa[i] = rand() % radix;
            bb[i] = rand() % radix;
            if (i == aa_s -1) {
                if (aa[i] == 0) aa[i]++;
                if (bb[i] == 0) bb[i]++;
            }
        }

        gettimeofday(&tpstart,NULL);
        fft.fast_prod(aa, aa_s, bb, bb_s, 0, radix, cc, cc_s); 
        gettimeofday(&tpend,NULL);

        double tm = ((tpend.tv_sec-tpstart.tv_sec)*1000000+(tpend.tv_usec-tpstart.tv_usec));
        int qi9_a=0, qi9_b=0, qi9_c = 0;
        qi9_a = qi9_check(aa, aa_s, radix);
        qi9_b = qi9_check(bb, bb_s, radix);
        qi9_c = qi9_check(cc, cc_s, radix);
        printf("len=%d time=%.4f, check=%d*%d=%d ~ %d\n", cc_s, tm/1000000, qi9_a, qi9_b, (qi9_b * qi9_a) % 9, qi9_c);
        if (0) //(qi9_b * qi9_a) % 9 != qi9_c)
        {
            printf("not match\n");
            p(aa, aa_s, radix);
            p(bb, bb_s, radix);
            p(cc, cc_s, radix);
            return -1;
        }
        //printf("\n%d\n", cc_s);
//        break;
    }
    return 0;
}
