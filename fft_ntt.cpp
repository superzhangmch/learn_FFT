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

FntMul fnt;
FftMul fft;

extern "C" {

    void fnt_init(int max_digit) {
        fnt.init(max_digit * 2);
    }

    void fft_init(int max_digit) {
        fft.init(max_digit * 2);
    }

    int do_fft(uint32_t * aa, int aa_s, uint32_t * bb, int bb_s, int radix, uint32_t * cc)
    {
        int cc_s;
        fft.fast_prod(aa, aa_s, bb, bb_s, radix, cc, cc_s);
        return cc_s;
    }
    int do_fnt(uint32_t * aa, int aa_s, uint32_t * bb, int bb_s, int radix, uint32_t * cc)
    {
        int cc_s;
        fnt.fast_prod(aa, aa_s, bb, bb_s, radix, cc, cc_s);
        return cc_s;
    }
}

