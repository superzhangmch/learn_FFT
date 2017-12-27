#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

void p(uint32_t * aa, int aa_s, int radix, bool show_all=false)
{
    for (int i = aa_s - 1; i>=0; --i) {
        bool s = true;
        if (show_all) {
            printf("%u", aa[i]);
        } else {
            if (i > aa_s - 100 && i > 100) printf("%u", aa[i]);
            else if (i == 100) {
                printf("...");
            } else if (i < 100) printf("%u", aa[i]);
            else {s = false;}
        }
        if (radix > 10 && s) {
            printf(",");
        }
    }
    printf("(%d)\n", aa_s);
}

// reverse the k-bit representation of t
int rev_bin(int t, int k)
{
    int ret = 0;
    for (int i = 0; i < k; ++i) {
        int a = (1<<i) & t;
        if (a != 0) {
            ret |= (1<<(k - i-1));
        }
    }
    return ret;
}

// 弃九法验证计算的正误
int qi9_check(uint32_t * aa, int aa_s, int radix, int show=0)
{
    int qi9_a=0;
    if (show)printf("\n");
    for (int i = aa_s - 1; i>=0; --i) {
        qi9_a += aa[i];
        qi9_a = qi9_a % 9;
        if(show)printf("%u", aa[i]);
    }
    if(show)printf("\n");
    return qi9_a;
}

