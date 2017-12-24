#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void p(int * aa, int aa_s, bool show_all=false)
{
    for (int i = aa_s - 1; i>=0; --i) {
        if (show_all) {
            printf("%d", aa[i]);
        } else {
            if (i > aa_s - 100 && i > 100) printf("%d", aa[i]);
            if (i == 100) printf("...");
            if (i < 100) printf("%d", aa[i]);
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
int qi9_check(int * aa, int aa_s, int radix, int show=0)
{
    int qi9_a=0;
    if (show)printf("\n");
    for (int i = aa_s - 1; i>=0; --i) {
        qi9_a += aa[i];
        qi9_a = qi9_a % 9;
        if(show)printf("%d", aa[i]);
    }
    if(show)printf("\n");
    return qi9_a;
}

