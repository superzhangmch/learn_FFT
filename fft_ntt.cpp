#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "complex.h"
#include "zp.h"
#include "fft_common.h"

template <typename TComplex>
class ProductFFT {
protected:
    int _max_result_len;
    int _max_k;
    TComplex * _omega_fft[32];
    int * _pi_k[32];

    void prepare_data(int k)
    {
        int n = (int)pow(2, k);
        TComplex * omega_fft = new TComplex[n];
        int * pi_k = new int[n];

        // if (k == 5) printf("%d--\n", k);
        // if (k == 5) printf("%d->%.5f+%.5fi\n", i, w_ret[i][0], w_ret[i][1]);
        TComplex::get_w_pow(n, omega_fft);

        for (int i = 0; i < n; ++i) {
            pi_k[i] = rev_bin(i, k);
            // if (k == 5) printf("%d->%d\n", i, pi_k_ret[i]);
        }
        _omega_fft[k] = omega_fft;
        _pi_k[k] = pi_k;
    }

public:
    ProductFFT() {
        _max_result_len = _max_k = 0;
        for (int i = 0; i < 32; ++i) {
            _omega_fft[i] = NULL;
            _pi_k[i] = NULL;
        }
    }
    virtual ~ProductFFT() {
        for (int i = 0; i < 32; ++i) {
            if (_omega_fft[i]) {
                delete [] _omega_fft[i];
                _omega_fft[i] = NULL;
            }
            if (_pi_k[i]) {
                delete [] _pi_k[i];
                _pi_k[i] = NULL;
            }
        }
    }

    void init(int max_len)
    {
        int k = int(ceil(log(1.L * max_len) / log(2.)));
        _max_result_len = int(round(2 * pow(2., k)));
        _max_k = k + 1;

        for (int i = 1; i < _max_k + 1; ++i) {
            printf("init %d/%d\n", i, _max_k);
            prepare_data(i);
        }
        printf("init done\n");
    }

    virtual void FFT(int k, int n, int * input, TComplex *input1, int input_size, TComplex *transform) = 0;
    virtual void rev_FFT(int k, int n, int * input, TComplex *input1, int input_size, TComplex *transform) = 0;

    int fast_prod(int * aa, int aa_length, int * bb, int bb_length, int radix, int * out, int &out_size) 
    {
        // aa or bb == 0
        if (aa_length == 0 || bb_length == 0) {
            out_size = 0;
            return 0;
        }
        // aa or bb is small
        if (aa_length == 1 || bb_length == 1) {
            int small_val, big_size, *big;
            if (aa_length == 1) {
                small_val = aa[0];
                big = bb;
                big_size = bb_length;
            } else {
                small_val = bb[0];
                big = aa;
                big_size = aa_length;
            }
            if (small_val == 0) {
                out_size = 0;
                return 0;
            }
            int out_len = big_size;
            int remain = 0;
            for (int i = 0; i < big_size; ++i) {
                int res = big[i] * small_val + remain;
                if (res >= radix) {
                    out[i] = res % radix;
                    remain = res / radix;
                } else {
                    out[i] = res;
                    remain = 0;
                }
            }
            if (remain) {
                out[big_size] = remain;
                out_len = big_size + 1;
            }
            out_size = out_len;
            return 0; 
        }

        // aa or bb both big
        int max_len = aa_length > bb_length ? aa_length : bb_length;
        int k = int(ceil(log(max_len) / log(2))) + 1;
        int n2 = int(round(pow(2, k)));
        if (k > _max_k) {
            out_size = 0;
            return -1;
        }
        //printf("--%d %d %d %d %d\n", k, n2, aa_length, bb_length, max_len);
        TComplex *buf = new TComplex[2*n2];
        TComplex *transform_aa = buf;
        TComplex *transform_bb = buf + n2;
        TComplex *transform_cc = buf + n2;
        TComplex *transform_dd = buf;
        FFT(k, n2, aa, NULL, aa_length, transform_aa);
        //printf("xxxx\n\n");
        // pc(transform_aa, n2); printf("-----------------000000000\n");
        FFT(k, n2, bb, NULL, bb_length, transform_bb);
        // pc(transform_bb, n2); printf("-----------------111111111\n");

        for (int i = 0; i < n2; i++) {
            transform_cc[i] = transform_aa[i] * transform_bb[i];
        }

        // pc(transform_cc, n2); printf("-----------------22222222\n");
        rev_FFT(k, n2, NULL, transform_cc, n2, transform_dd);

        int max_res_len = aa_length + bb_length;
        int remain = 0;
        for (int i = 0; i < max_res_len; ++i) {
            int o = transform_dd[i].to_int() + remain;
            if (o >= radix) {
                out[i] = o % radix;
                remain = o / radix;
            } else {
                out[i] = o;
                remain = 0;
            }
        }

        out_size = out[max_res_len-1] != 0 ? max_res_len : max_res_len - 1;
        delete [] buf;
        return 0;
    }
};

class FftMul: public ProductFFT<Complex> {

public:
    void FFT(int k, int n, int * input, Complex *input1, int input_size, Complex *transform)
    {
        Complex *omega = _omega_fft[k];
        int * pi_k = _pi_k[k];
        int * P = input;
        Complex * P1 = input1;

        if (P1 == NULL) {
            for (int t = 0; t < n-1; t+= 2) {
                int P_pi_k_t = pi_k[t] < input_size ? P[pi_k[t]] : 0;
                int P_pi_k_t_1 = pi_k[t+1] < input_size ? P[pi_k[t+1]] : 0;
                transform[t]   = P_pi_k_t + P_pi_k_t_1;
                transform[t+1] = P_pi_k_t - P_pi_k_t_1;
            }
        } else {
            for (int t = 0; t < n-1; t+= 2) {
                Complex P_pi_k_t = pi_k[t] < input_size ? P1[pi_k[t]] : 0;
                Complex P_pi_k_t_1 = pi_k[t+1] < input_size ? P1[pi_k[t+1]] : 0;

                transform[t]   = P_pi_k_t + P_pi_k_t_1;
                transform[t+1] = P_pi_k_t - P_pi_k_t_1;
            }
        }

        int m = n/2;
        int num = 2;
        int pow_d_2 = (int)round(pow(2, k-2));
        for (int d = k-2; d > -1; d--) {
            m /= 2;
            num *= 2;
            int Max = (pow_d_2-1)*num + 1;
            pow_d_2 /= 2;

            for (int t = 0; t < Max; t += num) {
                for (int j = 0; j < num / 2; ++j) {
                    Complex xPOdd     = omega[m*j] * transform[t+num/2+j];
                    Complex prevTrans = transform[t+j];
                    transform[t+j]       = prevTrans + xPOdd;
                    transform[t+num/2+j] = prevTrans - xPOdd;

                    // omega[m*j] == - omega[m*(j+num/2)], 故这里可节约一次乘法
                    // 如果如下面这样展开了(NTT种就没法节约一次乘法), 实测运行时间翻倍
                    /*
                    Complex xPOdd     = omega[m*j] * transform[t+num/2+j];
                    Complex xPOdd1    = omega[m*(j+num/2)] * transform[t+num/2+j];
                    Complex prevTrans = transform[t+j];
                    transform[t+j]       = prevTrans + xPOdd;
                    transform[t+num/2+j] = prevTrans + xPOdd1;
                    */
                }
            }
        }
    }

    void rev_FFT(int k, int n, int * input, Complex *input1, int input_size, Complex *transform)
    {
        FFT(k, n, input, input1, input_size, transform);
        for (int i = n-1; i > n/2-1; i--) {
            Complex a = transform[n - i] / n;
            transform[n-i] = transform[i] / n;
            transform[i] = a;
        }
        transform[0] /= n;
    }
};


class NttMul: public ProductFFT<Zp> {
    /*
     *  NTT 与 FFT 两种变换的框架是一样的。但是由于复根的特殊性可以利用，因此这里各实现一下
     *  这里对100万位的乘法实测得：NTT比FFT慢，NTT速度是FFT的1/5
     * */

public:
    void FFT(int k, int n, int * input, Zp *input1, int input_size, Zp *transform)
    {
        Zp *omega = _omega_fft[k];
        int * pi_k = _pi_k[k];
        int * P = input;
        Zp * P1 = input1;

        if (P1 == NULL) {
            for (int t = 0; t < n-1; t+= 2) {
                int P_pi_k_t = pi_k[t] < input_size ? P[pi_k[t]] : 0;
                int P_pi_k_t_1 = pi_k[t+1] < input_size ? P[pi_k[t+1]] : 0;
                transform[t]   = P_pi_k_t;
                transform[t+1] = P_pi_k_t_1;
            }
        } else {
            for (int t = 0; t < n-1; t+= 2) {
                Zp P_pi_k_t = pi_k[t] < input_size ? P1[pi_k[t]] : 0;
                Zp P_pi_k_t_1 = pi_k[t+1] < input_size ? P1[pi_k[t+1]] : 0;

                transform[t]   = P_pi_k_t;
                transform[t+1] = P_pi_k_t_1;
            }
        }

        int m = n/1;
        int num = 1;
        int pow_d_2 = (int)round(pow(2, k-1));
        for (int d = k-1; d > -1; d--) {
            m /= 2;
            num *= 2;
            int Max = (pow_d_2-1)*num + 1;
            pow_d_2 /= 2;
            for (int t = 0; t < Max; t += num) {
                for (int j = 0; j < num / 2; ++j) {
                    Zp xPOdd     = omega[m*j] * transform[t+num/2+j];
                    Zp xPOdd1    = omega[m*(j+num/2)] * transform[t+num/2+j];
                    Zp prevTrans = transform[t+j];
                    transform[t+j]       = prevTrans + xPOdd;
                    transform[t+num/2+j] = prevTrans + xPOdd1;
                }
            }
        }
    }

    void rev_FFT(int k, int n, int * input, Zp *input1, int input_size, Zp *transform)
    {
        FFT(k, n, input, input1, input_size, transform);
        // 必须除以mod P 意义上的n. 这里先求倒数，再作乘积
        Zp n_reciprocal = Zp(1) / Zp(n);
        for (int i = n-1; i > n/2-1; i--) {
            Zp a = transform[n - i] * n_reciprocal;
            transform[n-i] = transform[i] * n_reciprocal;
            transform[i] = a;
        }
        transform[0] *= n_reciprocal;
    }
};

//long Zp::P = 1004535809;
//long Zp::G = 3;
long Zp::P = 2113929217;
long Zp::G = 5;
int main()
{
    NttMul fft;
    int max_digit = 200000;
    fft.init(max_digit * 4);
    srand(time(NULL));

    int *aa = new int[max_digit];
    int *bb = new int[max_digit];
    int *cc = new int[max_digit * 2];
    int aa_s, bb_s, cc_s = 0;
    int radix = 1000; // for radix=100, max_digit = 200w not ok

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
        fft.fast_prod(aa, aa_s, bb, bb_s, radix, cc, cc_s); 
        gettimeofday(&tpend,NULL);

        double tm = ((tpend.tv_sec-tpstart.tv_sec)*1000000+(tpend.tv_usec-tpstart.tv_usec));
        int qi9_a=0, qi9_b=0, qi9_c = 0;
        qi9_a = qi9_check(aa, aa_s, radix);
        qi9_b = qi9_check(bb, bb_s, radix);
        qi9_c = qi9_check(cc, cc_s, radix);
        printf("len=%d time=%.4f, check=%d*%d=%d ~ %d\n", cc_s, tm/1000000, qi9_a, qi9_b, (qi9_b * qi9_a) % 9, qi9_c);
        if ((qi9_b * qi9_a) % 9 != qi9_c)
        {
            printf("not match\n");
            p(aa, aa_s, radix);
            p(bb, bb_s, radix);
            p(cc, cc_s, radix);
            return -1;
        }
        for (int i = cc_s - 1; i >= 0; --i) {
            // printf("%d ", cc[i]);
        }
        //printf("\n%d\n", cc_s);
    }
    return 0;
}
