#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"
#include "zp.h"
#include "fft_common.h"
#include "uint128t.h"
#include <stdint.h>

#ifndef FFT_NTT_H
#define FFT_NTT_H

typedef Zp<0> Zp0;
typedef Zp<1> Zp1;
typedef Zp<2> Zp2;

template <typename TComplex>
class ProductFFT {
protected:
    int _max_result_len;
    int _max_k;
    TComplex * _omega_fft[32];
    uint32_t * _pi_k[32];

    void prepare_data(int k)
    {
        int n = (int)pow(2, k);
        TComplex * omega_fft = new TComplex[n];
        uint32_t * pi_k = new uint32_t[n];

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

    virtual void FFT(int k, int n, uint32_t * input, TComplex *input1, int input_size, TComplex *transform) = 0;
    virtual void rev_FFT(int k, int n, uint32_t * input, TComplex *input1, int input_size, TComplex *transform) = 0;
    virtual int fft_ntt(int k, int n2, uint32_t * aa, int aa_length, uint32_t * bb, int bb_length, 
                        int radix, uint32_t * out, int &out_size) = 0;

    int fast_prod(uint32_t * aa, int aa_length, uint32_t * bb, int bb_length, int radix, uint32_t * out, int &out_size) 
    {
        // aa or bb == 0
        if (aa_length == 0 || bb_length == 0) {
            out_size = 0;
            return 0;
        }
        // aa or bb is small
        if (aa_length == 1 || bb_length == 1) {
            int big_size;
            uint32_t small_val, *big;
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
                uint64_t res = 1L * big[i] * small_val + remain;
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
        fft_ntt(k, n2, aa, aa_length, bb, bb_length, radix, out, out_size);
        //printf("--%d %d %d %d %d\n", k, n2, aa_length, bb_length, max_len);
        return 0;
    }
};

class FftMul: public ProductFFT<Complex> {

public:
    int fft_ntt(int k, int n2, uint32_t * aa, int aa_length, uint32_t * bb, int bb_length, 
                        int radix, uint32_t * out, int &out_size)
    {
        Complex *buf = new Complex[2*n2];
        Complex *transform_aa = buf;
        Complex *transform_bb = buf + n2;
        Complex *transform_cc = buf + n2;
        Complex *transform_dd = buf;
        FFT(k, n2, aa, NULL, aa_length, transform_aa);
        FFT(k, n2, bb, NULL, bb_length, transform_bb);

        for (int i = 0; i < n2; i++) {
            transform_cc[i] = transform_aa[i] * transform_bb[i];
        }

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

    void FFT(int k, int n, uint32_t * input, Complex *input1, int input_size, Complex *transform)
    {
        Complex *omega = _omega_fft[k];
        uint32_t * pi_k = _pi_k[k];
        uint32_t * P = input;
        Complex * P1 = input1;

        if (P1 == NULL) {
            for (int t = 0; t < n-1; t+= 2) {
                uint32_t P_pi_k_t = pi_k[t] < input_size ? P[pi_k[t]] : 0;
                uint32_t P_pi_k_t_1 = pi_k[t+1] < input_size ? P[pi_k[t+1]] : 0;
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

    void rev_FFT(int k, int n, uint32_t * input, Complex *input1, int input_size, Complex *transform)
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


template <typename TZp>
class NttMul: public ProductFFT<TZp> {
    /*
     *  NTT 与 FFT 两种变换的框架是一样的。但是由于复根的特殊性可以利用，因此这里各实现一下
     *  这里对100万位的乘法实测得：NTT比FFT慢，NTT速度是FFT的1/5
     * */

public:
    void FFT(int k, int n, uint32_t * input, TZp *input1, int input_size, TZp *transform)
    {
        TZp *omega = this->_omega_fft[k];
        uint32_t * pi_k = this->_pi_k[k];
        uint32_t * P = input;
        TZp * P1 = input1;

        if (P1 == NULL) {
            for (int t = 0; t < n-1; t+= 2) {
                uint32_t P_pi_k_t = pi_k[t] < input_size ? P[pi_k[t]] : 0;
                uint32_t P_pi_k_t_1 = pi_k[t+1] < input_size ? P[pi_k[t+1]] : 0;
                transform[t]   = P_pi_k_t;
                transform[t+1] = P_pi_k_t_1;
            }
        } else {
            for (int t = 0; t < n-1; t+= 2) {
                TZp P_pi_k_t = pi_k[t] < input_size ? P1[pi_k[t]] : 0;
                TZp P_pi_k_t_1 = pi_k[t+1] < input_size ? P1[pi_k[t+1]] : 0;

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
                    TZp xPOdd     = omega[m*j] * transform[t+num/2+j];
                    TZp xPOdd1    = omega[m*(j+num/2)] * transform[t+num/2+j];
                    TZp prevTrans = transform[t+j];
                    transform[t+j]       = prevTrans + xPOdd;
                    transform[t+num/2+j] = prevTrans + xPOdd1;
                }
            }
        }
    }

    void rev_FFT(int k, int n, uint32_t * input, TZp *input1, int input_size, TZp *transform)
    {
        FFT(k, n, input, input1, input_size, transform);
        // 必须除以mod P 意义上的n. 这里先求倒数，再作乘积
        TZp n_reciprocal = TZp(1) / TZp(n);
        for (int i = n-1; i > n/2-1; i--) {
            TZp a = transform[n - i] * n_reciprocal;
            transform[n-i] = transform[i] * n_reciprocal;
            transform[i] = a;
        }
        transform[0] *= n_reciprocal;
    }

    void do_fast_ntt(int k, int n2, uint32_t * aa, int aa_length, uint32_t * bb, int bb_length, 
                TZp *transform_aa, TZp *transform_bb, TZp *transform_out)
    {
        FFT(k, n2, aa, NULL, aa_length, transform_aa);
        FFT(k, n2, bb, NULL, bb_length, transform_bb);
        TZp *transform_cc = transform_aa;
        for (int i = 0; i < n2; i++) {
            transform_cc[i] = transform_aa[i] * transform_bb[i];
        }
        rev_FFT(k, n2, NULL, transform_cc, n2, transform_out);
    }
};

class FntMul: public NttMul<Zp0>, NttMul<Zp1>, NttMul<Zp2> {

    // 下面几个数字用于用中国剩余问题恢复overflow的数字
    // // Let P0, P1, P2 = 2113929217, 2013265921, 1811939329
    // // M0 = (P1*P2) * ((P1*P2)^(-1) mod P0)
    // // M1 = (P0*P2) * ((P0*P2)^(-1) mod P1)
    // // M2 = (P0*P1) * ((P0*P1)^(-1) mod P2)
    uint128_t M0; //(0x1d, 0x11e00082ec000093LU); // 0x1d11e00082ec000093
    uint128_t M1; //(0x18eabfd7LU, 0x1b97ff4a91ffff39LU); // 0x18eabfd71b97ff4a91ffff39
    uint128_t M2; //(0xc, 0x75600033e4000036LU); // 0xc75600033e4000036
    // // P012 = P0 * P1 * P2
    uint128_t P012; //(0x18eac000LU, 0xa2d8000162000001LU); // 0x18eac000a2d8000162000001

public:
    FntMul():M0(0x1d, 0x11e00082ec000093LU),
            M1(0x18eabfd7LU, 0x1b97ff4a91ffff39LU),
            M2(0xc, 0x75600033e4000036LU),
            P012(0x18eac000LU, 0xa2d8000162000001LU)
    {
    }
    void init(int max_len)
    {
        NttMul<Zp0>::init(max_len);
        NttMul<Zp1>::init(max_len);
        NttMul<Zp2>::init(max_len);
    }

    int fft_ntt(int k, int n2, uint32_t* aa, int aa_length, uint32_t * bb, int bb_length, 
                        int radix, uint32_t* out, int &out_size)
    {
        Zp0 *buf = new Zp0[5*n2];
        Zp0 *fnt_aa = buf;
        Zp0 *fnt_bb = buf + n2;
        Zp0 *fnt_out0 = buf + n2*2;
        Zp0 *fnt_out1 = buf + n2*3;
        Zp0 *fnt_out2 = buf + n2*4;

        NttMul<Zp0>::do_fast_ntt(k, n2, aa, aa_length, bb, bb_length, (Zp0*)fnt_aa, (Zp0*)fnt_bb, (Zp0*)fnt_out0);
        NttMul<Zp1>::do_fast_ntt(k, n2, aa, aa_length, bb, bb_length, (Zp1*)fnt_aa, (Zp1*)fnt_bb, (Zp1*)fnt_out1);
        NttMul<Zp2>::do_fast_ntt(k, n2, aa, aa_length, bb, bb_length, (Zp2*)fnt_aa, (Zp2*)fnt_bb, (Zp2*)fnt_out2);

        int max_res_len = aa_length + bb_length;

        if (0) { // 只做一个NTT
            int remain = 0;
            for (int i = 0; i < max_res_len; ++i) {
                int o = fnt_out1[i].to_int() + remain;
                if (o >= radix) {
                    out[i] = o % radix;
                    remain = o / radix;
                } else {
                    out[i] = o;
                    remain = 0;
                }
            }
        } else {
            uint128_t remain(0);
            for (int i = 0; i < max_res_len; ++i) {
                uint128_t d = (M0 * fnt_out0[i].n + M1 * fnt_out1[i].n + M2 * fnt_out2[i].n + remain) % P012;

                //uint128_t out_res = d % radix;
                //out[i] = out_res.LOWER;
                //remain = d / radix;

                remain = d.div_mod(radix, out[i]);
            }
        }
        out_size = out[max_res_len-1] != 0 ? max_res_len : max_res_len - 1;
        delete [] buf;
        return 0;
    }

    int fast_prod(uint32_t * aa, int aa_length, uint32_t * bb, int bb_length, 
                  int radix, uint32_t * out, int &out_size) 
    {
        return NttMul<Zp0>::fast_prod(aa, aa_length, bb, bb_length, radix, out, out_size);
    }
};

#endif
