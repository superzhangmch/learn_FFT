#encoding:utf8
import math
import sys

class ProductFFT(object):

    def __init__(self, max_len):
        print "init begin"
        k = int(math.ceil(math.log(max_len) / math.log(2)))
        self._max_result_len = 2 * 2 **k
        self._max_k = k + 1
        self._omega = {}
        self._pi_k = {}
        for i in xrange(1, self._max_k + 1, 1):
            self._omega[i], self._pi_k[i] = self.gen_omega_pi_k(i)
        print "init ok"

    def primitive_root(self, n):
        """ calculate primitive root for n """
        return math.cos(2*math.pi / n) + 1j * math.sin(2*math.pi / n)

    def gen_omega_pi_k(self, k):
        n = 2 ** k
        w = self.primitive_root(n)
        w_ret = [w ** i for i in xrange(n)]
        def rev_bin(t, k):
            b = bin(t)[2:]
            b1 = ("0"*(k - len(b))) + b
            b2 = b1[::-1]
            t1 = int("0b"+b2, 2)
            #print t, b, b1, b2, t1
            return t1
        pi_k_ret = [rev_bin(t, k) for t in xrange(n)]
        return w_ret, pi_k_ret

    def FFT(self, omega, pi_k, P, k):
        """ O(n*log(n)). from <<computer algorithms: introduction to design and analysis>>"""
        transform = [0.] * len(P)
        n = 2 ** k
        for t in xrange(0, n-1, 2):
            transform[t]   = P[pi_k[t]] + P[pi_k[t+1]]
            transform[t+1] = P[pi_k[t]] - P[pi_k[t+1]]
        m = n/2
        num = 2
        for d in xrange(k-2, -1, -1):
            m /= 2
            num *= 2
            for t in xrange(0, (2**d-1)*num + 1, num):
                for j in xrange(num/2):
                    xPOdd = omega[m*j] * transform[t+num/2+j]
                    prevTrans= transform[t+j]
                    transform[t+j] = prevTrans+xPOdd
                    transform[t+num/2+j] = prevTrans-xPOdd
        return transform

    def rev_FFT(self, omega, pi_k, P, k):
        """ from <<computer algorithms: introduction to design and analysis>>"""
        n = 2 ** k
        transform = self.FFT(omega, pi_k, P, k)
        for i in xrange(n-1, n/2 - 1, -1):
            transform[i], transform[n-i] = transform[n - i] / n, transform[i] / n
        transform[0] /= n
        return transform

    def slow_FFT(self, w_arr, P, r=False):
        """ O(n*n). FFT/rev_FFT就是个矩阵运算，所以可以按矩阵计算的方式进行
        另外numpy.fft.fft/numpy.fft.ifft 是numpy的FFT实现
        """
        n = len(w_arr)
        transform = [0.] * n
    
        if r == True:
            for i in xrange(n):
                for j in xrange(n):
                    transform[i] += w_arr[-(i*j) % n] * P[j]
                transform[i] *= 1. / n
        else:
            for i in xrange(n):
                for j in xrange(n):
                    transform[i] += w_arr[(i*j) % n] * P[j]
        return transform

    def fast_prod(self, aa, bb, result, radix, fast=True):
        # aa or bb == 0
        if aa.length == 0 or bb.length == 0:
            return result.set_val([], 0, 0)
        # aa or bb is small
        if aa.length == 1 or bb.length == 1:
            if aa.length == 1:
                small_val, big = aa.val[0], bb
            else:
                small_val, big = bb.val[0], aa
            out = [0] * (big.length + 1)
            out_len = big.length
            remain = 0
            for i in xrange(big.length):
                res = big.val[i] * small_val + remain
                #print i, '->', big.val[i], small_val, remain, res
                if res >= radix:
                    out[i] = res % radix
                    remain = res / radix
                else:
                    out[i] = res
                    remain = 0
            #print ""
            if remain:
                #print 'xxxxxxxxxxxxxxxxxxxxx', small_val
                out[big.length] = remain
                out_len = big.length + 1
            is_neg = False if aa.is_neg == bb.is_neg else True
            result.set_val(out, out_len, aa.exp_idx + bb.exp_idx, is_neg=is_neg)
            #print 'outoutout', out, out_len
            assert out[out_len-1] != 0
            return 

        # aa or bb both big
        def get_len(len_aa, len_bb):
            k =  int(math.ceil(math.log(max(len_aa, len_bb)) / math.log(2)))
            return 2 * 2 **k, k + 1
        n2, k = get_len(aa.length, bb.length)

        aa.extend_len(n2)
        bb.extend_len(n2)

        aa1 = aa.val
        bb1 = bb.val
    
        w_arr, pi_k_arr = self._omega[k], self._pi_k[k]

        # === FFT
        if not fast:
            aa2 = self.slow_FFT(w_arr, aa1)
        else:
            aa2 = self.FFT(w_arr, pi_k_arr, aa1, k)

        # === FFT
        if not fast:
            bb2 = self.slow_FFT(w_arr, bb1)
        else:
            bb2 = self.FFT(w_arr, pi_k_arr, bb1, k)

        # ==== rev FFT
        cc = [0.] * n2
        for i in xrange(n2):
            cc[i] = aa2[i]*bb2[i]
        if not fast:
            out = self.slow_FFT(w_arr, cc, True)
        else:
            out = self.rev_FFT(w_arr, pi_k_arr, cc, k)

        # ==== 
        max_res_len = aa.length + bb.length
        out = [int(round(o.real)) for o in out]
        out1 = [0] * n2
        remain = 0
        zero_cnt = 0
        not_count = 0
        for i in xrange(len(out)):
            o = out[i]
            o += remain
            if o >= radix:
                out1[i] = o % radix
                remain = o / radix
            else:
                out1[i] = o
                remain = 0

            if (not not_count) and out1[i] == 0:
                zero_cnt += 1
            if out1[i] != 0:
                not_count = 1

            if i >= max_res_len:
                break
        res_len = max_res_len if out1[max_res_len - 1] != 0 else max_res_len - 1
        is_neg = False if aa.is_neg == bb.is_neg else True
        result.set_val(out1, res_len, aa.exp_idx + bb.exp_idx, is_neg=is_neg)
        assert out1[res_len-1] != 0
        #if not_count == 0:
        #    result.set_val(out1, res_len, aa.exp_idx + bb.exp_idx, is_neg=is_neg)
        #else:
        #    result.set_val(out1[zero_cnt:], res_len-zero_cnt, aa.exp_idx + bb.exp_idx+zero_cnt, is_neg=is_neg)
 
