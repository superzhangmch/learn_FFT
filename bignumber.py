import math
import sys
import time
import random
from fft_product import ProductFFT

fft_prod = ProductFFT(100000 * 4)

class BigNumber(object):

    max_digit = 100

    def dec2other(self, num, radix):
        out = []
        while num:
            remain = num % radix
            num = num / radix
            out.append(int(remain))
        return out

    def __init__(self, val, length=None, exp_idx=0, is_neg=False, is_rev=True, radix=10):
        self.radix = radix
        if type(val) == type(1) or type(val) == long:
            if val == 0:
                self.val = []
            else:
                self.val = self.dec2other(abs(val), radix)
            self.length = len(self.val)
            self.exp_idx = exp_idx
            if val == 0:
                self.is_neg = True
            else:
                self.is_neg = True if val / abs(val) == -1 else False
            return

        if type(val) == type(""):
            val = [int(i) for i in val]

        if length is None:
            length = len(val)
        if is_rev:
            self.set_val(val, length, exp_idx, is_neg)
        else:
            self.set_normal_val(val, length, exp_idx, is_neg)

    def set_val(self, val, length, exp_idx, is_neg=False):
        assert len(val) >= length
        self.val = val
        self.length = length
        self.exp_idx = exp_idx
        self.is_neg = is_neg

        max_len = length
        for i in xrange(max_len-1, -1, -1):
            if self.val[i] == 0:
                self.length -= 1
            else:
                break

    def extend_len(self, size):
        if len(self.val) < size:
            self.val += [0] * (size - len(self.val))

    def set_normal_val(self, val, length, exp_idx=0, is_neg=False):
        self.length = length
        aa = val + []
        aa.reverse()
        self.val = aa
        self.exp_idx = exp_idx
        self.is_neg = is_neg

        max_len = length
        for i in xrange(max_len-1, -1, -1):
            if self.val[i] == 0:
                self.length -= 1
            else:
                break

    def get_normal_val(self):
        aa = self.val + []
        aa = aa[:self.length]
        aa.reverse()
        return aa

    def get_all_numbers(self):
        return "".join([str(i) for i in self.get_normal_val()])

    def get_string_val(self):
        if self.length == 0:
            return "0"
        c = int(round(math.log(self.radix) / math.log(10)))
        if self.radix == 10:
            res = "".join([str(i) for i in self.get_normal_val()])
        else:
            prefix = "0" * c
            res = "".join([(prefix + str(i))[-c:] for i in self.get_normal_val()])
        if len(res)> 600:
            res = res[:500] + "..." + res[-100:]
        res = "0."+res + "E" + str(c*(self.exp_idx + self.length)) + "("+str(self.length*c)+")"
        res = ("-" if self.is_neg else "") + res
        return res

    def __str__(self):
        return self.get_string_val()

    def adjust_precision(self, precision):
        if self.length <= precision:
            return
        else:
            self.exp_idx += (self.length - precision)
            self.length = precision
            self.val = self.val[-precision:]

    def adjust_exp_idx(self, new_exp_idx):
        if new_exp_idx == self.exp_idx:
            return
        if new_exp_idx < self.exp_idx:
            adjust_val = self.exp_idx - new_exp_idx
            self.val = ([0] * adjust_val) + self.val
            self.exp_idx = new_exp_idx
            self.length += adjust_val
            return
        else:
            adjust_val = new_exp_idx - self.exp_idx
            assert adjust_val <= self.length
            all_z = True
            for i in xrange(adjust_val):
                if self.val[i] != 0:
                    all_z = False
                    break
            assert all_z == True
            self.length -= adjust_val
            self.val = self.val[adjust_val:]
            self.exp_idx = new_exp_idx
        
    # =====================
    def __mul__(self, num1):
        if self.length:
            assert self.val[self.length - 1] != 0
        if num1.length:
            assert num1.val[num1.length - 1] != 0
        res = BigNumber(0)
        fft_prod.fast_prod(self, num1, res, self.radix, True)
        if res.length:
            assert res.val[res.length - 1] != 0
        return res

    def __add__(self, num1):
        if self.length:
            assert self.val[self.length - 1] != 0
        if num1.length:
            assert num1.val[num1.length - 1] != 0

        if self.is_neg != num1.is_neg:
            # sub: (+) + (-) || (-) + (+)
            res = self.do_sub(num1)
        # add: (+) + (+) || (-) + (-)
        else:
            res = self.do_add(num1)
        if res.length:
            assert res.val[res.length - 1] != 0
        return res

    def __sub__(self, num1):
        if self.length:
            assert self.val[self.length - 1] != 0
        if num1.length:
            assert num1.val[num1.length - 1] != 0

        if self.is_neg != num1.is_neg:
            # add: (+) - (-) => (+) || (-) - (+) => (-)
            res = self.do_add(num1)
        else: # sub: (-) - (-) || (+)-(+)
            res = self.do_sub(num1)
        if res.length:
            assert res.val[res.length - 1] != 0
        return res

    def do_add(self, num1):
        """ same sign, result.sign = self.sign """
        if self.length == 0:
            if self.is_neg == num1.is_neg:
                return num1
            else:
                val = num1.val+[]
                return BigNumber(val, num1.length, num1.exp_idx, is_neg=self.is_neg)
        if num1.length == 0:
            return self

        if self.exp_idx <= num1.exp_idx:
            a, b = self, num1
        else:
            a, b = num1, self
        # now a.exp_idx < b.exp_idx. b should be padded with zero at the begining

        exp_idx_diff = b.exp_idx - a.exp_idx
        if a.length >= b.length + exp_idx_diff:
            a_is_longer = True
        else:
            a_is_longer = False

        use_exp_idx = a.exp_idx
        min_length = min(a.length, b.length + exp_idx_diff)
        res_length = max(a.length, b.length + exp_idx_diff)
        radix = self.radix

        sum_arr = [0] * (1 + res_length)
        length_a = a.length
        a = a.val
        b = b.val
        #print "xx", exp_idx_diff, min_length, res_length
        #print "xx, a=", a
        #print "xx, b=", b
        remain = 0
        for i in xrange(res_length):
            if i < exp_idx_diff:
                if i < length_a:
                    sum_arr[i] = a[i]
            elif i < min_length:
                s = a[i] + b[i - exp_idx_diff] + remain
                if s >= radix:
                    sum_arr[i] = s % radix
                    remain = s / radix
                else:
                    sum_arr[i] = s
                    remain = 0
                #print "s=%d, remain=%d, radix=%d, res=%d" %( s, remain, radix, sum_arr[i])
            else:
                if a_is_longer:
                    s = a[i]
                else:
                    s = b[i - exp_idx_diff]
                if remain:
                    s += remain
                    sum_arr[i] = s % radix
                    remain = s / radix
                else:
                    sum_arr[i] = s
        if remain:
            sum_arr[res_length] = remain
            res_length += 1
        return BigNumber(sum_arr, res_length, use_exp_idx, is_neg=self.is_neg)

    def do_sub(self, num1):
        """ self.sign != num1.sign, result.sign = self.sign if (|self| > |num1| else (!self.sign)) """
        if self.length == 0:
            res_neg = not self.is_neg
            if res_neg == num1.is_neg:
                return num1
            else:
                val = num1.val+[]
                return BigNumber(val, num1.length, num1.exp_idx, is_neg=not self.is_neg)
        if num1.length == 0:
            return self

        is_a_self = False
        if self.exp_idx <= num1.exp_idx:
            a, b = self, num1
            is_a_self = True
        else:
            a, b = num1, self
        # now a.exp_idx < b.exp_idx. b should be padded with zero at the begining

        exp_idx_diff = b.exp_idx - a.exp_idx
        is_a_sub_b = False
        if a.length > b.length + exp_idx_diff:
            is_a_sub_b = True
        elif a.length < b.length + exp_idx_diff:
            is_a_sub_b = False
        else:
            all_equal = True
            for i in xrange(a.length - 1, exp_idx_diff - 1, -1):
                if a.val[i] > b.val[i - exp_idx_diff]:
                    is_a_sub_b = True
                    all_equal = False
                    break
                elif a.val[i] < b.val[i - exp_idx_diff]:
                    all_equal = False
                    break
            if all_equal:
                is_a_sub_b = True

        use_exp_idx = a.exp_idx
        min_length = min(a.length, b.length + exp_idx_diff)
        res_length = max(a.length, b.length + exp_idx_diff)
        radix = self.radix

        sub_arr = [0] * res_length
        length_a = a.length
        a = a.val
        b = b.val
        #print "xx", exp_idx_diff, min_length, res_length
        #print "xx, a=", a
        #print "xx, b=", b
        #print "is_a_sub_b", is_a_sub_b
        borrow = 0
        for i in xrange(res_length):
            if i < exp_idx_diff:
                if is_a_sub_b:
                    sub_arr[i] = a[i]
                else:
                    if i < length_a:
                        s = - a[i] - borrow
                        if s < 0:
                            s += radix
                            borrow = 1
                        else:
                            borrow = 0
                        sub_arr[i] = s
                    else:
                        s = - 0 - borrow
                        if s < 0:
                            s += radix
                            borrow = 1
                        else:
                            borrow = 0
                        sub_arr[i] = s
            elif i < min_length:
                if is_a_sub_b:
                    s = a[i] - b[i - exp_idx_diff] - borrow
                else:
                    s = b[i - exp_idx_diff] - a[i] - borrow
                if s < 0:
                    sub_arr[i] = s + radix
                    borrow = 1
                else:
                    sub_arr[i] = s
                    borrow = 0
                #print "s=%d, remain=%d, radix=%d, res=%d" %( s, remain, radix, sum_arr[i])
            else:
                if is_a_sub_b:
                    s = a[i] - borrow
                else:
                    s = b[i - exp_idx_diff] - borrow
                if s < 0:
                    sub_arr[i] = s + radix
                    borrow = 1
                else:
                    sub_arr[i] = s
                    borrow = 0

        for i in xrange(res_length - 1, -1, -1):
            if sub_arr[i] == 0:
                res_length -= 1
            else:
                break

        is_self_bigger = None
        if is_a_self:
            is_self_bigger = True if is_a_sub_b else False
        else:
            is_self_bigger = False if is_a_sub_b else True
        if is_self_bigger:
            res_neg = self.is_neg
        else:
            res_neg = not self.is_neg
        return BigNumber(sub_arr, res_length, use_exp_idx, is_neg=res_neg)

    def __div__(self, num):
        if self.length:
            assert self.val[self.length - 1] != 0
        if num.length:
            assert num.val[num.length - 1] != 0
        res = self * num.reciprocal()
        if res.length:
            assert res.val[res.length - 1] != 0
        return res

    def reciprocal(self):
        b = self
        assert b.length != 0
        if b.length == 1:
            b.adjust_exp_idx(b.exp_idx - 2)

        # print b.val, b.length, 11111111111111
        if b.val[b.length - 1] == 1:
            val = [1]
        else:
            val = [b.val[b.length - 1]]
        init_expidx = -b.length - b.exp_idx - len(val)
        
        def f(x):
            return x*(BigNumber(2) - b*x)
        aa = BigNumber(val, 1, init_expidx)
        return self.Newton_method(f, aa)

    def sqrt(self):
        b = self
        if b.length == 0:
            return BigNumber([], 0, 0)
        assert b.is_neg == False

        # print b.val, b.length, 11111111111111222222222222222
        #----
        first_num = b.val[b.length - 1]
        second_num = 0 if b.length == 1 else b.val[b.length - 2]
        if (b.length + b.exp_idx) % 2 == 0:
            init_expidx = -(b.length + b.exp_idx - 2) / 2
            val = first_num * self.radix + second_num
        else:
            init_expidx = -(b.length + b.exp_idx - 1) / 2
            val = first_num
        if val != 1:
            val = int(math.floor(self.radix * 1/math.sqrt(val)))
            init_expidx -= 1
        val = [val]

        #----
        def f(x):
            ret = x + x * (BigNumber(1) - b*x*x) * BigNumber([self.radix/2], 1, -1)
            return ret
        #----
        aa = BigNumber(val, 1, init_expidx)
        # print aa, val, init_expidx
        #print aa, aa.length, aa.exp_idx, aa.val
        x = self.Newton_method(f, aa)
        return x.reciprocal()

    def Newton_method(self, func, init_val, max_loop=100):
        f = func
        x = init_val
        last_2 = [-1, -1]
        has_first_2 = False
        digit_num = 0
        # print "init=", init_val.length, init_val.val, init_val.exp_idx
        for i in xrange(max_loop):
            x = f(x)
            x_len = x.length
            if not has_first_2:
                if x_len >= 2:
                    if last_2[0] == x.val[x_len-1] and last_2[1] == x.val[x_len-2]:
                        has_first_2 = True
                        x.val = x.val[x_len-2:x_len]
                        x.exp_idx += (x_len - 2)
                        x.length = 2
                        digit_num = 2
                    else:
                        last_2[0], last_2[1] = x.val[x_len-1], x.val[x_len-2]
                        if x.length > 4: 
                            x.val = x.val[x_len-4:x_len]
                            x.exp_idx += (x_len - 4)
                            x.length = 4
            else:
                digit_num *= 2
        
                x.val = x.val[x_len-digit_num:x_len]
                x.exp_idx += (x_len - digit_num)
                x.length = digit_num
                
            #print "round=%d precise=%d exp_idx=%d length=%d num=%s" % (i, digit_num, x.exp_idx, x.length, str(x)[:200])
            if digit_num >= self.max_digit:
                # print "vvvvvv", digit_num, self.max_digit, last_2
                break
        return x
