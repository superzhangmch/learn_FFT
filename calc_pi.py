#encoding: utf8
'''
测试验证。用pypy来加速运行。否则太慢了
'''
import math
import sys
import time
import random
from bignum_c import BigNum
#from bignum_py import BigNumber as BigNum
from pi import pi, match_num

time_start = time.time()

def calc_e():
    '''
    自然对数: e= sum(1/n!, n = 1 .. N)
    注意:公式中除以的数字都比较小，应该直接做这种除法，不应该用到FFT的乘除法。
    '''
    max_digit = 10000
    radix = 1000000000
    BigNum.init(max_digit, radix, True)

    sum1 = BigNum(0)
    x = BigNum(1)
    for i in xrange(10000):
        sum1 = sum1 + x
        x = x / BigNum(i + 1)
        if i % 1000 == 0:
            print i
    print sum1

def calc_pi_use_machin():
    '''
    Machin公式计算PI.
    Pi = 16*sum(1/(2*i+1)/5^(2*i+1)*(-1)^i)-4*sum(1/(2*i+1)/239^(2*i+1)*(-1)^i) 
    注意:公式中除以的数字都比较小，应该直接做这种除法，不应该用到FFT的乘除法。
    '''
    calc_digits = 5000
    radix = 1000000000

    max_digit = calc_digits / int(math.log(radix) / math.log(10)) + 20
    print "loop_cnt", max_digit
    BigNum.init(max_digit, radix, False)

    x = BigNum(5)
    y = BigNum(239)

    sum1 = BigNum(0)
    sum2 = BigNum(0)

    num5 = BigNum(5)
    num239 = BigNum(239)
    num_1 = BigNum([1], 1, 0, False)
    num_neg1 = BigNum([1], 1, 0, True)

    for i in xrange(calc_digits):
        tm = time.time()
        a = 2*i+1
        x = x / num5 / num5
        y = y / num239 / num239
        xx = BigNum(a, is_neg=True if i % 2 == 1 else False)
        sum1 = sum1 + x / xx
        sum2 = sum2 + y / xx
        if i % 1000 == 0:
            print "%d -> %.4f" % (i, time.time() - tm)
    PI = BigNum(16) * sum1 - BigNum(4) * sum2
    aa = PI.get_all_numbers().strip("0")
    print aa
    match_num(pi, aa)
    print "pi=", PI

def calc_pi():
    '''
    AGM 算术几何平均法
    '''
    def radix_len(radix):
        return int(round(math.log(radix)/ math.log(10)))
    radix10_digit = 1000000

    # FFT, 10进制下没法算到100万位，有溢出。更大进制更不行了
    # FNT 则无此限制
    radix = 1000000000
    use_fft = False if radix > 10 else True
    max_digit = radix10_digit / radix_len(radix)
    BigNum.init(max_digit, radix, use_fft)

    a = BigNum(1)
    b = (BigNum(1) / BigNum(2)).sqrt()
    t = BigNum(1) / BigNum(4)
    p = BigNum(1)

    round_cnt = int(math.ceil(math.log(radix10_digit) / math.log(2)))

    print "init ok, now to calc %d rounds:" % (round_cnt)

    for i in xrange(round_cnt):
        print "begin %d" % i
        tm = time.time()
        a1 = (a + b) / BigNum(2)
        b = (a*b).sqrt()
        diff = a - a1
        #t = t - p * diff * diff
        # 如果多线程计算乘法，则自乘也未必节约多少时间
        t = t - p * diff.pow2()
        p = BigNum(2)*p
        a = a1
        print "%d -> %.4f" % (i, time.time() - tm)

    #PI = (a + b) * (a + b) / BigNum(4) / t
    PI = (a + b).pow2() / BigNum(4) / t
    aa = PI.get_all_numbers().lstrip("0")
    match_num(pi, aa)
    print PI
calc_pi()
#calc_pi_use_machin()
#calc_e()
print time.time() - time_start
sys.exit(0)
