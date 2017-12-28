#encoding: utf8
'''
测试验证。用pypy来加速运行。否则太慢了
'''
import math
import sys
import time
import random
from bignumber import BigNumber
from pi import pi, match_num

time_start = time.time()

def calc_e():
    '''
    自然对数: e= sum(1/n!, n = 1 .. N)
    注意:公式中除以的数字都比较小，应该直接做这种除法，不应该用到FFT的乘除法。
    '''
    max_digit = 10000
    radix = 1000
    BigNumber.init(max_digit, radix, True)

    sum1 = BigNumber(0)
    x = BigNumber(1)
    for i in xrange(1000):
        sum1 = sum1 + x
        x = x / BigNumber(i + 1)
        print i
    print sum1

def calc_pi_use_machin():
    '''
    Machin公式计算PI.
    Pi = 16*sum(1/(2*i+1)/5^(2*i+1)*(-1)^i)-4*sum(1/(2*i+1)/239^(2*i+1)*(-1)^i) 
    注意:公式中除以的数字都比较小，应该直接做这种除法，不应该用到FFT的乘除法。
    '''
    calc_digits = 1000

    max_digit = calc_digits / int(math.log(BigNumber.radix) / math.log(10)) + 20
    radix = 1000
    BigNumber.init(max_digit, radix, True)

    x = BigNumber(5)
    y = BigNumber(239)

    sum1 = BigNumber(0)
    sum2 = BigNumber(0)
    for i in xrange(calc_digits):
        a = 2*i+1
        is_neg = True if i % 2 == 1 else False
        xx = BigNumber([1], 1, 0, is_neg) / BigNumber(a)
        x = x / BigNumber(5) / BigNumber(5)
        y = y / BigNumber(239) / BigNumber(239)
        sum1 = sum1 + xx * x
        sum2 = sum2 + xx * y
        if i % 10 == 0:
            print i
    PI = BigNumber(16) * sum1 - BigNumber(4) * sum2
    aa = PI.get_all_numbers().strip("0")
    print aa
    match_num(pi, aa)
    print "pi=", PI

def calc_pi():
    max_digit = 4000
    radix = 1000000000
    radix = 10
    BigNumber.init(max_digit, radix, True)

    a = BigNumber(1)
    b = BigNumber(1) / BigNumber(2).sqrt()
    t = BigNumber(1) / BigNumber(4)
    p = BigNumber(1)
    
    print "init ok"
    
    for i in xrange(14):
        a1 = (a + b) / BigNumber(2)
        b = (a*b).sqrt()
        diff = a - a1
        t = t - p * diff * diff
        p = BigNumber(2)*p
        a = a1
        print "-------------", i

        PI = (a + b) * (a +b) / BigNumber(4) / t
        aa = PI.get_all_numbers().lstrip("0")
        match_num(pi, aa)
        print PI
calc_pi()
#calc_pi_use_machin()
#calc_e()
print time.time() - time_start
sys.exit(0)
