import math
import sys
import time
import random
from bignumber import BigNumber, fft_prod
from pi import pi, match_num

time_start = time.time()

def calc_e():
    BigNumber.radix = 1000
    BigNumber.max_digit = 10000

    sum1 = BigNumber(0)
    x = BigNumber(1)
    for i in xrange(1000):
        sum1 = sum1 + x
        x = x / BigNumber(i + 1)
        print i
    print sum1

def calc_pi_use_machin():
    calc_digits = 1000

    BigNumber.radix = 1000
    BigNumber.max_digit = calc_digits / int(math.log(BigNumber.radix) / math.log(10)) + 20

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
    BigNumber.max_digit = 1000
    a = BigNumber(1)
    b = BigNumber(1) / BigNumber(2).sqrt()
    t = BigNumber(1) / BigNumber(4)
    p = BigNumber(1)
    
    print "init ok"
    
    num = 1000
    
    for i in xrange(10):
        #num *= 2
        #BigNumber.max_digit = num
        a1 = (a + b) / BigNumber(2)
        b = (a*b).sqrt()
        diff = a - a1
        t = t - p * diff * diff
        p = BigNumber(2)*p
        a = a1
        print "-------------", i

        PI = (a + b) * (a +b) / BigNumber(4) / t
        aa = PI.get_all_numbers()
        match_num(pi, aa)
        print PI
#calc_pi()
#calc_pi_use_machin()
calc_e()
print time.time() - time_start
sys.exit(0)
