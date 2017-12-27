from ctypes import cdll, c_int, byref, c_float
import time
import math
import random

libc = cdll.LoadLibrary("./fft_ntt.so")

def init_fft(max_digit_num):
    libc.fft_init(max_digit_num)

def init_fnt(max_digit_num):
    libc.fnt_init(max_digit_num)

def do_fnt(aa, bb, radix):
    aa_s = len(aa)
    bb_s = len(bb)
    k =  int(math.ceil(math.log(max(aa_s, bb_s)) / math.log(2)))
    cc_s = 2 **(k+1)
    AA = (c_int * aa_s)(*aa)
    BB = (c_int * bb_s)(*bb)
    CC = (c_int * cc_s)()
    cc_s = libc.do_fnt(byref(AA), c_int(aa_s), byref(BB), c_int(bb_s), radix, byref(CC))
    return list(CC)[:cc_s]

def do_fft(aa, bb, radix):
    aa_s = len(aa)
    bb_s = len(bb)
    k =  int(math.ceil(math.log(max(aa_s, bb_s)) / math.log(2)))
    cc_s = 2 **(k+1)
    AA = (c_int * aa_s)(*aa)
    BB = (c_int * bb_s)(*bb)
    CC = (c_int * cc_s)()
    cc_s = libc.do_fft(byref(AA), c_int(aa_s), byref(BB), c_int(bb_s), radix, byref(CC))
    print cc_s
    return list(CC)[:cc_s]

if __name__ == "__main__":
    aa_s = 1000000
    bb_s = 1000000
    radix = 1000000000
    aa = []
    bb = []
    for i in xrange(aa_s):
        aa.append(random.randint(0, radix - 1))
    for i in xrange(bb_s):
        bb.append(random.randint(0, radix - 1))
    tm = time.time()
    init_fnt(1000000)
    cc = do_fnt(aa, bb, radix)
    print len(cc)
    print ((sum(aa) % 9) * (sum(bb) % 9)) % 9
    print sum(cc) % 9
    print time.time() - tm
