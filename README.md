# learn_FFT

- FFT/FNT(快速傅里叶变换、快速数论变换) 做大整数的乘法，以及基于此的除法，开方。
  - 为了试探FFT/FNT到底能多快，发现在同一机器上，100 万位乘法FFT(基数=10)大约1.0秒, FNT(基数=10^9)大约1.3秒。相形之下，apfloat可以做到0.3秒。
  - FFT,FNT 可以多线程加速的。不过以上测试都没多线程加速。看apfloat在原生FNT基础上做了好多优化，所以那么快。
  
- python 和 C/C++分别实现了一遍。先用python 验证，然后用C/C++加速
  - c/c++实现封装成了so,仍然是被python调用使用
- 对于非c/c++部分，用pypy执行可以加速
