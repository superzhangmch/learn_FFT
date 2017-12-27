g++ fft_ntt_test.cpp -O2 fft_common.cpp 

g++ -O2 --shared -fPIC -o fft_ntt.so fft_ntt.cpp fft_common.cpp
