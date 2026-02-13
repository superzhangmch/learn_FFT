#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <iomanip>
#include "bigfloat.h"

// ---- Helpers templated on any BigFloat variant ----

template<typename BF>
BF sqrt_bf(const BF& S) {
    constexpr int TL = BF::BITS / BF::TBITS - 1; // top limb index
    using T2 = typename DoubleWidth<typename BF::value_type>::type;
    double top = (double)(uint64_t)(T2)S.mantissa_[TL];
    double s_approx = std::ldexp(top / std::ldexp(1.0, BF::TBITS - 1), S.exp_ + BF::TOP_BIT);
    BF x(std::sqrt(s_approx));
    BF half(0.5);
    for (int i = 0; i < 12; i++)
        x = (x + S / x) * half;
    return x;
}

template<typename BF>
int floor_int(const BF& val) {
    if (val.is_zero()) return 0;
    if (val.exp_ >= 0) return -1;
    int shift = -val.exp_;
    if (shift > BF::TOP_BIT) return 0;
    constexpr int NL = BF::BITS / BF::TBITS;
    using T2 = typename DoubleWidth<typename BF::value_type>::type;
    uint64_t result = 0;
    int lo = shift / BF::TBITS;
    int bo = shift % BF::TBITS;
    int rpos = 0;
    for (int k = 0; rpos < 31 && lo + k < NL; k++) {
        uint64_t v = (uint64_t)(T2)val.mantissa_[lo + k];
        if (k == 0) {
            v >>= bo;
            result |= v << rpos;
            rpos += BF::TBITS - bo;
        } else {
            result |= v << rpos;
            rpos += BF::TBITS;
        }
    }
    return (int)(result & 0x7FFFFFFF);
}

template<typename BF>
std::string to_digits(BF val, int ndigits) {
    std::string result;
    if (val.is_negative()) { result += '-'; val = -val; }
    int ipart = floor_int(val);
    result += std::to_string(ipart);
    result += '.';
    BF frac = val - BF((double)ipart);
    BF ten(10.0);
    for (int i = 0; i < ndigits; i++) {
        frac = frac * ten;
        int d = floor_int(frac);
        if (d < 0 || d > 9) d = 0;
        result += ('0' + d);
        frac = frac - BF((double)d);
    }
    return result;
}

template<typename BF>
std::string compute_pi(int ndigits, int nterms) {
    BF c640320(640320.0);
    BF C = c640320 * c640320 * c640320;
    BF sum(13591409.0);
    BF term(1.0);

    for (int k = 1; k <= nterms; k++) {
        long long k6 = 6LL * k, k3 = 3LL * k;
        BF num(1.0);
        for (int j = 0; j < 6; j++)
            num = num * BF((double)(k6 - j));
        BF den = BF((double)k3) * BF((double)(k3-1)) * BF((double)(k3-2));
        den = den * BF((double)k) * BF((double)k) * BF((double)k) * C;
        term = term * num / den;
        term = -term;
        sum = sum + term * BF(13591409.0 + 545140134.0 * k);
    }

    BF sqrt10005 = sqrt_bf(BF(10005.0));
    BF pi = BF(426880.0) * sqrt10005 / sum;
    return to_digits(pi, ndigits);
}

// Reference pi
static const char* PI_REF =
    "3."
    "14159265358979323846264338327950288419716939937510"
    "58209749445923078164062862089986280348253421170679"
    "82148086513282306647093844609550582231725359408128"
    "48111745028410270193852110555964462294895493038196"
    "44288109756659334461284756482337867831652712019091"
    "45648566923460348610454326648213393607260249141273";

int count_match(const std::string& a, const char* b) {
    int m = 0;
    for (size_t i = 0; i < a.size() && b[i]; i++) {
        if (a[i] == b[i]) m++; else break;
    }
    return m - 2; // subtract "3."
}

template<typename BF>
void run_test(const char* label, int ndigits, int nterms) {
    std::cout << label << " (" << BF::BITS << "-bit, "
              << BF::TBITS << "-bit limb x " << (BF::BITS / BF::TBITS) << ")" << std::endl;

    auto t0 = std::chrono::high_resolution_clock::now();
    std::string pi = compute_pi<BF>(ndigits, nterms);
    auto t1 = std::chrono::high_resolution_clock::now();

    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    int matched = count_match(pi, PI_REF);

    std::cout << "  pi = " << pi.substr(0, 52) << "..." << std::endl;
    std::cout << "  Correct digits: " << matched << " / " << ndigits
              << "  Time: " << std::fixed << std::setprecision(1) << ms << " ms" << std::endl << std::endl;
}

int main() {
    std::cout << "=== Computing pi with different base types ===" << std::endl << std::endl;

    run_test<Float1024>     ("uint64_t", 300, 24);
    run_test<Float1024_32>  ("uint32_t", 300, 24);
    run_test<Float1024_8>   ("uint8_t ", 300, 24);
    run_test<Float1024_4>   ("uint4_t ", 300, 24);
    run_test<Float1024_2>   ("uint2_t ", 300, 24);

    return 0;
}
