#include <iostream>
#include <string>
#include <cmath>
#include "float512.h"

// Square root via Newton's method: x_{n+1} = (x_n + S/x_n) / 2
Float512 sqrt512(const Float512& S) {
    double s_approx = std::ldexp(
        (double)S.mantissa_[7] / (double)(uint64_t(1) << 63),
        S.exp_ + 511);
    Float512 x(std::sqrt(s_approx));
    Float512 half(0.5);

    // Each iteration doubles precision: 53 -> 106 -> 212 -> 424 -> 848
    for (int i = 0; i < 10; i++)
        x = (x + S / x) * half;
    return x;
}

// Extract floor(val) for non-negative Float512, val must be < 2^31
int floor_int(const Float512& val) {
    if (val.is_zero()) return 0;
    if (val.exp_ >= 0) return -1; // too large
    int shift = -val.exp_;
    if (shift > 511) return 0;

    int limb = shift / 64;
    int bit  = shift % 64;
    uint64_t result = 0;
    if (limb < 8) result = val.mantissa_[limb] >> bit;
    if (bit > 0 && limb + 1 < 8) result |= val.mantissa_[limb + 1] << (64 - bit);
    return (int)(result & 0x7FFFFFFF);
}

// Extract decimal digits from Float512
std::string to_digits(Float512 val, int ndigits) {
    std::string result;
    if (val.is_negative()) { result += '-'; val = -val; }

    int ipart = floor_int(val);
    result += std::to_string(ipart);
    result += '.';

    Float512 frac = val - Float512((double)ipart);
    Float512 ten(10.0);

    for (int i = 0; i < ndigits; i++) {
        frac = frac * ten;
        int d = floor_int(frac);
        if (d < 0 || d > 9) d = 0;
        result += ('0' + d);
        frac = frac - Float512((double)d);
    }
    return result;
}

int main() {
    std::cout << "Chudnovsky algorithm with Float512 (479-bit mantissa)" << std::endl;
    std::cout << "Computing pi..." << std::endl << std::endl;

    // pi = 426880 * sqrt(10005) / sum
    // sum = sum_{k=0}^{N} T(k) * L(k)
    // L(k) = 13591409 + 545140134 * k
    // T(0) = 1
    // T(k) = T(k-1) * [-(6k)(6k-1)(6k-2)(6k-3)(6k-4)(6k-5)]
    //                 / [(3k)(3k-1)(3k-2) * k^3 * 640320^3]

    Float512 C = Float512(Int512("262537412640768000")); // 640320^3, exact

    Float512 sum(13591409.0); // k=0 term
    Float512 term(1.0);       // T(0) = 1

    for (int k = 1; k <= 14; k++) {
        long long k6 = 6LL * k;
        long long k3 = 3LL * k;

        // numerator: (6k)(6k-1)(6k-2)(6k-3)(6k-4)(6k-5)
        Float512 num(1.0);
        for (int j = 0; j < 6; j++)
            num = num * Float512((double)(k6 - j));

        // denominator: (3k)(3k-1)(3k-2) * k^3 * C
        Float512 den = Float512((double)k3)
                     * Float512((double)(k3 - 1))
                     * Float512((double)(k3 - 2));
        den = den * Float512((double)k) * Float512((double)k) * Float512((double)k);
        den = den * C;

        term = term * num / den;
        term = -term; // (-1)^k sign flip

        Float512 L(13591409.0 + 545140134.0 * k);
        sum = sum + term * L;
    }

    Float512 sqrt10005 = sqrt512(Float512(10005.0));
    Float512 pi = Float512(426880.0) * sqrt10005 / sum;

    // Display result
    int ndigits = 150;
    std::string pi_str = to_digits(pi, ndigits);

    std::cout << "pi = " << pi_str << std::endl;
    std::cout << std::endl;

    // Reference value for comparison
    std::string ref =
        "3."
        "14159265358979323846264338327950288419716939937510"
        "58209749445923078164062862089986280348253421170679"
        "82148086513282306647093844609550582231725359408128481117450284102701938521105559";

    std::cout << "ref= " << ref << std::endl;

    // Count matching digits
    int match = 0;
    for (size_t i = 0; i < pi_str.size() && i < ref.size(); i++) {
        if (pi_str[i] == ref[i]) match++;
        else break;
    }
    // Subtract 2 for "3."
    std::cout << std::endl << "Matching decimal digits: " << (match - 2) << std::endl;

    return 0;
}
