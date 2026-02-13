#pragma once
#include "bigfloat.h"
#include <string>

// ============================================================
// Utility helpers (from pi1024.cpp, templatized)
// ============================================================

template<typename BF>
int bf_floor_int(const BF& val) {
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
std::string bf_to_digits(BF val, int ndigits) {
    std::string result;
    if (val.is_negative()) { result += '-'; val = -val; }
    int ipart = bf_floor_int(val);
    result += std::to_string(ipart);
    result += '.';
    BF frac = val - BF((double)ipart);
    BF ten(10.0);
    for (int i = 0; i < ndigits; i++) {
        frac = frac * ten;
        int d = bf_floor_int(frac);
        if (d < 0 || d > 9) d = 0;
        result += ('0' + d);
        frac = frac - BF((double)d);
    }
    return result;
}

// ============================================================
// bf_abs
// ============================================================

template<typename BF>
BF bf_abs(const BF& x) {
    if (x.is_negative()) return -x;
    return x;
}

// ============================================================
// bf_sqrt — Newton iteration  (from pi1024.cpp)
// ============================================================

template<typename BF>
BF bf_sqrt(const BF& S) {
    if (S.is_zero()) return S;
    if (S.is_negative()) return BF::nan();
    constexpr int TL = BF::BITS / BF::TBITS - 1;
    using T2 = typename DoubleWidth<typename BF::value_type>::type;
    double top = (double)(uint64_t)(T2)S.mantissa_[TL];
    double s_approx = std::ldexp(top / std::ldexp(1.0, BF::TBITS - 1), S.exp_ + BF::TOP_BIT);
    BF x(std::sqrt(s_approx));
    BF half(0.5);
    for (int i = 0; i < 12; i++)
        x = (x + S / x) * half;
    return x;
}

// ============================================================
// bf_pi — Chudnovsky  (from pi1024.cpp)
// ============================================================

template<typename BF>
BF bf_pi() {
    constexpr int nterms = BF::BITS / 47 + 2; // ~14 digits per term
    BF c640320(640320.0);
    BF C = c640320 * c640320 * c640320;
    BF sum(13591409.0);
    BF term(1.0);

    for (int k = 1; k <= nterms; k++) {
        long long k6 = 6LL * k, k3 = 3LL * k;
        BF num(1.0);
        for (int j = 0; j < 6; j++)
            num = num * BF((double)(k6 - j));
        BF den = BF((double)k3) * BF((double)(k3 - 1)) * BF((double)(k3 - 2));
        den = den * BF((double)k) * BF((double)k) * BF((double)k) * C;
        term = term * num / den;
        term = -term;
        sum = sum + term * BF(13591409.0 + 545140134.0 * k);
    }

    BF sqrt10005 = bf_sqrt(BF(10005.0));
    return BF(426880.0) * sqrt10005 / sum;
}

// ============================================================
// bf_exp — argument reduction x/2^n, Taylor, square back
// ============================================================

template<typename BF>
BF bf_exp(const BF& x) {
    if (x.is_zero()) return BF(1.0);
    if (x.is_nan()) return BF::nan();
    if (x.is_inf()) return x.is_negative() ? BF(0.0) : BF::inf();

    // Handle negative: exp(-x) = 1/exp(x)
    if (x.is_negative()) return BF(1.0) / bf_exp(-x);

    // Argument reduction: divide by 2^n until |r| < 1
    BF r = x;
    int n = 0;
    BF one(1.0);
    BF half(0.5);
    while (r > one) {
        r = r * half;
        n++;
    }

    // Taylor series: exp(r) = sum r^k / k!
    constexpr int max_iter = BF::BITS / 2 + 50;
    BF sum_val(1.0);
    BF term(1.0);
    for (int k = 1; k <= max_iter; k++) {
        term = term * r / BF((double)k);
        BF prev = sum_val;
        sum_val = sum_val + term;
        if (sum_val == prev) break;
    }

    // Square back n times
    for (int i = 0; i < n; i++)
        sum_val = sum_val * sum_val;

    return sum_val;
}

// ============================================================
// bf_ln — Newton: y_{n+1} = y_n + x*exp(-y_n) - 1
// ============================================================

template<typename BF>
BF bf_ln(const BF& x) {
    if (x.is_zero()) return BF::inf(true); // -inf
    if (x.is_negative()) return BF::nan();
    if (x.is_nan()) return BF::nan();
    if (x.is_inf()) return BF::inf();

    BF one(1.0);
    // Initial guess from double
    using T2 = typename DoubleWidth<typename BF::value_type>::type;
    constexpr int TL = BF::BITS / BF::TBITS - 1;
    double top = (double)(uint64_t)(T2)x.mantissa_[TL];
    double x_approx = std::ldexp(top / std::ldexp(1.0, BF::TBITS - 1), x.exp_ + BF::TOP_BIT);
    BF y(std::log(x_approx));

    // Newton iterations: y_{n+1} = y_n + x * exp(-y_n) - 1
    for (int i = 0; i < 14; i++) {
        BF ey = bf_exp(-y);
        y = y + x * ey - one;
    }
    return y;
}

// ============================================================
// bf_sin / bf_cos — argument reduction mod 2pi, Taylor
// ============================================================

template<typename BF>
BF bf_sin(const BF& x) {
    if (x.is_zero()) return x;
    if (x.is_nan() || x.is_inf()) return BF::nan();

    // sin(-x) = -sin(x)
    if (x.is_negative()) return -bf_sin(-x);

    // Reduce x mod 2*pi
    BF pi = bf_pi<BF>();
    BF two_pi = pi * BF(2.0);
    BF half_pi = pi * BF(0.5);

    BF r = x;
    // Subtract multiples of 2*pi
    if (r > two_pi) {
        // Compute floor(r / two_pi) via double approximation
        using T2 = typename DoubleWidth<typename BF::value_type>::type;
        constexpr int TL = BF::BITS / BF::TBITS - 1;
        double rtop = (double)(uint64_t)(T2)r.mantissa_[TL];
        double r_approx = std::ldexp(rtop / std::ldexp(1.0, BF::TBITS - 1), r.exp_ + BF::TOP_BIT);
        double tp_approx = 2.0 * 3.14159265358979323846;
        int q = (int)(r_approx / tp_approx);
        r = r - two_pi * BF((double)q);
        while (r > two_pi) r = r - two_pi;
        while (r.is_negative()) r = r + two_pi;
    }

    // Now r in [0, 2*pi). Fold into [0, pi/2]
    bool negate = false;
    if (r > pi) {          // sin(pi+x) = -sin(x)
        r = r - pi;
        negate = true;
    }
    if (r > half_pi) {     // sin(pi-x) = sin(x)
        r = pi - r;
    }

    // Taylor series: sin(r) = r - r^3/3! + r^5/5! - ...
    // term *= -r^2 / ((2k)(2k+1))
    BF r2 = r * r;
    BF sum_val = r;
    BF term = r;
    constexpr int max_iter = BF::BITS / 2 + 50;
    for (int k = 1; k <= max_iter; k++) {
        term = term * r2 / BF((double)(2 * k) * (2 * k + 1));
        term = -term;
        BF prev = sum_val;
        sum_val = sum_val + term;
        if (sum_val == prev) break;
    }

    return negate ? -sum_val : sum_val;
}

template<typename BF>
BF bf_cos(const BF& x) {
    if (x.is_zero()) return BF(1.0);
    if (x.is_nan() || x.is_inf()) return BF::nan();

    // cos(-x) = cos(x)
    BF ax = bf_abs(x);

    // Reduce mod 2*pi
    BF pi = bf_pi<BF>();
    BF two_pi = pi * BF(2.0);
    BF half_pi = pi * BF(0.5);

    BF r = ax;
    if (r > two_pi) {
        using T2 = typename DoubleWidth<typename BF::value_type>::type;
        constexpr int TL = BF::BITS / BF::TBITS - 1;
        double rtop = (double)(uint64_t)(T2)r.mantissa_[TL];
        double r_approx = std::ldexp(rtop / std::ldexp(1.0, BF::TBITS - 1), r.exp_ + BF::TOP_BIT);
        double tp_approx = 2.0 * 3.14159265358979323846;
        int q = (int)(r_approx / tp_approx);
        r = r - two_pi * BF((double)q);
        while (r > two_pi) r = r - two_pi;
        while (r.is_negative()) r = r + two_pi;
    }

    // Fold into [0, pi/2]
    bool negate = false;
    if (r > pi) {          // cos(2pi-x) = cos(x)
        r = two_pi - r;
    }
    if (r > half_pi) {     // cos(pi-x) = -cos(x)
        r = pi - r;
        negate = true;
    }

    // Taylor series: cos(r) = 1 - r^2/2! + r^4/4! - ...
    // term *= -r^2 / ((2k-1)(2k))
    BF r2 = r * r;
    BF sum_val(1.0);
    BF term(1.0);
    constexpr int max_iter = BF::BITS / 2 + 50;
    for (int k = 1; k <= max_iter; k++) {
        term = term * r2 / BF((double)(2 * k - 1) * (2 * k));
        term = -term;
        BF prev = sum_val;
        sum_val = sum_val + term;
        if (sum_val == prev) break;
    }

    return negate ? -sum_val : sum_val;
}

// ============================================================
// bf_tan = sin / cos
// ============================================================

template<typename BF>
BF bf_tan(const BF& x) {
    return bf_sin(x) / bf_cos(x);
}

// ============================================================
// bf_atan — argument reduction + Taylor
// ============================================================

template<typename BF>
BF bf_atan(const BF& x) {
    if (x.is_zero()) return x;
    if (x.is_nan()) return BF::nan();
    if (x.is_inf()) {
        BF half_pi = bf_pi<BF>() * BF(0.5);
        return x.is_negative() ? -half_pi : half_pi;
    }

    // atan(-x) = -atan(x)
    if (x.is_negative()) return -bf_atan(-x);

    BF one(1.0);
    BF half(0.5);

    // If |x| > 1: atan(x) = pi/2 - atan(1/x)
    if (x > one) {
        BF half_pi = bf_pi<BF>() * half;
        return half_pi - bf_atan(one / x);
    }

    // If |x| > 0.5: half-angle formula
    // atan(x) = 2*atan(x / (1 + sqrt(1+x^2)))
    if (x > half) {
        BF x2 = x * x;
        BF denom = one + bf_sqrt(one + x2);
        return BF(2.0) * bf_atan(x / denom);
    }

    // Taylor series: atan(x) = x - x^3/3 + x^5/5 - ...
    // term *= -x^2 * (2k-1) / (2k+1)
    BF x2 = x * x;
    BF sum_val = x;
    BF term = x;
    constexpr int max_iter = BF::BITS / 2 + 50;
    for (int k = 1; k <= max_iter; k++) {
        term = term * x2 * BF((double)(2 * k - 1)) / BF((double)(2 * k + 1));
        term = -term;
        BF prev = sum_val;
        sum_val = sum_val + term;
        if (sum_val == prev) break;
    }

    return sum_val;
}

// ============================================================
// bf_asin — asin(x) = atan(x / sqrt(1 - x^2))
// ============================================================

template<typename BF>
BF bf_asin(const BF& x) {
    if (x.is_zero()) return x;
    if (x.is_nan()) return BF::nan();
    if (x.is_inf()) return BF::nan();

    // asin(-x) = -asin(x)
    if (x.is_negative()) return -bf_asin(-x);

    BF one(1.0);

    // |x| > 1 is out of domain
    if (x > one) return BF::nan();

    // asin(1) = pi/2
    if (x == one) return bf_pi<BF>() * BF(0.5);

    // asin(x) = atan(x / sqrt(1 - x^2))
    BF x2 = x * x;
    return bf_atan(x / bf_sqrt(one - x2));
}

// ============================================================
// bf_acos — acos(x) = pi/2 - asin(x)
// ============================================================

template<typename BF>
BF bf_acos(const BF& x) {
    if (x.is_nan()) return BF::nan();
    if (x.is_inf()) return BF::nan();

    BF one(1.0);
    BF ax = bf_abs(x);
    if (ax > one) return BF::nan();

    return bf_pi<BF>() * BF(0.5) - bf_asin(x);
}
