#include "float512.h"
#include <cstring>
#include <cmath>
#include <algorithm>

// ============================================================
// Mantissa helpers (512-bit unsigned operations)
// ============================================================

bool Float512::m_is_zero(const uint64_t m[8]) {
    for (int i = 0; i < 8; i++)
        if (m[i] != 0) return false;
    return true;
}

int Float512::m_highest_bit(const uint64_t m[8]) {
    for (int i = 7; i >= 0; i--) {
        if (m[i] != 0) {
            int bit = 63;
            while (bit >= 0 && !((m[i] >> bit) & 1)) bit--;
            return i * 64 + bit;
        }
    }
    return -1;
}

int Float512::m_compare(const uint64_t a[8], const uint64_t b[8]) {
    for (int i = 7; i >= 0; i--) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

void Float512::m_shift_right(uint64_t m[8], int shift) {
    if (shift <= 0) return;
    if (shift >= 512) { memset(m, 0, 64); return; }

    int limb_shift = shift / 64;
    int bit_shift = shift % 64;

    for (int i = 0; i < 8; i++) {
        int src = i + limb_shift;
        if (src < 8) {
            m[i] = m[src] >> bit_shift;
            if (bit_shift > 0 && src + 1 < 8)
                m[i] |= m[src + 1] << (64 - bit_shift);
        } else {
            m[i] = 0;
        }
    }
}

void Float512::m_shift_left(uint64_t m[8], int shift) {
    if (shift <= 0) return;
    if (shift >= 512) { memset(m, 0, 64); return; }

    int limb_shift = shift / 64;
    int bit_shift = shift % 64;

    for (int i = 7; i >= 0; i--) {
        int src = i - limb_shift;
        if (src >= 0) {
            m[i] = m[src] << bit_shift;
            if (bit_shift > 0 && src - 1 >= 0)
                m[i] |= m[src - 1] >> (64 - bit_shift);
        } else {
            m[i] = 0;
        }
    }
}

void Float512::m_add(const uint64_t a[8], const uint64_t b[8], uint64_t result[8], int& carry_out) {
    uint64_t carry = 0;
    for (int i = 0; i < 8; i++) {
        __uint128_t sum = (__uint128_t)a[i] + b[i] + carry;
        result[i] = (uint64_t)sum;
        carry = (uint64_t)(sum >> 64);
    }
    carry_out = (int)carry;
}

void Float512::m_sub(const uint64_t a[8], const uint64_t b[8], uint64_t result[8]) {
    uint64_t borrow = 0;
    for (int i = 0; i < 8; i++) {
        __uint128_t diff = (__uint128_t)a[i] - b[i] - borrow;
        result[i] = (uint64_t)diff;
        borrow = (diff >> 127) ? 1 : 0; // check if underflowed
    }
}

void Float512::m_mul_full(const uint64_t a[8], const uint64_t b[8], uint64_t result[16]) {
    memset(result, 0, 128);
    for (int i = 0; i < 8; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 8; j++) {
            __uint128_t prod = (__uint128_t)a[i] * b[j] + result[i + j] + carry;
            result[i + j] = (uint64_t)prod;
            carry = (uint64_t)(prod >> 64);
        }
        result[i + 8] = carry;
    }
}

void Float512::m_div(const uint64_t dividend[8], const uint64_t divisor[8],
                     uint64_t quotient[8], uint64_t remainder[8]) {
    memset(quotient, 0, 64);
    memset(remainder, 0, 64);

    // Find highest bit of dividend
    int high = -1;
    for (int i = 7; i >= 0; i--) {
        if (dividend[i] != 0) {
            int bit = 63;
            while (bit >= 0 && !((dividend[i] >> bit) & 1)) bit--;
            high = i * 64 + bit;
            break;
        }
    }
    if (high < 0) return; // dividend is 0

    for (int i = high; i >= 0; i--) {
        // Left shift remainder by 1
        m_shift_left(remainder, 1);
        // Set bit 0 of remainder to bit i of dividend
        int limb = i / 64;
        int bit = i % 64;
        remainder[0] |= (dividend[limb] >> bit) & 1;

        // If remainder >= divisor, subtract
        if (m_compare(remainder, divisor) >= 0) {
            m_sub(remainder, divisor, remainder);
            quotient[limb] |= (uint64_t(1) << bit);
        }
    }
}

// ============================================================
// Constructors
// ============================================================

Float512::Float512() : sign_(false), exp_(0), special_(ZERO) {
    memset(mantissa_, 0, sizeof(mantissa_));
}

Float512 Float512::zero(bool negative) {
    Float512 r;
    r.sign_ = negative;
    r.special_ = ZERO;
    return r;
}

Float512 Float512::inf(bool negative) {
    Float512 r;
    r.sign_ = negative;
    r.special_ = INF;
    return r;
}

Float512 Float512::nan() {
    Float512 r;
    r.special_ = NAN_;
    return r;
}

Float512::Float512(double val) : exp_(0), special_(NORMAL) {
    memset(mantissa_, 0, sizeof(mantissa_));

    if (std::isnan(val)) { special_ = NAN_; sign_ = false; return; }
    if (std::isinf(val)) { special_ = INF; sign_ = val < 0; return; }
    if (val == 0.0) { special_ = ZERO; sign_ = std::signbit(val); return; }

    sign_ = val < 0;
    double abs_val = std::fabs(val);

    // Decompose: abs_val = frac * 2^exp, where 0.5 <= frac < 1.0
    int exp;
    double frac = std::frexp(abs_val, &exp);
    // frac is in [0.5, 1.0), so frac * 2^64 gives us the top 64 bits of mantissa
    // We want mantissa with bit 511 as the leading 1.
    // frac = mantissa * 2^(-512), so mantissa = frac * 2^512
    // and the actual value = sign * mantissa * 2^(exp - 512)

    // Place the double's 53 significant bits into the top of mantissa_
    // frac * 2^64 gives us the leading 64 bits
    uint64_t top_bits = (uint64_t)(frac * (double)(uint64_t(1) << 63)) << 1;
    // frac is in [0.5, 1.0), so frac * 2^63 is in [2^62, 2^63)
    // top_bits = that value << 1, putting the leading 1 at bit 63
    // Actually let's be more precise:
    // frac * 2^64 = value in [2^63, 2^64)
    top_bits = (uint64_t)(frac * 18446744073709551616.0); // frac * 2^64
    mantissa_[7] = top_bits;
    // remaining bits are below double precision, leave as 0

    // exp_ = actual exponent such that value = mantissa * 2^exp_
    // where mantissa has its leading 1 at bit 511
    // frac = mantissa * 2^(-512), so val = frac * 2^exp = mantissa * 2^(exp - 512)
    exp_ = exp - 512;

    // Verify leading bit is at position 511
    normalize();
}

Float512::Float512(const Int512& val) : exp_(0), special_(NORMAL) {
    memset(mantissa_, 0, sizeof(mantissa_));

    if (val.is_zero()) { special_ = ZERO; sign_ = false; return; }

    sign_ = val.is_negative();
    Int512 abs_val = val.abs();

    // Copy absolute value into mantissa
    for (int i = 0; i < 8; i++)
        mantissa_[i] = abs_val.data[i];

    // Find highest bit
    int hb = m_highest_bit(mantissa_);
    // Shift so that the leading 1 is at bit 511
    exp_ = hb - 511;
    if (hb < 511) {
        m_shift_left(mantissa_, 511 - hb);
    } else if (hb > 511) {
        m_shift_right(mantissa_, hb - 511);
    }
}

// ============================================================
// Normalize
// ============================================================

void Float512::normalize() {
    if (special_ != NORMAL) return;

    if (m_is_zero(mantissa_)) {
        special_ = ZERO;
        return;
    }

    int hb = m_highest_bit(mantissa_);
    if (hb < 511) {
        int shift = 511 - hb;
        m_shift_left(mantissa_, shift);
        exp_ -= shift;
    } else if (hb > 511) {
        int shift = hb - 511;
        m_shift_right(mantissa_, shift);
        exp_ += shift;
    }
}

// ============================================================
// Unary negation
// ============================================================

Float512 Float512::operator-() const {
    Float512 r = *this;
    r.sign_ = !r.sign_;
    return r;
}

// ============================================================
// Addition / Subtraction
// ============================================================

Float512 Float512::add_impl(const Float512& rhs, bool subtract) const {
    bool rhs_sign = subtract ? !rhs.sign_ : rhs.sign_;

    // Handle special values
    if (is_nan() || rhs.is_nan()) return nan();
    if (is_inf()) {
        if (rhs.is_inf() && sign_ != rhs_sign) return nan(); // inf - inf
        return *this;
    }
    if (rhs.is_inf()) {
        Float512 r = rhs;
        r.sign_ = rhs_sign;
        return r;
    }
    if (is_zero() && rhs.is_zero()) {
        return zero(sign_ && rhs_sign); // -0 + -0 = -0
    }
    if (is_zero()) {
        Float512 r = rhs;
        r.sign_ = rhs_sign;
        return r;
    }
    if (rhs.is_zero()) return *this;

    // Align exponents
    Float512 a = *this;
    Float512 b = rhs;
    b.sign_ = rhs_sign;

    // Make 'a' the one with the larger exponent
    if (a.exp_ < b.exp_) std::swap(a, b);

    int exp_diff = a.exp_ - b.exp_;
    if (exp_diff > 512) {
        // b is negligible
        return a;
    }

    // Shift b's mantissa right to align
    uint64_t b_mantissa[8];
    memcpy(b_mantissa, b.mantissa_, sizeof(b_mantissa));
    m_shift_right(b_mantissa, exp_diff);

    Float512 result;
    result.exp_ = a.exp_;

    if (a.sign_ == b.sign_) {
        // Same sign: add mantissas
        int carry;
        m_add(a.mantissa_, b_mantissa, result.mantissa_, carry);
        result.sign_ = a.sign_;

        if (carry) {
            // Overflow: shift right by 1, set top bit
            m_shift_right(result.mantissa_, 1);
            result.mantissa_[7] |= (uint64_t(1) << 63);
            result.exp_++;
        }
    } else {
        // Different signs: subtract smaller from larger
        int cmp = m_compare(a.mantissa_, b_mantissa);
        if (cmp == 0) {
            return zero();
        } else if (cmp > 0) {
            m_sub(a.mantissa_, b_mantissa, result.mantissa_);
            result.sign_ = a.sign_;
        } else {
            m_sub(b_mantissa, a.mantissa_, result.mantissa_);
            result.sign_ = b.sign_;
        }
    }

    result.special_ = NORMAL;
    result.normalize();
    return result;
}

Float512 Float512::operator+(const Float512& rhs) const {
    return add_impl(rhs, false);
}

Float512 Float512::operator-(const Float512& rhs) const {
    return add_impl(rhs, true);
}

// ============================================================
// Multiplication
// ============================================================

Float512 Float512::operator*(const Float512& rhs) const {
    bool result_sign = sign_ != rhs.sign_;

    // Special values
    if (is_nan() || rhs.is_nan()) return nan();
    if (is_inf() || rhs.is_inf()) {
        if (is_zero() || rhs.is_zero()) return nan(); // 0 * inf
        return inf(result_sign);
    }
    if (is_zero() || rhs.is_zero()) return zero(result_sign);

    // Multiply mantissas (512x512 -> 1024 bits)
    uint64_t full[16];
    m_mul_full(mantissa_, rhs.mantissa_, full);

    // The product has its leading 1 at bit 1022 or 1023 (since each mantissa has leading 1 at bit 511)
    // 511 + 511 = 1022, or 1023 if there's no cancellation
    // We want the top 512 bits

    Float512 result;
    result.sign_ = result_sign;
    result.special_ = NORMAL;

    // Find highest bit in full[16]
    int hb = -1;
    for (int i = 15; i >= 0; i--) {
        if (full[i] != 0) {
            int bit = 63;
            while (bit >= 0 && !((full[i] >> bit) & 1)) bit--;
            hb = i * 64 + bit;
            break;
        }
    }

    if (hb < 0) return zero(result_sign);

    // Extract top 512 bits starting from hb
    // We want bits [hb..hb-511] placed at [511..0]
    int shift_down = hb - 511;
    if (shift_down >= 0) {
        // Shift the 1024-bit result right by shift_down
        // and take the lower 512 bits
        int limb_shift = shift_down / 64;
        int bit_shift = shift_down % 64;
        for (int i = 0; i < 8; i++) {
            int src = i + limb_shift;
            if (src < 16) {
                result.mantissa_[i] = full[src] >> bit_shift;
                if (bit_shift > 0 && src + 1 < 16)
                    result.mantissa_[i] |= full[src + 1] << (64 - bit_shift);
            }
        }
    } else {
        // hb < 511, need to shift left
        int shift_up = -shift_down;
        for (int i = 0; i < 8; i++)
            result.mantissa_[i] = full[i];
        m_shift_left(result.mantissa_, shift_up);
    }

    // value = (m_a * 2^exp_a) * (m_b * 2^exp_b) = m_a*m_b * 2^(exp_a+exp_b)
    // m_a*m_b has leading 1 at bit hb. We shifted right by (hb-511) to normalize.
    // So result_mantissa = m_a*m_b / 2^(hb-511), compensate in exponent:
    result.exp_ = exp_ + rhs.exp_ + (hb - 511);

    result.normalize();
    return result;
}

// ============================================================
// Division
// ============================================================

Float512 Float512::operator/(const Float512& rhs) const {
    bool result_sign = sign_ != rhs.sign_;

    // Special values
    if (is_nan() || rhs.is_nan()) return nan();
    if (is_inf() && rhs.is_inf()) return nan();
    if (is_inf()) return inf(result_sign);
    if (rhs.is_inf()) return zero(result_sign);
    if (rhs.is_zero()) {
        if (is_zero()) return nan(); // 0/0
        return inf(result_sign);     // x/0
    }
    if (is_zero()) return zero(result_sign);

    // Compute (mantissa_a << 512) / mantissa_b via binary long division.
    // Both mantissas have leading 1 at bit 511, so:
    //   mantissa_a / mantissa_b ∈ [0.5, 2.0)
    //   quotient = (mantissa_a << 512) / mantissa_b ∈ [2^511, 2^513)
    // The quotient needs at most 513 bits. We store it in 9 limbs.

    uint64_t q[9];   // 576-bit quotient (only ~513 bits used)
    uint64_t rem[9];  // 576-bit remainder (needs 513 bits: after shift, before subtract)
    memset(q, 0, sizeof(q));
    memset(rem, 0, sizeof(rem));

    // Extend divisor to 9 limbs for comparison
    uint64_t div9[9];
    memcpy(div9, rhs.mantissa_, 64);
    div9[8] = 0;

    // 1024-bit dividend = mantissa_ << 512
    // Bits 1023..512 are mantissa_[0..7], bits 511..0 are zero.
    for (int i = 1023; i >= 0; i--) {
        // Shift remainder left by 1 (9 limbs)
        for (int k = 8; k > 0; k--)
            rem[k] = (rem[k] << 1) | (rem[k-1] >> 63);
        rem[0] <<= 1;

        // Bring in bit i of the 1024-bit dividend
        if (i >= 512) {
            int pos = i - 512;
            rem[0] |= (mantissa_[pos / 64] >> (pos % 64)) & 1;
        }

        // Trial subtraction (compare 9-limb rem with 9-limb div9)
        bool ge = false;
        for (int k = 8; k >= 0; k--) {
            if (rem[k] > div9[k]) { ge = true; break; }
            if (rem[k] < div9[k]) { ge = false; break; }
            if (k == 0) ge = true; // equal
        }

        if (ge) {
            // Subtract div9 from rem (9 limbs)
            uint64_t borrow = 0;
            for (int k = 0; k < 9; k++) {
                __uint128_t diff = (__uint128_t)rem[k] - div9[k] - borrow;
                rem[k] = (uint64_t)diff;
                borrow = (diff >> 127) ? 1 : 0;
            }
            q[i / 64] |= (uint64_t(1) << (i % 64));
        }
    }

    // Find the leading 1 in the quotient (should be at bit 511 or 512)
    int q_hb = -1;
    for (int i = 8; i >= 0; i--) {
        if (q[i] != 0) {
            int bit = 63;
            while (bit >= 0 && !((q[i] >> bit) & 1)) bit--;
            q_hb = i * 64 + bit;
            break;
        }
    }

    if (q_hb < 0) return zero(result_sign);

    // Extract 512 bits with leading 1 at bit 511
    Float512 result;
    result.sign_ = result_sign;
    result.special_ = NORMAL;
    memset(result.mantissa_, 0, sizeof(result.mantissa_));

    int shift = q_hb - 511; // how many bits to shift right (0 or 1 typically)
    if (shift >= 0) {
        int limb_off = shift / 64;
        int bit_off = shift % 64;
        for (int i = 0; i < 8; i++) {
            int src = i + limb_off;
            if (src <= 8) {
                result.mantissa_[i] = q[src] >> bit_off;
                if (bit_off > 0 && src + 1 <= 8)
                    result.mantissa_[i] |= q[src + 1] << (64 - bit_off);
            }
        }
    } else {
        // q_hb < 511: shift left
        int up = -shift;
        for (int i = 0; i < 8; i++)
            result.mantissa_[i] = q[i];
        m_shift_left(result.mantissa_, up);
    }

    // value = mantissa_a * 2^exp_a / (mantissa_b * 2^exp_b)
    //       = (mantissa_a / mantissa_b) * 2^(exp_a - exp_b)
    //       = (quotient / 2^512) * 2^(exp_a - exp_b)
    //       = quotient * 2^(exp_a - exp_b - 512)
    // We normalized quotient to have leading 1 at bit 511 by shifting right by 'shift',
    // so result_mantissa = quotient >> shift, and:
    // value = result_mantissa * 2^(exp_a - exp_b - 512 + shift)
    result.exp_ = exp_ - rhs.exp_ - 512 + shift;

    result.normalize();
    return result;
}

// ============================================================
// Comparison
// ============================================================

bool Float512::operator==(const Float512& rhs) const {
    if (is_nan() || rhs.is_nan()) return false;
    if (is_zero() && rhs.is_zero()) return true;
    if (special_ != rhs.special_) return false;
    if (sign_ != rhs.sign_) return false;
    if (special_ == INF) return true;
    if (exp_ != rhs.exp_) return false;
    return m_compare(mantissa_, rhs.mantissa_) == 0;
}

bool Float512::operator!=(const Float512& rhs) const {
    if (is_nan() || rhs.is_nan()) return true;
    return !(*this == rhs);
}

bool Float512::operator<(const Float512& rhs) const {
    if (is_nan() || rhs.is_nan()) return false;
    if (is_zero() && rhs.is_zero()) return false;

    // -x < +x
    if (sign_ && !rhs.sign_) {
        if (is_zero() && rhs.is_zero()) return false;
        return true;
    }
    if (!sign_ && rhs.sign_) return false;

    // Both positive or both negative
    bool both_negative = sign_;

    // Handle infinities
    if (is_inf() && rhs.is_inf()) return false;
    if (is_inf()) return both_negative;  // +inf not less, -inf is less
    if (rhs.is_inf()) return !both_negative;

    if (is_zero()) return !both_negative; // 0 < positive, not less than negative
    if (rhs.is_zero()) return both_negative;

    // Compare by exponent first, then mantissa
    if (exp_ != rhs.exp_) {
        return both_negative ? (exp_ > rhs.exp_) : (exp_ < rhs.exp_);
    }

    int cmp = m_compare(mantissa_, rhs.mantissa_);
    return both_negative ? (cmp > 0) : (cmp < 0);
}

bool Float512::operator<=(const Float512& rhs) const {
    if (is_nan() || rhs.is_nan()) return false;
    return !(rhs < *this);
}

bool Float512::operator>(const Float512& rhs) const {
    return rhs < *this;
}

bool Float512::operator>=(const Float512& rhs) const {
    if (is_nan() || rhs.is_nan()) return false;
    return !(*this < rhs);
}

// ============================================================
// String output
// ============================================================

std::string Float512::to_string() const {
    if (is_nan()) return "NaN";
    if (is_inf()) return sign_ ? "-Inf" : "Inf";
    if (is_zero()) return sign_ ? "-0" : "0";

    // Convert to decimal: value = mantissa * 2^exp_
    // mantissa is a 512-bit integer with leading 1 at bit 511
    // We'll compute the integer and fractional parts separately

    std::string result;
    if (sign_) result += "-";

    // Effective exponent: exp_ means mantissa * 2^exp_
    // mantissa has leading 1 at bit 511, so value = mantissa * 2^exp_
    // = (mantissa / 2^511) * 2^(exp_ + 511)

    // Approximate as double for display.
    // mantissa has leading 1 at bit 511. Normalize to [1.0, 2.0):
    //   significand_dbl = mantissa / 2^511 (in [1.0, 2.0))
    //   value = significand_dbl * 2^(exp_ + 511)
    double sig_dbl = (double)mantissa_[7] / (double)(uint64_t(1) << 63);
    // sig_dbl ≈ mantissa_[7] / 2^63 ≈ mantissa / 2^(63 + 7*64) = mantissa / 2^511
    // This captures ~53 bits of precision.

    int total_exp2 = exp_ + 511; // value = sig_dbl * 2^total_exp2

    // Try direct double conversion if exponent is in range
    char buf[64];
    if (total_exp2 > -1000 && total_exp2 < 1000) {
        double val = std::ldexp(sig_dbl, total_exp2);
        snprintf(buf, sizeof(buf), "%.15g", val);
    } else {
        // For extreme exponents, use log10 approach
        double log10_val = std::log10(sig_dbl) + (double)total_exp2 * std::log10(2.0);
        int decimal_exp = (int)std::floor(log10_val);
        double significand = std::pow(10.0, log10_val - decimal_exp);
        if (significand >= 9.9999) { significand = 1.0; decimal_exp++; }
        if (significand < 1.0) { significand *= 10.0; decimal_exp--; }
        char tmp[32];
        snprintf(tmp, sizeof(tmp), "%.15g", significand);
        snprintf(buf, sizeof(buf), "%se%+d", tmp, decimal_exp);
    }
    result += buf;

    return result;
}

std::ostream& operator<<(std::ostream& os, const Float512& val) {
    return os << val.to_string();
}
