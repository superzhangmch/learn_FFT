#pragma once
#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>
#include <type_traits>

// ============================================================
// SmallUint<N>: sub-byte unsigned integer using N bits (1-7)
// Stored in a uint8_t, masked to N bits on every write.
// ============================================================
template<int N>
struct SmallUint {
    static_assert(N >= 1 && N <= 7, "SmallUint bit width must be 1-7");
    uint8_t val;
    static constexpr uint8_t MASK = (uint8_t)((1u << N) - 1);

    constexpr SmallUint() : val(0) {}
    constexpr SmallUint(int v)                : val((uint8_t)v & MASK) {}
    constexpr SmallUint(unsigned v)           : val((uint8_t)v & MASK) {}
    constexpr SmallUint(long v)               : val((uint8_t)v & MASK) {}
    constexpr SmallUint(unsigned long v)      : val((uint8_t)v & MASK) {}
    constexpr SmallUint(long long v)          : val((uint8_t)v & MASK) {}
    constexpr SmallUint(unsigned long long v) : val((uint8_t)v & MASK) {}
    constexpr SmallUint(double v)             : val((uint8_t)(int)v & MASK) {}
    template<int M>
    constexpr SmallUint(SmallUint<M> o)       : val(o.val & MASK) {}

    // Implicit widening to unsigned — lets all arithmetic "just work"
    constexpr operator unsigned() const { return val; }

    // Compound assignments (needed for in-place array element ops)
    SmallUint& operator|=(unsigned v)  { val = (uint8_t)(val | v) & MASK; return *this; }
    SmallUint& operator&=(unsigned v)  { val = (uint8_t)(val & v) & MASK; return *this; }
    SmallUint& operator<<=(int s)      { val = (uint8_t)(val << s) & MASK; return *this; }
    SmallUint& operator>>=(int s)      { val >>= s; return *this; }
};

using uint2_t = SmallUint<2>;
using uint4_t = SmallUint<4>;

// ============================================================
// LimbBits<T>: logical bit width of a limb type
//   Native types: sizeof(T)*8
//   SmallUint<N>: N
// ============================================================
template<typename T> struct LimbBits { static constexpr int value = sizeof(T) * 8; };
template<int N> struct LimbBits<SmallUint<N>> { static constexpr int value = N; };

// ============================================================
// DoubleWidth<T>: maps limb type to its double-width type
// ============================================================
template<typename T> struct DoubleWidth;
template<> struct DoubleWidth<SmallUint<2>> { using type = SmallUint<4>; };
template<> struct DoubleWidth<SmallUint<4>> { using type = uint8_t; };
template<> struct DoubleWidth<uint8_t>      { using type = uint16_t; };
template<> struct DoubleWidth<uint16_t>     { using type = uint32_t; };
template<> struct DoubleWidth<uint32_t>     { using type = uint64_t; };
template<> struct DoubleWidth<uint64_t>     { using type = __uint128_t; };

// ============================================================
// BigFloat<T, LIMBS>
//   T     = limb type (uint2_t .. uint64_t)
//   LIMBS = number of limbs
//   Mantissa = LIMBS * LimbBits<T> bits
// ============================================================
template<typename T, int LIMBS>
class BigFloat {
    using T2 = typename DoubleWidth<T>::type;

public:
    using value_type = T;
    static constexpr int TBITS   = LimbBits<T>::value;       // logical bits per limb
    static constexpr int BITS    = LIMBS * TBITS;             // total mantissa bits
    static constexpr int TOP_BIT = BITS - 1;
    static constexpr int BYTES   = LIMBS * (int)sizeof(T);    // physical storage

    enum Special { NORMAL, ZERO, INF, NAN_ };

    bool sign_;
    int32_t exp_;
    T mantissa_[LIMBS];
    Special special_;

    // --- Constructors ---
    BigFloat() : sign_(false), exp_(0), special_(ZERO) { memset(mantissa_, 0, BYTES); }

    BigFloat(double val) : exp_(0), special_(NORMAL) {
        memset(mantissa_, 0, BYTES);
        if (std::isnan(val)) { special_ = NAN_; sign_ = false; return; }
        if (std::isinf(val)) { special_ = INF; sign_ = val < 0; return; }
        if (val == 0.0) { special_ = ZERO; sign_ = std::signbit(val); return; }
        sign_ = val < 0;
        double abs_val = std::fabs(val);
        int exp;
        double frac = std::frexp(abs_val, &exp); // frac in [0.5, 1.0)
        // Extract 53 mantissa bits from double and place them individually
        uint64_t bits = (uint64_t)std::ldexp(frac, 53); // in [2^52, 2^53)
        for (int b = 52; b >= 0; b--) {
            if (bits & (1ULL << b)) {
                int mpos = TOP_BIT - (52 - b);
                if (mpos >= 0)
                    mantissa_[mpos / TBITS] |= T(T2(1) << (mpos % TBITS));
            }
        }
        exp_ = exp - BITS;
        normalize();
    }

    static BigFloat zero(bool neg = false) { BigFloat r; r.sign_ = neg; return r; }
    static BigFloat inf(bool neg = false)  { BigFloat r; r.sign_ = neg; r.special_ = INF; return r; }
    static BigFloat nan()                  { BigFloat r; r.special_ = NAN_; return r; }

    bool is_zero() const { return special_ == ZERO; }
    bool is_inf()  const { return special_ == INF; }
    bool is_nan()  const { return special_ == NAN_; }
    bool is_negative() const { return sign_; }

    BigFloat operator-() const { BigFloat r = *this; r.sign_ = !r.sign_; return r; }

    // --- Arithmetic ---
    BigFloat operator+(const BigFloat& rhs) const { return add_impl(rhs, false); }
    BigFloat operator-(const BigFloat& rhs) const { return add_impl(rhs, true); }

    BigFloat operator*(const BigFloat& rhs) const {
        bool rs = sign_ != rhs.sign_;
        if (is_nan() || rhs.is_nan()) return nan();
        if (is_inf() || rhs.is_inf()) {
            if (is_zero() || rhs.is_zero()) return nan();
            return inf(rs);
        }
        if (is_zero() || rhs.is_zero()) return zero(rs);

        T full[2 * LIMBS];
        m_mul_full(mantissa_, rhs.mantissa_, full);

        int hb = find_highest_bit(full, 2 * LIMBS);
        if (hb < 0) return zero(rs);

        BigFloat result;
        result.sign_ = rs;
        result.special_ = NORMAL;

        int sd = hb - TOP_BIT;
        if (sd >= 0) {
            extract_shifted(full, 2 * LIMBS, sd, result.mantissa_);
        } else {
            for (int i = 0; i < LIMBS; i++) result.mantissa_[i] = full[i];
            m_shift_left(result.mantissa_, -sd);
        }
        result.exp_ = exp_ + rhs.exp_ + (hb - TOP_BIT);
        result.normalize();
        return result;
    }

    BigFloat operator/(const BigFloat& rhs) const {
        bool rs = sign_ != rhs.sign_;
        if (is_nan() || rhs.is_nan()) return nan();
        if (is_inf() && rhs.is_inf()) return nan();
        if (is_inf()) return inf(rs);
        if (rhs.is_inf()) return zero(rs);
        if (rhs.is_zero()) { return is_zero() ? nan() : inf(rs); }
        if (is_zero()) return zero(rs);

        constexpr int QL = LIMBS + 1;
        T q[QL], rem[QL], div[QL];
        memset(q, 0, QL * sizeof(T));
        memset(rem, 0, QL * sizeof(T));
        memcpy(div, rhs.mantissa_, BYTES);
        div[LIMBS] = T(0);

        for (int i = 2 * BITS - 1; i >= 0; i--) {
            // Shift rem left by 1 (QL limbs)
            for (int k = QL - 1; k > 0; k--)
                rem[k] = T(((T2)rem[k] << 1) | ((T2)rem[k - 1] >> (TBITS - 1)));
            rem[0] <<= 1;

            if (i >= BITS) {
                int pos = i - BITS;
                rem[0] |= T(((T2)mantissa_[pos / TBITS] >> (pos % TBITS)) & 1);
            }

            // Compare rem >= div (QL limbs)
            bool ge = false;
            for (int k = QL - 1; k >= 0; k--) {
                if ((T2)rem[k] > (T2)div[k]) { ge = true; break; }
                if ((T2)rem[k] < (T2)div[k]) { ge = false; break; }
                if (k == 0) ge = true;
            }

            if (ge) {
                T borrow = T(0);
                for (int k = 0; k < QL; k++) {
                    T2 diff = (T2)rem[k] - (T2)div[k] - (T2)borrow;
                    rem[k] = (T)diff;
                    borrow = T((diff >> (2 * TBITS - 1)) ? 1 : 0);
                }
                q[i / TBITS] |= T(T2(1) << (i % TBITS));
            }
        }

        int q_hb = find_highest_bit(q, QL);
        if (q_hb < 0) return zero(rs);

        BigFloat result;
        result.sign_ = rs;
        result.special_ = NORMAL;

        int shift = q_hb - TOP_BIT;
        if (shift >= 0) {
            extract_shifted(q, QL, shift, result.mantissa_);
        } else {
            for (int i = 0; i < LIMBS; i++) result.mantissa_[i] = q[i];
            m_shift_left(result.mantissa_, -shift);
        }
        result.exp_ = exp_ - rhs.exp_ - BITS + shift;
        result.normalize();
        return result;
    }

    // --- Comparison ---
    bool operator==(const BigFloat& rhs) const {
        if (is_nan() || rhs.is_nan()) return false;
        if (is_zero() && rhs.is_zero()) return true;
        if (special_ != rhs.special_) return false;
        if (sign_ != rhs.sign_) return false;
        if (special_ == INF) return true;
        if (exp_ != rhs.exp_) return false;
        return m_compare(mantissa_, rhs.mantissa_) == 0;
    }
    bool operator!=(const BigFloat& rhs) const { if (is_nan() || rhs.is_nan()) return true; return !(*this == rhs); }
    bool operator<(const BigFloat& rhs) const {
        if (is_nan() || rhs.is_nan()) return false;
        if (is_zero() && rhs.is_zero()) return false;
        if (sign_ && !rhs.sign_) return true;
        if (!sign_ && rhs.sign_) return false;
        bool neg = sign_;
        if (is_inf() && rhs.is_inf()) return false;
        if (is_inf()) return neg;
        if (rhs.is_inf()) return !neg;
        if (is_zero()) return !neg;
        if (rhs.is_zero()) return neg;
        if (exp_ != rhs.exp_) return neg ? (exp_ > rhs.exp_) : (exp_ < rhs.exp_);
        int c = m_compare(mantissa_, rhs.mantissa_);
        return neg ? (c > 0) : (c < 0);
    }
    bool operator<=(const BigFloat& rhs) const { if (is_nan() || rhs.is_nan()) return false; return !(rhs < *this); }
    bool operator>(const BigFloat& rhs) const  { return rhs < *this; }
    bool operator>=(const BigFloat& rhs) const { if (is_nan() || rhs.is_nan()) return false; return !(*this < rhs); }

    // --- Display ---
    std::string to_string() const {
        if (is_nan()) return "NaN";
        if (is_inf()) return sign_ ? "-Inf" : "Inf";
        if (is_zero()) return sign_ ? "-0" : "0";
        std::string result;
        if (sign_) result += "-";
        // Read up to 53 bits from the top of the mantissa for display
        double sig = 0;
        int bits_read = 0;
        for (int i = LIMBS - 1; i >= 0 && bits_read < 53; i--) {
            sig += std::ldexp((double)(uint64_t)(T2)mantissa_[i], bits_read == 0 ? 0 : -bits_read);
            bits_read += TBITS;
        }
        sig /= std::ldexp(1.0, TBITS - 1);
        int total_exp2 = exp_ + TOP_BIT;
        char buf[64];
        if (total_exp2 > -1000 && total_exp2 < 1000) {
            snprintf(buf, sizeof(buf), "%.15g", std::ldexp(sig, total_exp2));
        } else {
            double l = std::log10(sig) + (double)total_exp2 * std::log10(2.0);
            int de = (int)std::floor(l);
            double s = std::pow(10.0, l - de);
            if (s >= 9.9999) { s = 1.0; de++; }
            char tmp[32]; snprintf(tmp, sizeof(tmp), "%.15g", s);
            snprintf(buf, sizeof(buf), "%se%+d", tmp, de);
        }
        result += buf;
        return result;
    }
    friend std::ostream& operator<<(std::ostream& os, const BigFloat& v) { return os << v.to_string(); }

private:
    void normalize() {
        if (special_ != NORMAL) return;
        if (m_is_zero(mantissa_)) { special_ = ZERO; return; }
        int hb = m_highest_bit(mantissa_);
        if (hb < TOP_BIT) { int s = TOP_BIT - hb; m_shift_left(mantissa_, s); exp_ -= s; }
        else if (hb > TOP_BIT) { int s = hb - TOP_BIT; m_shift_right(mantissa_, s); exp_ += s; }
    }

    BigFloat add_impl(const BigFloat& rhs, bool subtract) const {
        bool rhs_sign = subtract ? !rhs.sign_ : rhs.sign_;
        if (is_nan() || rhs.is_nan()) return nan();
        if (is_inf()) {
            if (rhs.is_inf() && sign_ != rhs_sign) return nan();
            return *this;
        }
        if (rhs.is_inf()) { BigFloat r = rhs; r.sign_ = rhs_sign; return r; }
        if (is_zero() && rhs.is_zero()) return zero(sign_ && rhs_sign);
        if (is_zero()) { BigFloat r = rhs; r.sign_ = rhs_sign; return r; }
        if (rhs.is_zero()) return *this;

        BigFloat a = *this, b = rhs;
        b.sign_ = rhs_sign;
        if (a.exp_ < b.exp_) std::swap(a, b);
        int ed = a.exp_ - b.exp_;
        if (ed > BITS) return a;

        T bm[LIMBS];
        memcpy(bm, b.mantissa_, BYTES);
        m_shift_right(bm, ed);

        BigFloat result;
        result.exp_ = a.exp_;
        if (a.sign_ == b.sign_) {
            int carry;
            m_add(a.mantissa_, bm, result.mantissa_, carry);
            result.sign_ = a.sign_;
            if (carry) {
                m_shift_right(result.mantissa_, 1);
                result.mantissa_[LIMBS - 1] |= T(T2(1) << (TBITS - 1));
                result.exp_++;
            }
        } else {
            int cmp = m_compare(a.mantissa_, bm);
            if (cmp == 0) return zero();
            if (cmp > 0) { m_sub(a.mantissa_, bm, result.mantissa_); result.sign_ = a.sign_; }
            else { m_sub(bm, a.mantissa_, result.mantissa_); result.sign_ = b.sign_; }
        }
        result.special_ = NORMAL;
        result.normalize();
        return result;
    }

    // --- Generic helpers ---
    static int find_highest_bit(const T arr[], int N) {
        for (int i = N - 1; i >= 0; i--)
            if ((T2)arr[i]) {
                int b = TBITS - 1;
                while (b >= 0 && !(((T2)arr[i] >> b) & 1)) b--;
                return i * TBITS + b;
            }
        return -1;
    }

    static void extract_shifted(const T src[], int src_n, int shift, T dst[]) {
        memset(dst, 0, BYTES);
        int lo = shift / TBITS, bo = shift % TBITS;
        for (int i = 0; i < LIMBS; i++) {
            int s = i + lo;
            if (s < src_n) {
                T2 v = (T2)src[s] >> bo;
                if (bo > 0 && s + 1 < src_n)
                    v |= (T2)src[s + 1] << (TBITS - bo);
                dst[i] = T(v);
            }
        }
    }

    // --- Mantissa helpers ---
    static bool m_is_zero(const T m[]) {
        for (int i = 0; i < LIMBS; i++) if ((T2)m[i]) return false;
        return true;
    }
    static int m_highest_bit(const T m[]) { return find_highest_bit(m, LIMBS); }

    static int m_compare(const T a[], const T b[]) {
        for (int i = LIMBS - 1; i >= 0; i--) {
            if ((T2)a[i] < (T2)b[i]) return -1;
            if ((T2)a[i] > (T2)b[i]) return 1;
        }
        return 0;
    }
    static void m_shift_right(T m[], int shift) {
        if (shift <= 0) return;
        if (shift >= BITS) { memset(m, 0, BYTES); return; }
        int ls = shift / TBITS, bs = shift % TBITS;
        for (int i = 0; i < LIMBS; i++) {
            int src = i + ls;
            if (src < LIMBS) {
                T2 v = (T2)m[src] >> bs;
                if (bs > 0 && src + 1 < LIMBS) v |= (T2)m[src + 1] << (TBITS - bs);
                m[i] = T(v);
            } else m[i] = T(0);
        }
    }
    static void m_shift_left(T m[], int shift) {
        if (shift <= 0) return;
        if (shift >= BITS) { memset(m, 0, BYTES); return; }
        int ls = shift / TBITS, bs = shift % TBITS;
        for (int i = LIMBS - 1; i >= 0; i--) {
            int src = i - ls;
            if (src >= 0) {
                T2 v = (T2)m[src] << bs;
                if (bs > 0 && src - 1 >= 0) v |= (T2)m[src - 1] >> (TBITS - bs);
                m[i] = T(v);
            } else m[i] = T(0);
        }
    }
    static void m_add(const T a[], const T b[], T r[], int& carry_out) {
        T carry = T(0);
        for (int i = 0; i < LIMBS; i++) {
            T2 s = (T2)a[i] + (T2)b[i] + (T2)carry;
            r[i] = (T)s;
            carry = T(s >> TBITS);
        }
        carry_out = (int)(T2)carry;
    }
    static void m_sub(const T a[], const T b[], T r[]) {
        T borrow = T(0);
        for (int i = 0; i < LIMBS; i++) {
            T2 d = (T2)a[i] - (T2)b[i] - (T2)borrow;
            r[i] = (T)d;
            borrow = T((d >> (2 * TBITS - 1)) ? 1 : 0);
        }
    }
    static void m_mul_full(const T a[], const T b[], T r[]) {
        memset(r, 0, 2 * BYTES);
        for (int i = 0; i < LIMBS; i++) {
            T carry = T(0);
            for (int j = 0; j < LIMBS; j++) {
                T2 p = (T2)a[i] * (T2)b[j] + (T2)r[i + j] + (T2)carry;
                r[i + j] = (T)p;
                carry = T(p >> TBITS);
            }
            r[i + LIMBS] = carry;
        }
    }
};

// === Convenience typedefs ===

// 64-bit limbs (fastest, needs __uint128_t)
using Float512  = BigFloat<uint64_t, 8>;
using Float1024 = BigFloat<uint64_t, 16>;
using Float2048 = BigFloat<uint64_t, 32>;

// 32-bit limbs (portable, no __uint128_t)
using Float512_32  = BigFloat<uint32_t, 16>;
using Float1024_32 = BigFloat<uint32_t, 32>;

// 8-bit limbs
using Float512_8  = BigFloat<uint8_t, 64>;
using Float1024_8 = BigFloat<uint8_t, 128>;

// 4-bit limbs
using Float512_4  = BigFloat<uint4_t, 128>;
using Float1024_4 = BigFloat<uint4_t, 256>;

// 2-bit limbs
using Float512_2  = BigFloat<uint2_t, 256>;
using Float1024_2 = BigFloat<uint2_t, 512>;
