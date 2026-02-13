#pragma once
#include <cstdint>
#include <string>
#include <iostream>
#include "int512.h"

class Float512 {
public:
    enum Special { NORMAL, ZERO, INF, NAN_ };

    bool sign_;             // true = negative
    int32_t exp_;           // unbiased exponent
    uint64_t mantissa_[8];  // 512-bit mantissa, bit 511 is the implicit leading 1 position
    Special special_;

    // Constructors
    Float512();
    Float512(double val);
    Float512(const Int512& val);

    // Named constructors for special values
    static Float512 zero(bool negative = false);
    static Float512 inf(bool negative = false);
    static Float512 nan();

    // Queries
    bool is_zero() const { return special_ == ZERO; }
    bool is_inf() const { return special_ == INF; }
    bool is_nan() const { return special_ == NAN_; }
    bool is_negative() const { return sign_; }

    // Arithmetic
    Float512 operator+(const Float512& rhs) const;
    Float512 operator-(const Float512& rhs) const;
    Float512 operator*(const Float512& rhs) const;
    Float512 operator/(const Float512& rhs) const;

    Float512 operator-() const;

    // Comparison
    bool operator==(const Float512& rhs) const;
    bool operator!=(const Float512& rhs) const;
    bool operator<(const Float512& rhs) const;
    bool operator<=(const Float512& rhs) const;
    bool operator>(const Float512& rhs) const;
    bool operator>=(const Float512& rhs) const;

    // Output
    std::string to_string() const;
    friend std::ostream& operator<<(std::ostream& os, const Float512& val);

private:
    // Normalize: shift mantissa so bit 511 is 1, adjust exponent
    void normalize();

    // Mantissa helpers (512-bit unsigned operations on mantissa_[8])
    static void m_shift_right(uint64_t m[8], int shift);
    static void m_shift_left(uint64_t m[8], int shift);
    static int  m_compare(const uint64_t a[8], const uint64_t b[8]);
    static void m_add(const uint64_t a[8], const uint64_t b[8], uint64_t result[8], int& carry_out);
    static void m_sub(const uint64_t a[8], const uint64_t b[8], uint64_t result[8]); // a >= b assumed
    static int  m_highest_bit(const uint64_t m[8]);
    static bool m_is_zero(const uint64_t m[8]);

    // Full 1024-bit multiply for mantissa
    static void m_mul_full(const uint64_t a[8], const uint64_t b[8], uint64_t result[16]);

    // 512-bit / 512-bit division giving 512-bit quotient
    static void m_div(const uint64_t dividend[8], const uint64_t divisor[8],
                      uint64_t quotient[8], uint64_t remainder[8]);

    // Addition helper (handles sign and exponent alignment)
    Float512 add_impl(const Float512& rhs, bool subtract) const;
};
