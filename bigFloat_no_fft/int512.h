#pragma once
#include <cstdint>
#include <string>
#include <iostream>

class Int512 {
public:
    // 8 x 64-bit limbs, little-endian (data[0] = least significant)
    uint64_t data[8];

    // Constructors
    Int512();
    Int512(int val);
    Int512(int64_t val);
    explicit Int512(const char* decimal_str);
    explicit Int512(const std::string& decimal_str);

    // Check sign (two's complement: top bit of data[7])
    bool is_negative() const;
    bool is_zero() const;

    // Unary
    Int512 operator-() const;
    Int512 operator~() const;

    // Arithmetic
    Int512 operator+(const Int512& rhs) const;
    Int512 operator-(const Int512& rhs) const;
    Int512 operator*(const Int512& rhs) const;
    Int512 operator/(const Int512& rhs) const;
    Int512 operator%(const Int512& rhs) const;

    Int512& operator+=(const Int512& rhs);
    Int512& operator-=(const Int512& rhs);
    Int512& operator*=(const Int512& rhs);
    Int512& operator/=(const Int512& rhs);
    Int512& operator%=(const Int512& rhs);

    // Shift
    Int512 operator<<(int shift) const;
    Int512 operator>>(int shift) const;

    // Comparison
    bool operator==(const Int512& rhs) const;
    bool operator!=(const Int512& rhs) const;
    bool operator<(const Int512& rhs) const;
    bool operator<=(const Int512& rhs) const;
    bool operator>(const Int512& rhs) const;
    bool operator>=(const Int512& rhs) const;

    // Conversion
    std::string to_string() const;
    friend std::ostream& operator<<(std::ostream& os, const Int512& val);

    // Absolute value (returns unsigned interpretation helper)
    Int512 abs() const;

    // Get bit at position (0-511)
    int get_bit(int pos) const;
    void set_bit(int pos, int val);

    // Find highest set bit (-1 if zero)
    int highest_bit() const;

private:
    // Unsigned division: divides |dividend| by |divisor|, returns quotient and remainder
    static void unsigned_divmod(const Int512& dividend, const Int512& divisor,
                                Int512& quotient, Int512& remainder);
    // Unsigned comparison (treats both as positive)
    static int unsigned_compare(const Int512& a, const Int512& b);
};
