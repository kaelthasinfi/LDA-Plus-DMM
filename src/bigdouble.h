//#ifndef __BIGDOUBLE_H
//#define __BIGDOUBLE_H

#include <cassert>
#include <iostream>
#include <cmath>

class BigDouble {
public:
    BigDouble(): num(0), exp(0) {}
    BigDouble(double val): num(val), exp(0) { rescale(); }
    BigDouble(double num, int exp): num(num), exp(exp) { rescale(); }

    bool operator>(BigDouble rhs) const {
        if (num != 0 && rhs.num == 0) {
            return true;
        } else if (num == 0) {
            return false;
        } else if (exp > rhs.exp) {
            return true;
        } else if (exp < rhs.exp) {
            return false;
        } else {
            return num > rhs.num;
        }
    }

    BigDouble &operator+=(BigDouble rhs) {
        if (rhs.num == 0) return *this;
        if (num == 0) {
            num = rhs.num;
            exp = rhs.exp;
            return *this;
        }
        int expdiff = exp - rhs.exp;
        while (expdiff > 0) {
            rhs.num /= LARGE_DOUBLE;
            expdiff--;
        }
        while (expdiff < 0) {
            num /= LARGE_DOUBLE;
            expdiff++;
            exp++;
        }
        num += rhs.num;
        return *this;
    }

    BigDouble &operator*=(BigDouble rhs) {
        num *= rhs.num;
        exp += rhs.exp;
        if (num > LARGE_DOUBLE) {
            num /= LARGE_DOUBLE;
            exp++;
        }
        return *this;
    }

    BigDouble &operator/=(BigDouble rhs) {
        num /= rhs.num;
        exp -= rhs.exp;
        if (num < 1) {
            num *= LARGE_DOUBLE;
            exp--;
        }
        return *this;
    }

    BigDouble operator*(BigDouble rhs) const {
        BigDouble result = *this;
        result *= rhs;
        return result;
    }

    BigDouble operator/(BigDouble rhs) const {
        BigDouble result = *this;
        result /= rhs;
        return result;
    }

    explicit operator double() const {
        double val = num;
        int e = exp;
        while (e > 0) {
            val *= LARGE_DOUBLE;
            e--;
        }
        while (e < 0) {
            val /= LARGE_DOUBLE;
            e++;
        }
        return val;
    }

    void print() const {
        printf("(%f,%d)", num, exp);
    }

//private:
    constexpr static double LARGE_DOUBLE = 1e100;
    constexpr static double SMALL_DOUBLE = 1 / LARGE_DOUBLE;
    double num;
    int exp;

    void rescale() {
        assert(!std::isnan(num) && !std::isinf(num) && num >= 0);
        if (num == 0) exp = 0;
        while (num >= LARGE_DOUBLE) {
            num /= LARGE_DOUBLE;
            exp++;
        }
        while (num != 0 && num < 1) {
            num *= LARGE_DOUBLE;
            exp--;
        }
    }

    //friend inline std::ostream &operator<<(std::ostream &out, BigDouble d);
};
/*
inline std::ostream &operator<<(std::ostream &out, BigDouble d) {
    out << "(" << d.num << "," << d.exp << ")";
    return out;
}
*/
//#endif // __BIGDOUBLE_H
