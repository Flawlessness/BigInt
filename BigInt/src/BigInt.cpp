#include "BigInt.h"
#include <string_view>
#include <charconv>

namespace
{
    const int32_t SIZE_NUM_DIGIT = 9;
    const int32_t MAX_NUM = 1000000000;
}

BigInt BigInt::convert_dec_to_bin_mod(BigInt decimal_number)
{
    BigInt binary_number(0);
    BigInt factor(1);
    decimal_number.sign = 1;
    const BigInt number_zero = BigInt(0);
    while (decimal_number > number_zero)
    {
        binary_number += (decimal_number % BigInt(2)) * factor;
        factor *= 10;
        decimal_number /= 2;
    }
    return binary_number;
}

BigInt BigInt::convert_bin_to_dec_mod(BigInt binary_number)
{
    BigInt decimal_number(0);
    BigInt factor(1);
    binary_number.sign = 1;
    const BigInt number_zero = BigInt(0);
    while (binary_number > number_zero)
    {
        decimal_number += (binary_number % BigInt(2)) * factor;
        factor *= 2;
        binary_number /= 10;
    }
    return decimal_number;
}

BigInt BigInt::plus_one(BigInt binary_number)
{
    binary_number = convert_bin_to_dec_mod(binary_number);
    binary_number++;
    binary_number = convert_dec_to_bin_mod(binary_number);
    return binary_number;
}

BigInt BigInt::Not(BigInt binary_number)
{
    for (int64_t i = (binary_number.number.size()-1); i >= 0; --i)
    {
        int32_t factor = 1;
        for (long long int j = 0; j < SIZE_NUM_DIGIT; ++j)
        {
            if ((binary_number.number[i] % (factor * 10)) / factor == 1)
            {
                binary_number.number[i] -= factor;
            }
            else
            {
                binary_number.number[i] += factor;
            }
            factor *= 10;
        }
    }
    return binary_number;
}

BigInt::BigInt()
{
    number.push_back(0);
}

BigInt::BigInt(int32_t input_number)
{
    if (input_number < 0)
    {
        sign = -1;
        input_number *= sign;
    }
    if (input_number == 0)
    {
        number.push_back(0);
    }
    while (input_number > 0)
    {
        number.push_back(input_number % MAX_NUM);
        input_number /= MAX_NUM;
    }
}

BigInt::BigInt(std::string input_number)
{
    std::string_view input_number_view{input_number};
    if (input_number_view.empty())
    {
        throw std::invalid_argument("Empty string.");
    }
    if (input_number_view[0] == '-')
    {
        if ((input_number_view.length() == 2) && ((input_number_view[1]) == '0'))
        {
            sign = 1;
        }
        else
        {
            sign = -1;
        }
        input_number_view.remove_prefix(1);
    }
    else if (input_number_view[0] == '+')
    {
        input_number_view.remove_prefix(1);
    }
    const std::string valid_characters = "1234567890";
    if (input_number_view.find_first_not_of(valid_characters) != input_number_view.npos)
    {
        throw std::invalid_argument("Invalid syntax.");
    }
    size_t str_len = input_number_view.length();
    while (str_len > 0)
    {
        if (str_len >= 9)
        {
            int tmp = 0;
            str_len -= 9;
            std::from_chars(input_number_view.data() + str_len, input_number_view.data() + str_len + 9, tmp);
            number.push_back(tmp);
        }
        else
        {
            int tmp = 0;
            std::from_chars(input_number_view.data(), input_number_view.data() + str_len, tmp);
            number.push_back(tmp);
            break;
        }
    }
    size_t size_vector = number.size();
    size_t real_size = size_vector;
    size_vector--;
    while ((number[size_vector] == 0) && (real_size > 1))
    {
        real_size = size_vector;
        size_vector--;
    }
    number.resize(real_size);
}

BigInt::BigInt(const BigInt &other)
{
    number = other.number;
    sign = other.sign;
}

BigInt &BigInt::operator=(const BigInt &other)
{
    number = other.number;
    sign = other.sign;
    return *this;
}

bool BigInt::operator==(const BigInt &other) const
{
    return (number == other.number) && (sign == other.sign);
}

bool BigInt::operator!=(const BigInt &other) const
{
    return !(*this == other);
}

bool BigInt::operator<(const BigInt &other) const
{
    if (sign != other.sign)
    {
        return sign < other.sign;
    }
    if (number.size() != other.number.size())
    {
        return number.size() * sign < other.number.size() * sign;
    }
    for (int64_t i = static_cast<int64_t> (number.size() - 1); i >= 0; --i)
    {
        if (number[i] != other.number[i])
        {
            return number[i] * sign < other.number[i] * sign;
        }
    }
    return false;
}

bool BigInt::operator>(const BigInt &other) const
{
    return (other < *this);
}

bool BigInt::operator<=(const BigInt &other) const
{
    return !(other < *this);
}

bool BigInt::operator>=(const BigInt &other) const
{
    return !(*this < other);
}

BigInt BigInt::operator+() const
{
    return *this;
}

BigInt BigInt::operator-() const
{
    BigInt tmp = *this;
    tmp.sign *= -1;
    if (*this == BigInt(0))
    {
        tmp.sign = 1;
    }
    return tmp;
}

BigInt &BigInt::operator++()
{
    *this += 1;
    return *this;
}

BigInt &BigInt::operator--()
{
    *this -= 1;
    return *this;
}

const BigInt BigInt::operator++(int32_t)
{
    BigInt tmp = *this;
    *this += 1;
    return tmp;
}

const BigInt BigInt::operator--(int32_t)
{
    BigInt tmp = *this;
    *this -= 1;
    return tmp;
}

BigInt &BigInt::operator+=(const BigInt &other)
{
    if (sign == other.sign)
    {
        if (((*this < other) && (sign == 1)) || ((*this > other) && (sign == -1)))
        {
            size_t size_vector_max = std::max(number.size(), other.number.size());
            number.resize(size_vector_max, 0);
            size_vector_max--;
            for (size_t i = 0; i < size_vector_max; ++i)
            {
                number[i] += other.number[i];
                if (number[i] / MAX_NUM != 0)
                {
                    number[i] %= MAX_NUM;
                    number[i + 1]++;
                }
            }
            number[size_vector_max] += other.number[size_vector_max];
            if (number[size_vector_max] / MAX_NUM != 0)
            {
                number[size_vector_max] %= MAX_NUM;
                number.push_back(1);
            }
        }
        else
        {
            size_t size_vector = other.number.size();
            size_vector--;
            for (size_t i = 0; i < size_vector; ++i)
            {
                number[i] += other.number[i];
                if (number[i] / MAX_NUM != 0)
                {
                    number[i] %= MAX_NUM;
                    number[i + 1]++;
                }
            }
            number[size_vector] += other.number[size_vector];
            if (number[size_vector] / MAX_NUM != 0)
            {
                number[size_vector] %= MAX_NUM;
                number.push_back(1);
            }
        }
    }
    else
    {
        *this -= -other;
    }
    size_t real_size = number.size();
    size_t size_vector = real_size - 1;
    number.resize(real_size);
    if (-*this == BigInt(0))
    {
        sign = 1;
    }
    return *this;
}

BigInt &BigInt::operator-=(const BigInt &other)
{
    if (sign == other.sign)
    {
        if (((*this < other) && (sign == 1)) || ((*this > other) && (sign == -1)))
        {
            size_t size_vector = other.number.size();
            number.resize(size_vector, 0);
            for (size_t i = 0; i < size_vector; ++i)
            {
                number[i] = other.number[i] - number[i];
                if (number[i] < 0)
                {
                    number[i] += MAX_NUM;
                    number[i + 1]++;
                }
            }
            sign *= -1;
        }
        else
        {
            size_t size_vector = other.number.size();
            for (size_t i = 0; i < size_vector; ++i)
            {
                number[i] -= other.number[i];
                if (number[i] < 0)
                {
                    number[i] += MAX_NUM;
                    number[i + 1]--;
                }
            }
        }
    }
    else
    {
        *this += -other;
    }

    size_t real_size = number.size();
    size_t size_vector = real_size - 1;
    while ((number[size_vector] == 0) && (real_size > 1))
    {
        real_size = size_vector;
        size_vector--;
    }
    number.resize(real_size);
    if (-*this == BigInt(0))
    {
        sign = 1;
    }
    return *this;
}

BigInt &BigInt::operator*=(const BigInt &other)
{
    size_t this_size_vector = number.size();
    size_t other_size_vector = static_cast<int64_t> (other.number.size());
    size_t size_vector = this_size_vector + other_size_vector;
    std::vector<int64_t> tmp;
    tmp.resize(size_vector, 0);
    number.resize(size_vector, 0);
    for (size_t i = 0; i < this_size_vector; ++i)
    {
        for (size_t j = 0; j < other_size_vector; ++j)
        {
            tmp[i + j] += static_cast<int64_t> (number[i]) * static_cast<int64_t> (other.number[j]);
            if (tmp[i + j] / MAX_NUM != 0)
            {
                tmp[i + j + 1] += tmp[i + j] / MAX_NUM;
                tmp[i + j] %= MAX_NUM;
            }
        }
    }
    for (size_t i = 0; i < size_vector; ++i)
    {
        number[i] = static_cast<int32_t> (tmp[i]);
    }
    size_t real_size = size_vector;
    size_vector--;
    while ((number[size_vector] == 0) && (real_size > 1))
    {
        real_size = size_vector;
        size_vector--;
    }
    number.resize(real_size);
    sign *= other.sign;
    if (-*this == BigInt(0))
    {
        sign = 1;
    }
    return *this;
}

BigInt &BigInt::operator/=(const BigInt &other)
{
    if (other == BigInt(0))
    {
        throw std::invalid_argument("Division by zero.");
    }
    size_t this_vector_size = number.size();
    int8_t sign_tmp = sign * other.sign;
    BigInt result(0), current(0), other_copy(other);
    other_copy.sign = 1;
    result.number.resize(this_vector_size);
    for (int64_t i = (this_vector_size - 1); i >= 0; --i)
    {
        current *= MAX_NUM;
        current += number[i];
        int32_t left = 0, right = MAX_NUM, mid = 0, digit = 0;
        while (left <= right)
        {
            mid = (left + right) >> 1;
            BigInt tmp(static_cast<int32_t> (mid));
            tmp *= other_copy;
            if (tmp < current)
            {
                digit = mid;
                left = mid + 1;
            }
            else if (tmp > current)
            {
                right = mid - 1;
            }
            else
            {
                digit = mid;
                break;
            }

        }
        result.number[i] = digit;
        BigInt tmp(digit);
        tmp *= other_copy;
        current -= tmp;
    }

    size_t real_size = this_vector_size;
    this_vector_size--;
    while ((result.number[this_vector_size] == 0) && (real_size > 1))
    {
        real_size = this_vector_size;
        this_vector_size--;
    }
    result.sign = sign * other.sign;
    result.number.resize(real_size);
    *this = result;
    if (-*this == BigInt(0))
    {
        sign = 1;
    }
    return *this;
}

BigInt &BigInt::operator%=(const BigInt &other)
{
    if (other == BigInt(0))
    {
        throw std::invalid_argument("Division by zero.");
    }
    int64_t this_vector_size = static_cast<int64_t> (number.size());
    BigInt current(0), other_copy(other);
    other_copy.sign = 1;
    for (int64_t i = (this_vector_size - 1); i >= 0; --i)
    {
        current *= MAX_NUM;
        current += number[i];
        int32_t left = 0, right = MAX_NUM, mid = 0, digit = 0;
        while (left <= right)
        {
            mid = (left + right) >> 1;
            BigInt tmp(static_cast<int32_t> (mid));
            tmp *= other_copy;
            if (tmp < current)
            {
                digit = mid;
                left = mid + 1;
            }
            else if (tmp > current)
            {
                right = mid - 1;
            }
            else
            {
                digit = mid;
                break;
            }

        }
        BigInt tmp(digit);
        tmp *= other_copy;
        current -= tmp;
    }

    if ((sign == -1) && (current != BigInt(0)))
    {
        current -= other_copy;
        current.sign = 1;
    }
    *this = current;
    return *this;
}

BigInt &BigInt::operator^=(const BigInt &other)
{
    BigInt this_bin(0);
    BigInt other_bin(0);
    BigInt result(0);
    this_bin = convert_dec_to_bin_mod(*this);
    other_bin = convert_dec_to_bin_mod(other);
    size_t this_size = this_bin.number.size();
    size_t other_size = other_bin.number.size();
    if (sign == -1)
    {
        this_bin = Not(this_bin);
        this_bin = plus_one(this_bin);
    }
    if (other.sign == -1)
    {
        other_bin = Not(other_bin);
        other_bin = plus_one(other_bin);
    }
    if (this_size < other_size)
    {
        if (sign == -1)
        {
            this_bin.number.resize(other_size, 111111111);
        }
        else
        {
            this_bin.number.resize(other_size, 0);
        }
    }
    else if (this_size > other_size)
    {
        if (other.sign == -1)
        {
            other_bin.number.resize(this_size, 111111111);
        }
        else
        {
            other_bin.number.resize(this_size, 0);
        }
    }
    size_t max_size = std::max(this_size, other_size);
    result.number.resize(max_size, 0);
    for (size_t i = 0; i < max_size; ++i)
    {
        int32_t factor = 1;
        for (size_t j = 0; j < SIZE_NUM_DIGIT; ++j)
        {
            if (((this_bin.number[i] % (factor * 10)) / factor) != ((other_bin.number[i] % (factor * 10)) / factor))
                result.number[i] += factor;
            factor *= 10;
        }
    }
    if (sign != other.sign)
    {
        result = Not(result);
        result = plus_one(result);
    }
    result = convert_bin_to_dec_mod(result);
    if (sign != other.sign)
    {
        result.sign = -1;
    }
    *this = result;
    return *this;
}

BigInt &BigInt::operator&=(const BigInt &other)
{
    BigInt this_bin(0);
    BigInt other_bin(0);
    BigInt result(0);
    this_bin = convert_dec_to_bin_mod(*this);
    other_bin = convert_dec_to_bin_mod(other);
    size_t this_size = this_bin.number.size();
    size_t other_size = other_bin.number.size();
    if (sign == -1)
    {
        this_bin = Not(this_bin);
        this_bin = plus_one(this_bin);
    }
    if (other.sign == -1)
    {
        other_bin = Not(other_bin);
        other_bin = plus_one(other_bin);
    }
    if (this_size < other_size)
    {
        if (sign == -1)
        {
            this_bin.number.resize(other_size, 111111111);
        }
        else
        {
            this_bin.number.resize(other_size, 0);
        }
    }
    else if (this_size > other_size)
    {
        if (other.sign == -1)
        {
            other_bin.number.resize(this_size, 111111111);
        }
        else
        {
            other_bin.number.resize(this_size, 0);
        }
    }
    size_t max_size = std::max(this_size, other_size);
    result.number.resize(max_size, 0);
    for (size_t i = 0; i < max_size; ++i)
    {
        int32_t factor = 1;
        for (size_t j = 0; j < SIZE_NUM_DIGIT; ++j)
        {
            if (((this_bin.number[i] % (factor * 10)) / factor == 1) &&
                (other_bin.number[i] % (factor * 10)) / factor == 1)
                result.number[i] += factor;
            factor *= 10;
        }
    }
    if ((sign == -1) && (other.sign == -1))
    {
        result = Not(result);
        result = plus_one(result);
    }
    result = convert_bin_to_dec_mod(result);
    if ((sign == -1) && (other.sign == -1))
    {
        result.sign = -1;
    }
    *this = result;
    return *this;
}

BigInt &BigInt::operator|=(const BigInt &other)
{
    BigInt this_bin(0);
    BigInt other_bin(0);
    BigInt result(0);
    this_bin = convert_dec_to_bin_mod(*this);
    other_bin = convert_dec_to_bin_mod(other);
    size_t this_size = this_bin.number.size();
    size_t other_size = other_bin.number.size();
    if (sign == -1)
    {
        this_bin = Not(this_bin);
        this_bin = plus_one(this_bin);
    }
    if (other.sign == -1)
    {
        other_bin = Not(other_bin);
        other_bin = plus_one(other_bin);
    }
    if (this_size < other_size)
    {
        if (sign == -1)
        {
            this_bin.number.resize(other_size, 111111111);
        }
        else
        {
            this_bin.number.resize(other_size, 0);
        }
    }
    else if (this_size > other_size)
    {
        if (other.sign == -1)
        {
            other_bin.number.resize(this_size, 111111111);
        }
        else
        {
            other_bin.number.resize(this_size, 0);
        }
    }
    size_t max_size = std::max(this_size, other_size);
    result.number.resize(max_size, 0);
    for (size_t i = 0; i < max_size; ++i)
    {
        int32_t factor = 1;
        for (size_t j = 0; j < SIZE_NUM_DIGIT; ++j)
        {
            if (((this_bin.number[i] % (factor * 10)) / factor == 1) ||
                (other_bin.number[i] % (factor * 10)) / factor == 1)
                result.number[i] += factor;
            factor *= 10;
        }
    }
    if ((sign == -1) || (other.sign == -1))
    {
        result = Not(result);
        result = plus_one(result);
    }
    result = convert_bin_to_dec_mod(result);
    if ((sign == -1) || (other.sign == -1))
    {
        result.sign = -1;
    }
    *this = result;
    return *this;
}

BigInt BigInt::operator~() const
{
    BigInt bin_number(0);
    bin_number = convert_dec_to_bin_mod(*this);
    if (sign == -1)
    {
        bin_number = Not(bin_number);
        bin_number = plus_one(bin_number);
    }
    bin_number = Not(bin_number);
    if (sign == 1)
    {
        bin_number = Not(bin_number);
        bin_number = plus_one(bin_number);
    }
    bin_number = convert_bin_to_dec_mod(bin_number);
    bin_number.sign *= sign * (-1);
    return bin_number;
}

BigInt::operator int() const
{
    return sign * number[0];
}

BigInt::operator std::string() const
{
    //std::string result;
    size_t size = number.size();
    size_t len = size * 9 + 1;
    size_t iter = 0;
    std::string result(len, '0');
    if (sign == -1)
    {
        result[iter] = '-';
        iter++;
    }
    int32_t tmp = MAX_NUM;
    bool null = true;
    for (int32_t j = 0; j < SIZE_NUM_DIGIT; ++j)
    {
        result[iter] = ((number[size-1] % tmp) / (tmp / 10)) + 48;
        if ((result[iter] == '0') && (null))
        {
            tmp /= 10;
            continue;
        }
        else
        {
            null = false;
        }
        iter++;
        tmp /= 10;
    }
    for (int64_t i = (size - 2); i >= 0; i--)
    {
        int32_t tmp_1 = MAX_NUM;
        for (int32_t j = 0; j < SIZE_NUM_DIGIT; ++j)
        {
            result[iter] = ((number[i] % tmp_1) / (tmp_1 / 10)) + 48;
            iter++;
            tmp_1 /= 10;
        }
    }
    if (iter != 0)
    {
        result.resize(iter);
    }
    else
    {
        result.resize(1);
    }
    return result;
}

size_t BigInt::size() const
{
    return number.size() * 4;
}

BigInt operator+(const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp += b;
    return tmp;
}

BigInt operator-(const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp -= b;
    return tmp;
}

BigInt operator*(const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp *= b;
    return tmp;
}

BigInt operator/(const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp /= b;
    return tmp;
}

BigInt operator^(const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp ^= b;
    return tmp;
}

BigInt operator%(const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp %= b;
    return tmp;
}

BigInt operator&(const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp &= b;
    return tmp;
}

BigInt operator|(const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp |= b;
    return tmp;
}

std::ostream &operator<<(std::ostream &output, const BigInt &number)
{
    output << std::string(number);
    return output;
}

BigInt pow_bigint(const BigInt &base,const BigInt &power)
{
    BigInt result = base;
    if (power == BigInt(1))
    {
        return base;
    }
    BigInt tmp = result;
    BigInt i = 2;
    BigInt sqr_i = 4;
    result *= result;
    while (sqr_i <= power)
    {
        result *= result;
        i += i;
        sqr_i = 1;
        sqr_i *= i;
        sqr_i *= i;
    }
    for (; i < power; ++i)
    {
        result *= tmp;
    }
    return result;
}
