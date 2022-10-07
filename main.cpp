#include <iostream>
#include <vector>

class BigInt;
const int32_t size_num_digit = 9;
const int32_t max_num_digit = 1000000000;

BigInt operator+(const BigInt&, const BigInt&);
BigInt operator-(const BigInt&, const BigInt&);
BigInt operator*(const BigInt&, const BigInt&);
BigInt operator/(const BigInt&, const BigInt&);
BigInt operator^(const BigInt&, const BigInt&);
BigInt operator%(const BigInt&, const BigInt&);

class BigInt
{
    friend BigInt convert_dec_to_bin_mod(BigInt decimal_number);
    friend BigInt convert_bin_to_dec_mod(BigInt binary_number);
    friend BigInt plus_one(BigInt binary_number);
    friend BigInt Not(BigInt binary_number);
    friend std::ostream& operator<<(std::ostream& o, const BigInt& i);
    friend BigInt operator & (const BigInt &a, const BigInt &b);
    friend BigInt operator | (const BigInt &a, const BigInt &b);
private:
    std :: vector <int32_t> number;
    char sign = 1;
public:
    BigInt ()
    {
        number.push_back(0);
    }

    BigInt (int32_t input_number)
    {
        if (input_number < 0)
        {
            sign = -1;
            input_number *= sign;
        }
        if (input_number == 0)
            number.push_back(0);
        while (input_number > 0)
        {
            number.push_back(input_number % max_num_digit);
            input_number /= max_num_digit;
        }
    }

    BigInt (std :: string input_number)
    {
        std :: string valid_characters = "-1234567890";
        if (input_number.find_first_not_of(valid_characters) != input_number.npos)
        {
            throw std::invalid_argument("Invalid syntax.");
        }
        if (input_number[0] == '-')
        {
            sign = -1;
            input_number = input_number.substr(1).c_str();
        }
        unsigned long long int str_len = input_number.length();
        while (str_len > 0)
        {
            if (str_len >= 9)
            {
                str_len -= 9;
                number.push_back(atoi(input_number.substr(str_len, 9).c_str()));
            }
            else
            {
                number.push_back(atoi(input_number.substr(0, str_len).c_str()));
                break;
            }
        }
    }

    BigInt (const BigInt& other)
    {
        this->number = other.number;
        this->sign = other.sign;
    }

    BigInt & operator = (const BigInt & other)
    {
        this->number = other.number;
        this->sign = other.sign;
        return *this;
    }

    bool operator == (const BigInt & other) const
    {
        return (this->number == other.number) && (this->sign == other.sign);
    }

    bool operator != (const BigInt & other) const
    {
        return !(*this == other);
    }

    bool operator < (const BigInt & other) const
    {
        if (this->sign != other.sign)
            return this->sign < other.sign;
        if (this->number.size() != other.number.size())
            return this->number.size()*this->sign < other.number.size() * this->sign;
        for (int32_t i = ((int32_t)this->number.size()- 1); i >= 0; --i)
            if (this->number[i] != other.number[i])
                return this->number[i]*this->sign < other.number[i]* this->sign;
        return false;
    }

    bool operator > (const BigInt & other) const
    {
        return (other < *this);
    }

    bool operator <= (const BigInt & other) const
    {
        return !(other < *this);
    }

    bool operator >= (const BigInt & other) const
    {
        return !(*this < other);
    }

    BigInt operator + () const
    {
        BigInt tmp = *this;
        tmp.sign = 1;
        return tmp;
    }

    BigInt operator - () const
    {
        BigInt tmp = *this;
        tmp.sign *= -1;
        return tmp;
    }

    BigInt& operator++()
    {
        *this += 1;
        return  *this;
    }

    BigInt& operator--()
    {
        *this -= 1;
        return  *this;
    }

    const BigInt operator ++ (int32_t)
    {
        BigInt tmp = *this;
        *this += 1;
        return tmp;
    }

    const BigInt operator -- (int32_t)
    {
        BigInt tmp = *this;
        *this -= 1;
        return tmp;
    }

    BigInt & operator += (const BigInt & other)
    {
        if (this->sign == other.sign)
        {
            if ((*this < other) && (this->sign == 1) || (*this > other) && (this->sign == -1))
            {
                unsigned long long int size_vector_max = std::max(this->number.size(), other.number.size());
                this->number.resize(size_vector_max, 0);
                size_vector_max--;

                for (int32_t i = 0; i < size_vector_max; ++i)
                {
                    this->number[i] += other.number[i];
                    if (this->number[i] / max_num_digit != 0)
                    {
                        this->number[i] %= max_num_digit;
                        this->number[i + 1]++;
                    }
                }
                this->number[size_vector_max] += other.number[size_vector_max];
                if (this->number[size_vector_max] / max_num_digit != 0)
                {
                    this->number[size_vector_max] %= max_num_digit;
                    this->number.push_back(1);
                }
            }
            else
            {
                unsigned long long int size_vector = other.number.size();
                size_vector--;

                for (int32_t i = 0; i < size_vector; ++i)
                {
                    this->number[i] += other.number[i];
                    if (this->number[i] / max_num_digit != 0)
                    {
                        this->number[i] %= max_num_digit;
                        this->number[i + 1]++;
                    }
                }
                this->number[size_vector] += other.number[size_vector];
                if (this->number[size_vector] / max_num_digit != 0)
                {
                    this->number[size_vector] %= max_num_digit;
                    this->number.push_back(1);
                }
            }
        }
        else
        {
            *this -= -other;
        }
        return *this;
    }

    BigInt & operator -=(const BigInt & other)
    {
        if (this->sign == other.sign)
        {
            if ((*this < other) && (this->sign == 1) || (*this > other) && (this->sign == -1))
            {
                unsigned long long int size_vector_max = other.number.size();
                this->number.resize(size_vector_max, 0);

                for (int32_t i = 0; i < size_vector_max; ++i)
                {
                    this->number[i] = other.number[i] - this->number[i];
                    if (this->number[i] < 0)
                    {
                        this->number[i] += max_num_digit;
                        this->number[i + 1]--;
                    }
                }
                this->sign *= -1;
            }
            else
            {
                unsigned long long int size_vector = other.number.size();
                for (int32_t i = 0; i < size_vector; ++i)
                {
                    this->number[i] -= other.number[i];
                    if (this->number[i] < 0)
                    {
                        this->number[i] += max_num_digit;
                        this->number[i + 1]--;
                    }
                }
            }
        }
        else
        {
            *this += -other;
        }
        return *this;
    }

    BigInt & operator *= (const BigInt & other)
    {
        unsigned long long int this_size_vector = this->number.size(), other_size_vector = other.number.size();
        unsigned long long int size_vector = this_size_vector + other_size_vector;
        std :: vector <unsigned long long int> tmp;
        tmp.resize(size_vector, 0);
        this->number.resize(size_vector, 0);
        for (int32_t i = 0; i < this_size_vector; ++i)
        {
            for (int32_t j = 0; j < other_size_vector; ++j)
            {
                tmp[i + j] += (unsigned long long int)this->number[i] * (unsigned long long int)other.number[j];
                if (tmp[i+j] / max_num_digit != 0)
                {
                    tmp[i + j + 1] += tmp[i + j] / max_num_digit;
                    tmp[i + j] %= max_num_digit;
                }
            }
        }
        for (int32_t i = 0; i < size_vector; ++i)
        {
            this->number[i] = (int32_t)tmp[i];
        }
        unsigned long long int real_size = size_vector;
        size_vector--;
        while ((number[size_vector] == 0) && (real_size > 1))
        {
            real_size = size_vector;
            size_vector--;
        }
        number.resize(real_size);
        sign *= other.sign;

        return *this;
    }

    BigInt & operator /= (const BigInt &other)
    {
        unsigned long long int this_vector_size = number.size();
        char sign_tmp = sign * other.sign;
        BigInt result(0), current(0), other_copy(other);
        other_copy.sign = 1;
        result.number.resize(this_vector_size);
        for (int32_t i = (this_vector_size-1); i >= 0; --i)
        {
            current *= max_num_digit;
            current += number[i];
            unsigned long long int left = 0, right = max_num_digit, mid = 0, digit = 0;
            while (left <= right)
            {
                mid = (left + right) >> 1;
                BigInt tmp((int32_t)mid);
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
            BigInt tmp((int32_t)digit);
            tmp *= other_copy;
            current -= tmp;
        }

        unsigned long long int real_size = this_vector_size;
        this_vector_size--;
        while ((result.number[this_vector_size] == 0) && (real_size > 1))
        {
            real_size = this_vector_size;
            this_vector_size--;
        }
        result.sign = sign * other.sign;
        result.number.resize(real_size);
        *this = result;
        return *this;
    }

    BigInt & operator %= (const BigInt &other)
    {
        unsigned long long int this_vector_size = number.size();
        BigInt current(0), other_copy(other);
        other_copy.sign = 1;
        for (int32_t i = (this_vector_size-1); i >= 0; --i)
        {
            current *= max_num_digit;
            current += number[i];
            unsigned long long int left = 0, right = max_num_digit, mid = 0, digit = 0;
            while (left <= right)
            {
                mid = (left + right) >> 1;
                BigInt tmp((int32_t)mid);
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
            BigInt tmp((int32_t)digit);
            tmp *= other_copy;
            current -= tmp;
        }

        unsigned long long int real_size = current.number.size();
        real_size--;
        while ((current.number[real_size] == 0) && (real_size > 0))
        {
            real_size --;
        }
        real_size++;
        current.number.resize(real_size);
        if ((sign == -1) && (current != BigInt(0)))
        {
            current -= other_copy;
            current.sign = 1;
        }
        *this = current;
        return *this;
    }

    BigInt& operator ^= (const BigInt &other)
    {
        if (other == BigInt(1))
            return *this;

        BigInt base = *this, i = 2, sqr_i = 4;
        *this *= *this;

        while (sqr_i <= other)
        {
            *this *= *this;
            i += i;
            sqr_i = 1;
            sqr_i *= i;
            sqr_i *= i;
        }
        for (; i < other; ++i)
            *this *= base;
        return *this;
    }

    BigInt operator~() const
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

    operator int() const
    {
        return number[0];
    }

    operator std::string() const
    {
        std :: string result;
        unsigned long long int size = number.size();
        if (sign == -1)
            result.push_back('-');
        result += std :: to_string(number[size - 1]);
        for (int32_t i = (size-2); i >= 0; i--)
        {
            int32_t tmp = max_num_digit;
            for (int32_t j = 0; j < size_num_digit; ++j)
            {
                result += std :: to_string((number[i] % tmp)/(tmp / 10));
                tmp /= 10;
            }
        }
        return result;
    }

    size_t size() const
    {
        return this->number.size() * 4;
    }

    ~BigInt()
    {
        std :: vector <int32_t> ().swap(this->number);
    }
};

BigInt convert_dec_to_bin_mod(BigInt decimal_number)
{
    BigInt binary_number(0), factor(1);
    decimal_number.sign = 1;
    while (decimal_number > BigInt(0))
    {
        binary_number += (decimal_number % BigInt(2)) * factor;
        factor *= 10;
        decimal_number /= 2;
    }
    return binary_number;
}

BigInt convert_bin_to_dec_mod(BigInt binary_number)
{
    BigInt decimal_number(0), factor(1);
    binary_number.sign = 1;
    while (binary_number > BigInt(0))
    {
        decimal_number += (binary_number % BigInt(2)) * factor;
        factor *= 2;
        binary_number /= 10;
    }
    return decimal_number;
}

BigInt plus_one(BigInt binary_number)
{
    binary_number = convert_bin_to_dec_mod(binary_number);
    binary_number++;
    binary_number = convert_dec_to_bin_mod(binary_number);
    return binary_number;
}

BigInt Not(BigInt binary_number)
{
    unsigned long long int vector_size = binary_number.number.size();
    vector_size--;
    for (long long int i = vector_size; i >= 0; --i)
    {
        unsigned long long int factor = 1;
        for (long long int j = 0; j < size_num_digit; ++j)
        {
            if ((binary_number.number[i] % (factor * 10)) / factor == 1)
                binary_number.number[i] -= factor;
            else
                binary_number.number[i] += factor;
            factor *= 10;
        }
    }
    return binary_number;
}

BigInt operator + (const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp += b;
    return tmp;
}

BigInt operator - (const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp -= b;
    return tmp;
}

BigInt operator * (const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp *= b;
    return tmp;
}

BigInt operator / (const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp /= b;
    return tmp;
}

BigInt operator ^ (const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp ^= b;
    return tmp;
}

BigInt operator % (const BigInt &a, const BigInt &b)
{
    BigInt tmp = a;
    tmp %= b;
    return tmp;
}

BigInt operator & (const BigInt &a, const BigInt &b)
{
    BigInt a_bin(0), b_bin(0), result(0);
    a_bin = convert_dec_to_bin_mod(a);
    b_bin = convert_dec_to_bin_mod(b);
    int32_t a_size = a_bin.number.size(), b_size = b_bin.number.size();
    if (a.sign == -1)
    {
        a_bin = Not(a_bin);
        a_bin = plus_one(a_bin);
    }
    if (b.sign == -1)
    {
        b_bin = Not(b_bin);
        b_bin = plus_one(b_bin);
    }
    if (a_size < b_size)
    {
        if (a.sign == -1)
            a_bin.number.resize(b_size, 111111111);
        else
            a_bin.number.resize(b_size, 0);
    }
    else if (a_size > b_size)
    {
        if (b.sign == -1)
            b_bin.number.resize(a_size, 111111111);
        else
            b_bin.number.resize(a_size, 0);
    }
    unsigned long long int max_size = std::max(a_size, b_size);
    result.number.resize(max_size, 0);
    for (long long int i = 0; i < max_size; ++i)
    {
        int32_t factor = 1;
        for (long long int j = 0; j < size_num_digit; ++j)
        {
            if (((a_bin.number[i] % (factor * 10)) / factor == 1) && (b_bin.number[i] % (factor * 10)) / factor == 1)
                result.number[i] += factor;
            factor *= 10;
        }
    }
    if ((a.sign == -1) && (b.sign == -1))
    {
        result = Not(result);
        result = plus_one(result);
    }
    result = convert_bin_to_dec_mod(result);
    if ((a.sign == -1) && (b.sign == -1))
        result.sign = -1;
    return result;
}

BigInt operator | (const BigInt &a, const BigInt &b)
{
    BigInt a_bin(0), b_bin(0), result(0);
    a_bin = convert_dec_to_bin_mod(a);
    b_bin = convert_dec_to_bin_mod(b);
    int32_t a_size = a_bin.number.size(), b_size = b_bin.number.size();
    if (a.sign == -1)
    {
        a_bin = Not(a_bin);
        a_bin = plus_one(a_bin);
    }
    if (b.sign == -1)
    {
        b_bin = Not(b_bin);
        b_bin = plus_one(b_bin);
    }
    if (a_size < b_size)
    {
        if (a.sign == -1)
            a_bin.number.resize(b_size, 111111111);
        else
            a_bin.number.resize(b_size, 0);
    }
    else if (a_size > b_size)
    {
        if (b.sign == -1)
            b_bin.number.resize(a_size, 111111111);
        else
            b_bin.number.resize(a_size, 0);
    }
    unsigned long long int max_size = std::max(a_size, b_size);
    result.number.resize(max_size, 0);
    for (long long int i = 0; i < max_size; ++i)
    {
        int32_t factor = 1;
        for (long long int j = 0; j < size_num_digit; ++j)
        {
            if (((a_bin.number[i] % (factor * 10)) / factor == 1) || (b_bin.number[i] % (factor * 10)) / factor == 1)
                result.number[i] += factor;
            factor *= 10;
        }
    }
    if ((a.sign == -1) || (b.sign == -1))
    {
        result = Not(result);
        result = plus_one(result);
    }
    result = convert_bin_to_dec_mod(result);
    if ((a.sign == -1) || (b.sign == -1))
        result.sign = -1;
    return result;
}

std::ostream& operator<<(std::ostream& output, const BigInt& number)
{
    unsigned long long int size = number.number.size();
    if (number.sign == -1)
        output << "-";
    std :: cout << number.number[size - 1];
    for (int32_t i = (size-2); i >= 0; i--)
    {
        int32_t tmp = max_num_digit;
        for (int32_t j = 0; j < size_num_digit; ++j)
        {
            output << (number.number[i] % tmp)/(tmp / 10);
            tmp /= 10;
        }
    }
    output << std :: endl;
    return output;
}

int32_t main()
{
    BigInt a("7656532");
    BigInt b("-23254564");
    BigInt c;
    c = a | b;
    std :: cout << c;
    return 0;
}
