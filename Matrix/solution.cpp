#include <iostream>
#include <vector>
#include <string>

class BigInteger{
public:
  std::vector<int> digits;
  int sign;
  BigInteger();
  
  BigInteger(int number) {
    if (number > 0) sign = 1;
    if (number == 0) sign = 0;
    if (number < 0) {
      sign = -1;
      number *= -1;
    }
    while (number != 0) {
      digits.push_back(number % 10);
      number /= 10;
    }
  }
  
  BigInteger(long long number) {
    if (number > 0) sign = 1;
    if (number == 0) sign = 0;
    if (number < 0) {
      sign = -1;
      number *= -1;
    }
    while (number != 0) {
      digits.push_back(number % 10);
      number /= 10;
    }
  }
  
  BigInteger(size_t number) {
    if (number > 0) sign = 1;
    if (number == 0) sign = 0;
    while (number != 0) {
      digits.push_back(number % 10);
      number /= 10;
    }
  }
  
  BigInteger(unsigned long long number) {
    if (number > 0) sign = 1;
    if (number == 0) sign = 0;
    while (number != 0) {
      digits.push_back(number % 10);
      number /= 10;
    }
  }
  
  BigInteger(const BigInteger& number): digits(number.digits), sign(number.sign) {
    remove_non_sign_zeros();
  }
  
  BigInteger& operator+=(const BigInteger&);
  BigInteger& operator-=(const BigInteger&);
  BigInteger& operator*=(const BigInteger&);
  BigInteger& operator*=(int);
  BigInteger& operator/=(const BigInteger&);
  BigInteger& operator%=(const BigInteger&);
  
  
  BigInteger& operator=(BigInteger) ;
  BigInteger& operator=(int) ;
  explicit operator bool();
  explicit operator int();
  
  friend bool operator<(const BigInteger &, const BigInteger &);
  friend bool operator>(const BigInteger &, const BigInteger &);
  friend bool operator==(const BigInteger &, const BigInteger &);
  friend bool operator<=(const BigInteger &, const BigInteger &);
  friend bool operator>=(const BigInteger &, const BigInteger &);
  friend bool operator!=(const BigInteger &, const BigInteger &);
  
  friend std::ostream& operator<<(std::ostream &, const BigInteger&);
  friend std::istream& operator>>(std::istream&, BigInteger&);
  
  BigInteger operator-();
  BigInteger& operator++(); //prefix
  BigInteger operator++(int);
  BigInteger& operator--();
  BigInteger operator--(int);
  void bitwise_subtraction_less_from_more(const BigInteger&);
  void bitwise_subtraction_more_from_less(const BigInteger&);
  void bitwise_add(const BigInteger&);
  
  void set_sign(int x) {
    sign = x;
  }
  void remove_non_sign_zeros();
  const std::string toString() const;
  void right();
  ~BigInteger();
};


BigInteger operator ""_bi(unsigned long long x) {
  return BigInteger(x);
}

BigInteger operator%(const BigInteger& a, const BigInteger& b) {
  BigInteger copy = a;
  copy %= b;
  copy.remove_non_sign_zeros();
  return copy;
}


BigInteger operator/(const BigInteger& a, const BigInteger& b) {
  BigInteger copy = a;
  copy /= b;
  copy.remove_non_sign_zeros();
  return copy;
}

BigInteger operator*(const BigInteger& a, const BigInteger& b) {
  BigInteger copy = a;
  copy *= b;
  copy.remove_non_sign_zeros();
  return copy;
}


BigInteger operator+(const BigInteger& a, const BigInteger& b) {
  BigInteger copy = a;
  copy += b;
  copy.remove_non_sign_zeros();
  return copy;
}

BigInteger operator-(const BigInteger& a, const BigInteger& b) {
  BigInteger copy = a;
  copy -= b;
  copy.remove_non_sign_zeros();
  return copy;
}


BigInteger BigInteger::operator--(int) {
  
  BigInteger copy = *this;
  *this -= 1;
  remove_non_sign_zeros();
  return copy;
}


BigInteger& BigInteger::operator--() {
  *this -= 1;
  remove_non_sign_zeros();
  return *this;
}

BigInteger BigInteger::operator++(int) {
  BigInteger copy = *this;
  *this += 1;
  remove_non_sign_zeros();
  return copy;
}


BigInteger& BigInteger::operator++() {
  *this += 1;
  remove_non_sign_zeros();
  return *this;
}


void div2(BigInteger& a) {
  for(int i = (int)a.digits.size() - 1; i > -1; --i) {
    if (i > 0) {
      a.digits[i - 1] += (a.digits[i] % 2) * 10;
    }
    a.digits[i] /= 2;
  }
  a.remove_non_sign_zeros();
}

int ost2(BigInteger& a) {
  return a.digits[0] % 2;
}

void multiply2(BigInteger& a) {
  int carry = 0;
  for (int i = 0; i < (int)a.digits.size(); ++i) {
    a.digits[i] *= 2;
    a.digits[i] += carry;
    if (a.digits[i] >= 10) {
      a.digits[i] -= 10;
      carry = 1;
    }
    else {
      carry = 0;
    }
  }
  
  if (carry) a.digits.push_back(carry);
}

BigInteger gcd(const BigInteger& a, const BigInteger& b){
  BigInteger k = 1;
  BigInteger a1 = a, b1 = b;
  while ((a1 != 0) && (b1 != 0)) {
    while ((ost2(a1) == 0) && ((ost2(b1) == 0))){
      div2(a1);
      div2(b1);
      multiply2(k);
    }
    while (ost2(a1) == 0)
      div2(a1);
    while (ost2(b1) == 0)
      div2(b1);
    if (a1 >= b1)
      a1 -= b1;
    else
      b1 -= a1;
  }
  return b1 * k;
}

const std::string BigInteger::toString() const{
  std::string ans = "";
  if (sign == -1) ans += '-';
  for (size_t i = digits.size(); i > 0; --i) {
    ans += std::to_string(digits[i - 1]);
  }
  return ans;
}

BigInteger::operator bool(){
  remove_non_sign_zeros();
  if ((digits.size() == 1 && digits[0] == 0) || sign == 0 || digits.size() == 0) return 0;
  return 1;
}

BigInteger::operator int(){
  remove_non_sign_zeros();
  int ans = 0;
  int t = 1;
  for(size_t i = 0; i < digits.size(); ++i){
    ans += digits[i] * t;
    t *= 10;
  }
  return sign * ans;
}

void BigInteger::remove_non_sign_zeros() {
  while (this->digits.size() > 1 && this->digits.back() == 0) {
    this->digits.pop_back();
  }
  if (this->digits.size() == 1 && this->digits[0] == 0)
    this->sign = 0;
}

BigInteger::~BigInteger() {
  digits.clear();
}


BigInteger::BigInteger(){
  std::vector<int> digits;
  std::string num_str;
  sign = 0;
}

std::istream& operator>>(std::istream & is, BigInteger & input_number) {
  std::ios_base::sync_with_stdio(0);
  std::cin.tie(0);
  std::string input_str;
  is >> input_str;
  if (input_str.size() == 0) {
    input_number.sign = 0;
    return is;
  }
  input_number.digits.clear();
  input_number.sign = 1;
  bool fl_0 = 1;
  for (size_t i = input_str.size(); i > 0; --i) {
    if (i == 1 && input_str[0] == '-') {
      input_number.sign = -1;
      continue;
    }
    if (input_str[i] != '0') {
      fl_0 = 0;
    }
    input_number.digits.push_back(input_str[i - 1] - '0');
  }
  if (fl_0) input_number.sign = 0;
  input_number.remove_non_sign_zeros();
  return is;
}

std::ostream& operator<<(std::ostream &os, const BigInteger &number) {
  std::ios_base::sync_with_stdio(0);
  std::cout.tie(0);
  if (number.digits.size() == 0) {
    os << 0;
    return os;
  }
  if (number.sign == -1) {
    os << '-';
  }
  for (size_t i = number.digits.size(); i > 0; --i) {
    os << number.digits[i - 1];
  }
  return os;
}


BigInteger BigInteger::operator -(){
  BigInteger copy = *this;
  copy.sign = -1 * copy.sign;
  remove_non_sign_zeros();
  return copy;
}

void BigInteger::bitwise_subtraction_less_from_more(const BigInteger& number) {
  int carry = 0;
  for (size_t i = 0; i < number.digits.size() || carry; ++i) {
    if (i == digits.size())
      digits.push_back (0);
    if (i < number.digits.size()) {
      digits[i] -= carry + number.digits[i];
    }
    else {
      digits[i] -= carry;
    }
    carry = 0;
    if (digits[i] < 0)
      carry = 1;
    if (carry)
      digits[i] += 10;
  }
}

void BigInteger::bitwise_subtraction_more_from_less(const BigInteger& number) {
  int carry = 0;
  for (size_t i = 0; i < number.digits.size() || carry; ++i) {
    if (i == digits.size())
      digits.push_back (0);
    if (i < number.digits.size()) {
      digits[i] = number.digits[i] - digits[i] - carry;
    }
    else {
      digits[i] = -digits[i] - carry;
    }
    carry = 0;
    if (digits[i] < 0)
      carry = 1;
    if (carry)
      digits[i] += 10;
  }
}
void BigInteger::bitwise_add(const BigInteger& number) {
  int carry = 0;
  for (size_t i = 0; i < std::max(digits.size(), number.digits.size()) || carry; ++i) {
    if (i == digits.size())
      digits.push_back (0);
    if (i < number.digits.size())
      digits[i] += carry + number.digits[i];
    else
      digits[i] += carry;
    carry = 0;
    if (digits[i] >= 10)
      carry = 1;
    if (carry)
      digits[i] -= 10;
  }
}

BigInteger& BigInteger::operator-=(const BigInteger& number){
  if(number.sign == 0)
    return *this;
  if(sign == 0)
    sign = 1;
  if(sign == number.sign) {
    if (sign >= 0) {
      if (*this >= number) {
        bitwise_subtraction_less_from_more(number);
      }
      else {
        sign = -1;
        bitwise_subtraction_more_from_less(number);
      }
    }
    else {
      if (*this * sign > number * number.sign) {
        bitwise_subtraction_less_from_more(number);
      }
      else {
        sign = 1;
        bitwise_subtraction_more_from_less(number);
      }
    }
  }
  else {
    bitwise_add(number);
  }
  this->remove_non_sign_zeros();
  return *this;
}


BigInteger& BigInteger::operator+=(const BigInteger& number){
  if (sign * number.sign != -1) {
    if (number > 0) sign = 1;
    if (number < 0) sign = -1;
    bitwise_add(number);
  }
  else {
    if (*this * sign > number * number.sign) {
      bitwise_subtraction_less_from_more(number);
    }
    else {
      sign *= -1;
      bitwise_subtraction_more_from_less(number);
    }
  }
  this->remove_non_sign_zeros();
  return *this;
}

BigInteger& BigInteger::operator*=(int number) {
  if (number < 0) {
    sign *= -1;
    number *= -1;
  }
  int carry = 0;
  for (size_t i = 0; i < digits.size() || carry; ++i) {
    if (i == digits.size())
      digits.push_back(0);
    long long cur = carry + 1ll * digits[i] * number;
    digits[i] = int (cur % 10);
    carry = int (cur / 10);
  }
  remove_non_sign_zeros();
  return *this;
}

BigInteger& BigInteger::operator*=(const BigInteger& number) {
  sign = sign * number.sign;
  std::vector<int> y;
  for (size_t i = 0; i < digits.size(); ++i){
    int carry = 0;
    for (size_t j = 0; j < number.digits.size() || carry; ++j) {
      int cur;
      if (y.size() <= i + j) {
        y.push_back(0);
      }
      if (j < number.digits.size()) {
        cur = y[i + j] + digits[i] * number.digits[j] + carry;
      }
      else {
        cur = y[i + j] + carry;
      }
      y[i + j] = cur % 10;
      carry = cur / 10;
    }
  }
  digits = y;
  this->remove_non_sign_zeros();
  return *this;
}

void BigInteger::right() {
  if (this->digits.size() == 0) {
    this->digits.push_back(0);
    return;
  }
  digits.push_back(digits[digits.size() - 1]);
  for (size_t i = digits.size() - 2; i > 0; --i) {
    digits[i] = digits[i - 1];
  }
  digits[0] = 0;
}


BigInteger& BigInteger::operator /= (const BigInteger& number) {
  sign = sign * number.sign;
  BigInteger result, current;
  BigInteger copy_s1 = number;
  copy_s1.sign = 1;
  result.digits.resize(digits.size());
  result.sign = sign;
  current.sign = 1;
  for (size_t i = digits.size(); i > 0; --i) {
    current.right();
    current.digits[0] = digits[i - 1];
    current.sign = 1;
    current.remove_non_sign_zeros();
    int new_digit = 0;
    while (new_digit < 10 && copy_s1 * (new_digit + 1) <= current) {
      new_digit++;
    }
    result.digits[i - 1] = new_digit;
    current = current - copy_s1 * new_digit;
  }
  *this = result;
  this->remove_non_sign_zeros();
  return *this;
}


BigInteger& BigInteger::operator%=(const BigInteger& number) {
  *this = *this - number * (*this / number);
  return *this;
}


bool operator==(const BigInteger &st1, const BigInteger &st2) {
  if (st1.sign == st2.sign && st1.sign == 0) return 1;
  if (st1.sign != st2.sign) return 0;
  if (st1.digits.size() != st2.digits.size()) return 0;
  for (size_t i = 0; i < st1.digits.size(); ++i) {
    if (st1.digits[i] != st2.digits[i]) return 0;
  }
  return 1;
}

bool operator != (const BigInteger &st1, const BigInteger &st2) {
  return !(st1 == st2);
}

bool operator < (const BigInteger &st1, const BigInteger &st2) {
  if (st1 == st2) {
    return 0;
  }
  if (st1.sign < st2.sign) {
    return 1;
  }
  if (st1.sign > st2.sign) {
    return 0;
  }
  if (st1.digits.size() < st2.digits.size()) {
    if (st1.sign > 0)
      return 1;
    return 0;
  }
  if (st1.digits.size() > st2.digits.size()) {
    if (st1 > 0)
      return 0;
    return 1;
  }
  for (size_t i = st1.digits.size(); i > 0; --i) {
    if (st1.digits[i - 1] < st2.digits[i - 1]) {
      if (st1.sign == -1) return 0;
      return 1;
    }
    if (st1.digits[i - 1] > st2.digits[i - 1]) {
      if (st1.sign == -1) return 1;
      return 0;
    }
  }
  return 0;
}

bool operator > (const BigInteger &st1, const BigInteger &st2) {
    return st2 < st1;
}

bool operator >= (const BigInteger &st1, const BigInteger &st2) {
    return (st1 > st2) || (st1 == st2);
}

bool operator <= (const BigInteger &st1, const BigInteger &st2) {
    return (st1 < st2) || (st1 == st2);
}


BigInteger& BigInteger::operator=(int number) {
  if (number > 0) sign = 1;
  if (number == 0) sign = 0;
  if (number < 0) sign = -1;
  digits.clear();
  if (number == 0) {
    digits.push_back(0);
  }
  if (number < 0) number *= -1;
  while (number != 0) {
    digits.push_back(number % 10);
    number /= 10;
  }
  remove_non_sign_zeros();
  return *this;
}

BigInteger& BigInteger::operator=(BigInteger number) {
  digits = number.digits;
  sign = number.sign;
  remove_non_sign_zeros();
  return *this;
}

const BigInteger operator / (const BigInteger& a, int b){
  int carry = 0;
  BigInteger copy_a = a;
  for (int i = (int)copy_a.digits.size() - 1; i >= 0; --i) {
    long long cur = copy_a.digits[i] + carry * 1ll * 10;
    copy_a.digits[i] = int (cur / b);
    carry = int (cur % b);
  }
  copy_a.remove_non_sign_zeros();
  return copy_a;
}


//--------------------------------------------------------------------------------------

class Rational{
private:
  BigInteger p;
  BigInteger q;

public:
  Rational();
  Rational(int);
  Rational(const BigInteger&);
  
  Rational& operator+=(const Rational&);
  Rational& operator-=(const Rational&);
  Rational& operator*=(const Rational&);
  Rational& operator/=(const Rational&);
  
  friend std::istream& operator>>(std::istream&, Rational&);

  Rational& operator=(Rational);
  Rational& operator=(int);
  explicit operator double();
  Rational operator-();
  
  friend bool operator<(const Rational &, const Rational &);
  friend bool operator>(const Rational &, const Rational &);
  friend bool operator==(const Rational &, const Rational &);
  friend bool operator<=(const Rational &, const Rational &);
  friend bool operator>=(const Rational &, const Rational &);
  friend bool operator!=(const Rational &, const Rational &);
  
  const std::string toString() const;
};

Rational::Rational(): p(0), q(1)
{
}

Rational::Rational(const BigInteger& t): p(t), q(1)
{
}

Rational::Rational(int x): p(x), q(1)
{
}

const std::string Rational::toString() const{
  std::string ans = "";
  ans += p.toString();
  if (q == 1) {
    return ans;
  }
  ans += '/';
  ans += q.toString();
  return ans;
}

Rational& Rational::operator+=(const Rational& t) {
  p = p * t.q + t.p * q;
  q = q * t.q;
  bool p_is_negative = 0;
  bool q_is_negative = 0;
  if (p < 0) {
    p_is_negative = 1;
    p.set_sign(1);
  }
  if (q < 0) {
    q_is_negative = 1;
    q.set_sign(1);
  }
  BigInteger gcd_p_q = gcd(p, q);
  p /= gcd_p_q;
  q /= gcd_p_q;
  if (p_is_negative) {
    p.set_sign(-1);
  }
  if (q_is_negative) {
    q.set_sign(-1);
  }
  if (p > 0 && q < 0) {
    p.set_sign(-1);
    q.set_sign(1);
  }
  if (p < 0 && q < 0) {
    p.set_sign(1);
    q.set_sign(1);
  }
  if (p == 0) {
    p.set_sign(0);
    q = 1;
  }
  return *this;
}

Rational Rational::operator -(){
  Rational copy = *this;
  copy.p *= -1;
  return copy;
}

Rational& Rational::operator-=(const Rational& t) {
  p = p * t.q - t.p * q;
  q = q * t.q;
  bool p_is_negative = 0;
  bool q_is_negative = 0;
  if (p.sign == -1 && q.sign == -1) {
    p.sign = 1;
    q.sign = 1;
  }
  if (p < 0) {
    p_is_negative = 1;
    p.set_sign(1);
  }
  if (q < 0) {
    q_is_negative = 1;
    q.set_sign(1);
  }
  BigInteger gcd_p_q = gcd(p, q);
  p /= gcd_p_q;
  q /= gcd_p_q;
  if (p_is_negative) {
    p.set_sign(-1);
  }
  if (q_is_negative) {
    q.set_sign(-1);
  }
  if (p > 0 && q < 0) {
    p.set_sign(-1);
    q.set_sign(1);
  }
  if (p == 0) {
    p.set_sign(0);
    q = 1;
  }
  return *this;
}

Rational& Rational::operator*=(const Rational& t){
  p *= t.p;
  q *= t.q;
  bool p_is_negative = 0;
  bool q_is_negative = 0;
  if (p.sign == -1 && q.sign == -1) {
    p.sign = 1;
    q.sign = 1;
  }
  if (p < 0) {
    p_is_negative = 1;
    p.set_sign(1);
  }
  if (q < 0) {
    q_is_negative = 1;
    q.set_sign(1);
  }
  BigInteger gcd_p_q = gcd(p, q);
  p /= gcd_p_q;
  q /= gcd_p_q;
  if (p_is_negative) {
    p.set_sign(-1);
  }
  if (q_is_negative) {
    q.set_sign(-1);
  }
  if (p > 0 && q < 0) {
    p.set_sign(-1);
    q.set_sign(1);
  }
  if (p == 0) {
    p.set_sign(0);
    q = 1;
  }
  return *this;
}

Rational& Rational::operator/=(const Rational& t) {
  p *= t.q;
  q *= t.p;
  bool p_is_negative = 0;
  bool q_is_negative = 0;
  if (p.sign == -1 && q.sign == -1) {
    p.sign = 1;
    q.sign = 1;
  }
  if (p < 0) {
    p_is_negative = 1;
    p.set_sign(1);
  }
  if (q < 0) {
    q_is_negative = 1;
    q.set_sign(1);
  }
  BigInteger gcd_p_q = gcd(p, q);
  p /= gcd_p_q;
  q /= gcd_p_q;
  if (p_is_negative) {
    p.set_sign(-1);
  }
  if (q_is_negative) {
    q.set_sign(-1);
  }
  if (p > 0 && q < 0) {
    p.set_sign(-1);
    q.set_sign(1);
  }
  if (p == 0) {
    p.set_sign(0);
    q = 1;
  }
  return *this;
}

Rational operator+(const Rational& a, const Rational& b) {
  Rational copy = a;
  copy += b;
  return copy;
}

Rational operator-(const Rational& a, const Rational& b) {
  Rational copy = a;
  copy -= b;
  return copy;
}

Rational operator/(const Rational& a, const Rational& b) {
  Rational copy = a;
  copy /= b;
  return copy;
}

Rational operator*(const Rational& a, const Rational& b) {
  Rational copy = a;
  copy *= b;
  return copy;
}

bool operator < (const Rational& a, const Rational& b) {
  if (a.p * b.q < b.p * a.q) {
    return 1;
  }
  return 0;
}

bool operator > (const Rational& a, const Rational& b) {
  return b < a;
}

bool operator == (const Rational& a, const Rational& b) {
  if (a.p == b.p && a.q == b.q) return 1;
  return 0;
}

bool operator <= (const Rational& a, const Rational& b) {
  return (a < b) || (a == b);
}

bool operator >= (const Rational& a, const Rational& b) {
  return (a > b) || (a == b);
}

bool operator != (const Rational& a, const Rational& b) {
  return !(a == b);
}

std::istream& operator>>(std::istream & is, Rational & input_number) {
  std::ios_base::sync_with_stdio(0);
  std::cin.tie(0);
  int num;
  is >> num;
  input_number.p.digits.clear();
  input_number.p.sign = 1;
  input_number.q = 1;
  input_number.p = abs(num);
  if (num < 0) {
    input_number.p.sign = -1;
  }
  input_number.p.remove_non_sign_zeros();
  return is;
}

Rational::operator double(){
  if (p.toString().size() == 0) return 0;
  return stod(p.toString()) / stod(q.toString());
}

Rational& Rational::operator=(Rational number) {
  p = number.p;
  q = number.q;
  return *this;
}

Rational& Rational::operator=(int number) {
  p = number;
  q = 1;
  return *this;
}



template<long long N>
class Finite{

public:
  long long n;
  
  Finite(): n(0)
  {
  }
  Finite(long long y) {
    n = y % N;
    while(n < 0) n += N;
  }
  
  explicit operator long long() {
    return n % N;
  }
  
  Finite(const Finite& finite2): n(finite2.n)
  {
  }
  
  Finite& operator=(long long finite2) {
    n = finite2;
    n %= N;
    while(n < 0) n += N;
    return *this;
  }
  Finite& operator=(const Finite& finite2) {
    n = finite2.n;
    n %= N;
    return *this;
  }
  
  bool operator == (const Finite& finite2) const {
    return n == finite2.n;
  }
  
  bool operator == (long long s) const{
    if (s % N == n) return 1;
    return 0;
  }
  
  bool operator != (const Finite& finite2) const {
    return !(*this == finite2);
  }
  
  Finite operator + (const Finite& finite2) const {
    return Finite((n + finite2.n) % N);
  }
  
  Finite& operator += (const Finite& finite2) {
    n += finite2.n;
    n %= N;
    return *this;
  }
  
  Finite& operator *= (const Finite& finite2) {
    n *= finite2.n;
    n %= N;
    return *this;
  }
  
  Finite& operator++(){ //prefix
    n++;
    n %= N;
    return *this;
  }
  Finite operator++(int) {
    Finite copy = *this;
    n++;
    n %= N;
    return copy;
  }
  
  Finite& operator -= (const Finite& finite2) {
    n = n - finite2.n + N;
    n %= N;
    return *this;
  }
  
  Finite operator - (const Finite& finite2) const {
    return Finite((N + n - finite2.n) % N);
  }
  Finite operator * (const Finite& finite2) const {
    return Finite((n * finite2.n) % N);
  }
  Finite operator / (const Finite& finite2) const {
    for (int i = 0; i < N; ++i) {
      if ((finite2.n * i) % N == n) return Finite(i);
    }
    return Finite(0);
    //assert(0, "Finite / ");
  }
  Finite& operator /= (const Finite& finite2) {
    for (int i = 0; i < N; ++i) {
      if ((finite2.n * i) % N == n) {
        n = i;
        return *this;
      }
    }
    n = 0;
    return *this;
    //assert(0, "Finite / ");
  }
};

template<size_t N, size_t M, typename Field = Rational>
class Matrix{
public:
  std::vector<std::vector<Field > > matrix;

  template<typename T>
  Matrix(std::vector<std::vector<T>>& a){
    matrix.resize(N);
    for (size_t i = 0; i < N; ++i) {
      matrix[i].resize(M);
    }
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        matrix[i][j] = a[i][j];
      }
    }
  }
  
  Matrix() {
    static_assert(N == M, "create matrix n*m default constructor");
    matrix.resize(N);
    for (size_t i = 0; i < N; ++i) {
      matrix[i].resize(N, 0);
    }
    for (size_t i = 0; i < N; ++i) {
      matrix[i][i] = 1;
    }
  }
  std::vector<Field>& operator[](size_t idx){
    return matrix[idx];
  }
  const std::vector<Field>& operator[](size_t idx) const{
    return matrix[idx];
  }
  
  Matrix& operator=(const Matrix& matrix2) {
    matrix = matrix2.matrix;
    return *this;
  }
  
  std::vector<Field>& getRow(size_t idx){
    return matrix[idx];
  }
  
  std::vector<Field> getColumn(size_t idx){
    std::vector<Field> ans_col;
    for (size_t i = 0; i < N; ++i) {
      ans_col.push_back(matrix[i][idx]);
    }
    return ans_col;
  }
  
  Matrix<M, N, Field> transposed() const{
    std::vector<std::vector<Field> > for_resize;
    for_resize.resize(M);
    for (size_t i = 0; i < M; ++i) {
      for_resize[i].resize(N, 0);
    }
    Matrix<M, N, Field> ans(for_resize);
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        ans[i][j] = matrix[j][i];
      }
    }
    return ans;
  }
  
  void to_triangle_matrix(size_t& fixed_col, size_t& fixed_row, std::vector<std::vector<Field>>& copy, int idx_func, Field& ans_det, int& sign_det) const { //idx_func = 0 - det, 1 - rank, 2 - inverted
    size_t N_invert = N;
    if (idx_func == 2) N_invert *= 2;
    while(fixed_col < N && fixed_row < N) {
      size_t idx_row_max = 0;
      bool is_non_zero_in_col = 0, fl_full_0 = 1;
      for (auto u : copy[fixed_row]) {
        if (u != 0) {
          fl_full_0 = 0;
          break;
        }
      }
      if (fl_full_0 == 1) {
        fixed_row++;
        ans_det = 0;
        continue;
      }
      for (size_t i = fixed_row; i < N; ++i) {
        if (copy[i][fixed_col] != 0) {
          idx_row_max = i;
          is_non_zero_in_col = 1;
          break;
        }
      }
      if (is_non_zero_in_col == 0) {
        fixed_col++;
        if (idx_func == 0)
          ans_det = 0;
        continue;
      }
      if (fixed_row != idx_row_max) {
        swap(copy[fixed_row], copy[idx_row_max]);
        if (idx_func == 0)
          sign_det *= -1;
      }
      for (size_t i = fixed_row; i < N; ++i) {
        if (copy[i][fixed_col] == 0) continue;
        if (idx_func == 0)
          ans_det *= copy[i][fixed_col];
        for (size_t j = fixed_col + 1; j < N_invert; ++j) {
          copy[i][j] /= copy[i][fixed_col];
        }
        copy[i][fixed_col] = 1;
      }
      for (size_t i = fixed_row + 1; i < N; ++i) {
        if (copy[i][fixed_col] == 0) continue;
        for (size_t j = fixed_col; j < N_invert; ++j) {
          copy[i][j] -= copy[fixed_row][j];
        }
      }
      fixed_row++;
      fixed_col++;
    }
  }
  Field det() const{
    Field ans = 1;
    int sign = 1;
    std::vector<std::vector<Field> > copy;
    for (auto row : matrix) {
      copy.push_back(row);
    }
    size_t fixed_col = 0;
    size_t fixed_row = 0;
    to_triangle_matrix(fixed_col, fixed_row, copy, 0, ans, sign);
    ans *= sign;
    for (size_t i = 0; i < N; ++i) {
      ans *= copy[i][i];
    }
    return ans;
  }
  
  size_t rank() const {
    std::vector<std::vector<Field> > copy;
    for (auto u : matrix) {
      copy.push_back(u);
    }
    size_t ans = N;
    size_t fixed_col = 0;
    size_t fixed_row = 0;
    Field tmp1_for_triangle = 0;
    int tmp2_for_triangle = 0;
    
    to_triangle_matrix(fixed_col, fixed_row, copy, 1, tmp1_for_triangle, tmp2_for_triangle);
    for (size_t i = 0; i < N; ++i) {
      bool is_row_full_zero = 1;
      for (size_t j = 0; j < M; ++j) {
        if (copy[i][j] != 0) {
          is_row_full_zero = 0;
          break;
        }
      }
      if (is_row_full_zero == 1) {
        ans--;
      }
    }
    return ans;
  }
  
  Field trace() const {
    static_assert(N == M, "trace to non square matrix");
    Field ans = 0;
    for (size_t i = 0; i < N; ++i) {
      ans += matrix[i][i];
    }
    return ans;
  }
  
  Matrix<N, N, Field> inverted() const {
    std::vector<std::vector<Field> > copy;
    copy.resize(N);
    for (size_t i = 0; i < N; ++i) {
      copy[i].resize(2 * N, 0);
    }
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        copy[i][j] = matrix[i][j];
      }
    }
    for (size_t i = 0; i < N; ++i) {
      copy[i][i + N] = 1;
    }
    size_t fixed_col = 0;
    size_t fixed_row = 0;
    
    Field tmp1_for_triangle = 0;
    int tmp2_for_triangle = 0;

    to_triangle_matrix(fixed_col, fixed_row, copy, 2, tmp1_for_triangle, tmp2_for_triangle);
    for (int i = N - 1; i > 0; --i) {
      bool is_elem1 = 0;
      int idx = -1;
      for (size_t j = 0; j < N; ++j) {
        if (copy[i][j] == 1) is_elem1 = 1;
        if (is_elem1 == 1) {
          idx = (int)j;
          break;
        }
      }
      for (int k = 0; k < i; ++k) {
        auto coeff_str_multiply = copy[k][idx];
        for (size_t l = 0; l < 2 * N; ++l) {
          copy[k][l] -= (copy[i][l] * coeff_str_multiply);
        }
      }
    }
    Matrix<N, N, Field> ans;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        ans[i][j] = copy[i][j + N];
      }
    }
    return ans;
  }
  
  void invert() {
    Matrix<N, N, Field> copy = *this;
    *this = copy.inverted();
  }
  ~Matrix() {
    matrix.clear();
  }
};

template <size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;

template<size_t N, size_t M, typename Field = Rational>
Matrix<N, M, Field>& operator +=(Matrix<N, M, Field>& a, const Matrix<N, M, Field>& b) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      a[i][j] += b[i][j];
    }
  }
  return a;
}

template<size_t N, size_t M, typename Field = Rational>
Matrix<N, M, Field>& operator -=(Matrix<N, M, Field>& a, const Matrix<N, M, Field>& b) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      a[i][j] -= b[i][j];
    }
  }
  return a;
}

template<size_t N, size_t M = N, typename Field = Rational>
Matrix<N, M, Field> operator+(const Matrix<N, M, Field>& a, const Matrix<N, M, Field>& b) {
  Matrix<N, M, Field> copy = a;
  copy += b;
  return copy;
}

template<size_t N, size_t M = N, typename Field = Rational>
Matrix<N, M, Field> operator-(const Matrix<N, M, Field>& a, const Matrix<N, M, Field>& b) {
  Matrix<N, M, Field> copy = a;
  copy -= b;
  return copy;
}

template<size_t N, size_t M, size_t K, typename Field = Rational>
Matrix<N, K, Field> operator * (const Matrix<N, M, Field>& x, const Matrix<M, K, Field>& y) {
  std::vector<std::vector<Field> > for_resize;
  for_resize.resize(N);
  for (size_t i = 0; i < N; ++i) {
    for_resize[i].resize(K, 0);
  }
  Matrix<N, K, Field> ans(for_resize);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < K; ++j) {
      for (size_t k = 0; k < M; k++)
        ans[i][j] += x[i][k] * y[k][j];
    }
  }
  return ans;
}

template<size_t N, typename Field = Rational>
Matrix<N, N, Field> operator *= (Matrix<N, N, Field>& x, const Matrix<N, N, Field>& y) {
  x = x * y;
  return x;
}

template<size_t N, size_t M, typename Field = Rational, typename T>
Matrix<N, M, Field>& operator * (const Matrix<N, M, Field>& x, const T& y) {
  std::vector<std::vector<Field> > for_resize;
  for_resize.resize(N);
  for (size_t i = 0; i < N; ++i) {
    for_resize[i].resize(M, 0);
  }
  Matrix<N, M, Field> ans(for_resize);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      ans[i][j] = x[i][j] * y;
    }
  }
  return ans;
}

template<size_t N, size_t M, typename Field = Rational, typename T>
Matrix<N, M, Field> operator * (const T& y, const Matrix<N, M, Field>& x) {
  std::vector<std::vector<Field> > for_resize;
  for_resize.resize(N);
  for (size_t i = 0; i < N; ++i) {
    for_resize[i].resize(M, 0);
  }
  Matrix<N, M, Field> ans(for_resize);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      ans[i][j] = y * x[i][j];
    }
  }
  return ans;
}


template<size_t N1, size_t M1, typename Field1 = Rational, size_t N2, size_t M2, typename Field2 = Rational>
bool operator == (const Matrix<N1, M1, Field1>& t1, const Matrix<N2, M2, Field2>& t2) {
  if (N1 != N2 || M1 != M2) {
    return 0;
  }
  for (size_t i = 0; i < t1.matrix.size(); ++i) {
    for (size_t j = 0; j < t1[i].size(); ++j) {
      if (t1[i][j] != t2[i][j]) {
        return 0;
      }
    }
  }
  return 1;
}

template<size_t N1, size_t M1, typename Field1 = Rational, size_t N2, size_t M2, typename Field2 = Rational>
bool operator != (const Matrix<N1, M1, Field1>& t1, const Matrix<N2, M2, Field2>& t2) {
  return !(t1 == t2);
}

