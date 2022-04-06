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
  
  void set_sign(int x) {
    sign = x;
  }
  void remove_non_sign_zeros();
  void bitwise_subtraction(const BigInteger&, int);
  void bitwise_addition(const BigInteger&);
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

void ymn2(BigInteger& a) {
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
      ymn2(k);
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
  std::string tmp = "";
  if (sign == -1) tmp += '-';
  for (size_t i = digits.size(); i > 0; --i) {
    tmp += std::to_string(digits[i - 1]);
  }
  return tmp;
}

BigInteger::operator bool(){
  remove_non_sign_zeros();
  if ((digits.size() == 1 && digits[0] == 0) || sign == 0 || digits.size() == 0) return 0;
  return 1;
}

BigInteger::operator int(){
  remove_non_sign_zeros();
  int ans = 0;
  int degree10 = 1;
  for(size_t i = 0; i < digits.size(); ++i){
    ans += digits[i] * degree10;
    degree10 *= 10;
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


BigInteger::BigInteger(): digits(), sign(0)
{
}

std::istream& operator>>(std::istream & is, BigInteger & number_input) {
  std::ios_base::sync_with_stdio(0);
  std::cin.tie(0);
  std::string tmp;
  is >> tmp;
  if (tmp.size() == 0) {
    number_input.sign = 0;
    return is;
  }
  number_input.digits.clear();
  number_input.sign = 1;
  bool fl_0 = 1;
  for (size_t i = tmp.size(); i > 0; --i) {
    if (i == 1 && tmp[0] == '-') {
      number_input.sign = -1;
      continue;
    }
    if (tmp[i] != '0') {
      fl_0 = 0;
    }
    number_input.digits.push_back(tmp[i - 1] - '0');
  }
  if (fl_0) number_input.sign = 0;
  number_input.remove_non_sign_zeros();
  return is;
}

std::ostream& operator<<(std::ostream &os, const BigInteger &s1) {
  std::ios_base::sync_with_stdio(0);
  std::cout.tie(0);
  if (s1.digits.size() == 0) {
    os << 0;
    return os;
  }
  if (s1.sign == -1) {
    os << '-';
  }
  for (size_t i = s1.digits.size(); i > 0; --i) {
    os << s1.digits[i - 1];
  }
  return os;
}


BigInteger BigInteger::operator -(){
  BigInteger copy = *this;
  copy.sign = -1 * copy.sign;
  remove_non_sign_zeros();
  return copy;
}


void BigInteger::bitwise_subtraction(const BigInteger& s1, int flag_more) {
  int carry = 0;
  for (size_t i = 0; i < s1.digits.size() || carry != 0; ++i) {
    if (i == digits.size())
      digits.push_back (0);
    if (i < s1.digits.size()) {
      if (flag_more == 1)
        digits[i] -= carry + s1.digits[i];
      else
        digits[i] = s1.digits[i] - digits[i] - carry;
    }
    else {
      if (flag_more == 1)
        digits[i] -= carry;
      else
        digits[i] = -digits[i] - carry;
      
    }
    if (digits[i] < 0) {
      carry = 1;
      digits[i] += 10;
    }
    else {
      carry = 0;
    }
  }
}

void BigInteger::bitwise_addition(const BigInteger& s1) {
  int carry = 0;
  for (size_t i = 0; i < std::max(digits.size(), s1.digits.size()) || carry != 0; ++i) {
    if (i == digits.size())
      digits.push_back (0);
    if (i < s1.digits.size()) {
      digits[i] += carry + s1.digits[i];
    }
    else {
      digits[i] += carry;
    }
    if (digits[i] >= 10) {
      carry = 1;
      digits[i] -= 10;
    }
    else {
      carry = 0;
    }
  }
}

BigInteger& BigInteger::operator-=(const BigInteger& s1){
  if(s1.sign == 0)
    return *this;
  if(sign == 0)
    sign = 1;
  if(sign == s1.sign) {
    if (*this * sign > s1 * s1.sign) {
      bitwise_subtraction(s1, 1);
    }
    else {
      sign *= -1;
      bitwise_subtraction(s1, 0);
    }
  }
  else {
    bitwise_addition(s1);
  }
  this->remove_non_sign_zeros();
  return *this;
}


BigInteger& BigInteger::operator+=(const BigInteger& s1){
  if (sign * s1.sign != -1) {
    if (s1 > 0) sign = 1;
    if (s1 < 0) sign = -1;
    bitwise_addition(s1);
  }
  else {
    if (*this * sign > s1 * s1.sign) {
      bitwise_subtraction(s1, 1);
    }
    else {
      sign *= -1;
      bitwise_subtraction(s1, 0);
    }
  }
  this->remove_non_sign_zeros();
  return *this;
}

BigInteger& BigInteger::operator*=(int s1) {
  if (s1 < 0) {
    sign *= -1;
    s1 *= -1;
  }
  int carry = 0;
  for (size_t i = 0; i < digits.size() || carry; ++i) {
    if (i == digits.size())
      digits.push_back(0);
    long long cur = carry + 1ll * digits[i] * s1;
    digits[i] = int (cur % 10);
    carry = int (cur / 10);
  }
  remove_non_sign_zeros();
  return *this;
}

BigInteger& BigInteger::operator*=(const BigInteger& s1) {
  sign = sign * s1.sign;
  std::vector<int> y;
  for (size_t i = 0; i < digits.size(); ++i){
    int carry = 0;
    for (size_t j = 0; j < s1.digits.size() || carry; ++j) {
      int cur;
      if (y.size() <= i + j) {
        y.push_back(0);
      }
      if (j < s1.digits.size()) {
        cur = y[i + j] + digits[i] * s1.digits[j] + carry;
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


BigInteger& BigInteger::operator /= (const BigInteger& s1) {
  sign = sign * s1.sign;
  BigInteger result, current;
  BigInteger tmp = s1;
  tmp.sign = 1;
  result.digits.resize(digits.size());
  result.sign = sign;
  current.sign = 1;
  for (size_t i = digits.size(); i > 0; --i) {
    current.right();
    current.digits[0] = digits[i - 1];
    current.sign = 1;
    current.remove_non_sign_zeros();
    int new_digit = 0;
    while (new_digit < 10 && tmp * (new_digit + 1) <= current) {
      new_digit++;
    }
    result.digits[i - 1] = new_digit;
    current = current - tmp * new_digit;
  }
  *this = result;
  this->remove_non_sign_zeros();
  return *this;
}

BigInteger& BigInteger::operator%=(const BigInteger& s1) {
  *this = *this - s1 * (*this / s1);
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
  BigInteger tmp = a;
  for (int i = (int)tmp.digits.size() - 1; i >= 0; --i) {
    long long cur = tmp.digits[i] + carry * 1ll * 10;
    tmp.digits[i] = int (cur / b);
    carry = int (cur % b);
  }
  tmp.remove_non_sign_zeros();
  return tmp;
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
  
  std::string asDecimal(size_t);
  const std::string toString() const;
};

Rational::Rational(): p(0), q(1)
{
}

Rational::Rational(const BigInteger& number): p(number), q(1)
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

Rational& Rational::operator+=(const Rational& rational2) {

  p = p * rational2.q + rational2.p * q;
  q = q * rational2.q;
  bool fl_p = 0;
  bool fl_q = 0;
  if (p < 0) {
    fl_p = 1;
    p.set_sign(1);
  }
  if (q < 0) {
    fl_q = 1;
    q.set_sign(1);
  }
  BigInteger tmp = gcd(p, q);
  p /= tmp;
  q /= tmp;
  if (fl_p) {
    p.set_sign(-1);
  }
  if (fl_q) {
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

Rational& Rational::operator-=(const Rational& rational2) {

  p = p * rational2.q - rational2.p * q;
  q = q * rational2.q;
  bool fl_p = 0;
  bool fl_q = 0;
  if (p < 0) {
    fl_p = 1;
    p.set_sign(1);
  }
  if (q < 0) {
    fl_q = 1;
    q.set_sign(1);
  }
  BigInteger tmp = gcd(p, q);
  p /= tmp;
  q /= tmp;
  if (fl_p) {
    p.set_sign(-1);
  }
  if (fl_q) {
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

Rational& Rational::operator*=(const Rational& rational2){

  p *= rational2.p;
  q *= rational2.q;
  bool fl_p = 0;
  bool fl_q = 0;
  if (p < 0) {
    fl_p = 1;
    p.set_sign(1);
  }
  if (q < 0) {
    fl_q = 1;
    q.set_sign(1);
  }
  BigInteger tmp = gcd(p, q);
  p /= tmp;
  q /= tmp;
  if (fl_p) {
    p.set_sign(-1);
  }
  if (fl_q) {
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

Rational& Rational::operator/=(const Rational& rational2) {

  p *= rational2.q;
  q *= rational2.p;
  bool fl_p = 0;
  bool fl_q = 0;
  if (p < 0) {
    fl_p = 1;
    p.set_sign(1);
  }
  if (q < 0) {
    fl_q = 1;
    q.set_sign(1);
  }
  BigInteger tmp = gcd(p, q);
  p /= tmp;
  q /= tmp;
  if (fl_p) {
    p.set_sign(-1);
  }
  if (fl_q) {
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

std::istream& operator>>(std::istream & is, Rational & number_input) {
  std::ios_base::sync_with_stdio(0);
  std::cin.tie(0);
  int tmp;
  is >> tmp;
  number_input.p.digits.clear();
  number_input.p.sign = 1;
  if (tmp < 0) {
    number_input.p.sign = -1;
  }
  number_input.q = 1;
  number_input.p = abs(tmp);
  number_input.p.remove_non_sign_zeros();
  return is;
}


std::string Rational::asDecimal(size_t precision = 0) {
  std::string a = "", b = "";
  bool fl = 0;
  if (p < 0) {
    fl = 1;
    p *= -1;
  }
  for (size_t i = 0; i < precision; i++) {
    p *= 10;
  }
  a += (p / q).toString();
  if (precision == 0) {
    return a;
  }
  while(a.size() < precision) {
    a = '0' + a;
  }
  for (size_t i = 0; i < a.size() - precision; i++) {
    b += a[i];
  }
  if (b == "") {
    b += '0';
  }
  b += '.';
  for (size_t i = a.size() - precision; i < a.size(); i++) {
    b += a[i];
  }
  if (fl) {
    b = '-' + b;
    return b;
  }
  if (fl) p *= -1;
  return b;
}

Rational::operator double(){
  double a = stod(p.toString());
  double b = stod(q.toString());
  return a / b;
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

