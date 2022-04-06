#include <cstring>
#include <iostream>

class String{
private:
  size_t size;
  size_t capacity;
  char* str;
public:
  String();
  String(char);
  String(const char*);
  String(size_t, char);
  String(const String &);
  
  char& operator[](size_t);
  const char& operator[](size_t) const;
  String& operator+=(const String&);
  String & operator=(String) ;
  String & operator=(const char *) ;
  friend bool operator<(const String &, const String &);
  friend bool operator>(const String &, const String &);
  friend bool operator==(const String &, const String &);
  
  friend std::ostream& operator<<(std::ostream &, const String&);
  friend std::istream& operator>>(std::istream&, String&);
  
  size_t find(const String&) const;
  size_t rfind(const String&) const;
  size_t length() const;
  void push_back(char);
  void pop_back();
  char& front();
  const char& front() const;
  char& back();
  const char& back() const;
  size_t find(String &);
  size_t rfind(String &);
  String substr(size_t, size_t) const;
  bool empty() const;
  void clear();
  void new_sz_str(size_t);
  size_t get_size() const{
    return size;
  }
  void swap(String&);
  const char* get_str() const;
  ~String();
};


void String::clear() {
  
  str = new char[1];
  str[0] = '\0';
  size = 0;
  capacity = 1;
}

char& String::front(){
  return str[0];
}

const char& String::front() const{
  return str[0];
}

char& String::back() {
  return str[size - 1];
}

const char& String::back() const{
  return str[size - 1];
}

void String::swap(String& s) {
  std::swap(size, s.size);
  std::swap(capacity, s.capacity);
  std::swap(str, s.str);
  
}
String::~String() {
  delete[] str;
}

const char* String::get_str() const {
  return str;
}


String String::substr(size_t start, size_t count) const {
  String tmp = "";
  for(size_t i = start; i < start + count; ++i) {
    tmp += str[i];
  }
  return tmp;
}


size_t String::length() const {
  return size;
}

String::String(){
  size = 0;
  capacity = 1;
  str = new char[1];
  str[0] = '\0';
}

String::String(size_t n, char c = '\0') : size(n), capacity(n + 1), str(new char[n + 1]) {
  memset(str, c, n);
  str[n] = '\0';
}

String::String(char st){
  str = new char[2];
  size = 1;
  capacity = 2;
  str[0] = st;
  str[1] = '\0';
}

String::String(const char* st){
  size_t len = std::strlen(st);
  str = new char[len + 1];
  size = len;
  capacity = len + 1;
  memcpy(str, st, len);
  str[len] = '\0';
}

String::String(const String &s){
  size_t len = s.length();
  size = len;
  str = new char[len + 1];
  capacity = len + 1;
  memcpy(str, s.str, size);
  str[size] = '\0';
  
}

char& String::operator[](size_t i){
  return str[i];
}

const char& String::operator[](size_t i) const {
  return str[i];
}

String& String::operator+=(const String& s1) {
  for (size_t i = 0; i < s1.size; ++i)
    push_back(s1[i]);
  return *this;
}

String operator+(const String& s1, const String& s2) {
  String copy = s1;
  copy += s2;
  return copy;
}


std::istream& operator>>(std::istream & is, String & st) {
  std::ios_base::sync_with_stdio(0);
  std::cin.tie(0);
  char tmp[10000];
  is >> tmp;
  st = tmp;
  return is;
}

std::ostream& operator<<(std::ostream &os, const String &s1) {
  std::ios_base::sync_with_stdio(0);
  std::cout.tie(0);
  for (size_t i = 0; i < s1.size; ++i) {
    os << s1[i];
  }
  return os;
}

bool String::empty() const {
  if (size == 0) return 1;
  return 0;
}

bool operator<(const String &st1, const String &st2) {
    return (std::strcmp(st1.str, st2.str) < 0);
}

bool operator>(const String &st1, const String &st2) {
    return st2 < st1;
}

bool operator==(const String &st1, const String &st2) {
  return (std::strcmp(st1.str, st2.str) == 0);
}

String & String::operator=(String st) {
  swap(st);
  return *this;
}

String & String::operator=(const char * s) {
  size = std::strlen(s);
  capacity = size + 1;
  str = new char[size + 1];
  memcpy(str, s, size);
  return *this;
}


void String::new_sz_str(size_t new_size) {
  char* newVec = new char[new_size];
  memcpy(newVec, str, size);
  std::swap(str, newVec);
  capacity = new_size;
  delete[] newVec;
}

void String::push_back(char t) {
  if (size == capacity || size == capacity - 1) {
    size_t new_size = 2 * capacity;
    if (new_size == 0) {
      new_size = 2;
    }
    new_sz_str(new_size);
  }
  str[size] = t;
  size++;
  str[size] = '\0';
}

void String::pop_back() {
  if (size == 0) return;
  str[size - 1] = '\0';
  size--;
  if (size < capacity / 2 - 1) {
    size_t new_size = capacity / 2;
    new_sz_str(new_size);
    str[size] = '\0';
  }
}

  
int equal_symbols(const String& main, const String& substring, size_t& i, int func) {
  size_t j = 0;
  size_t i1 = i - func;
  while(j < substring.get_size() && i1 < main.get_size() && main.get_str()[i1] == substring[j]){
    ++j;
    ++i1;
  }
  if (j == substring.get_size()) {
    return (int)i - func;
  }
  return -1;
}

size_t String::find(const String& substring) const{
  for(size_t i = 0; i < size; ++i) {
    int ans = equal_symbols(*this, substring, i, 0);
    if (ans != -1) return ans;
  }
  return size;
}

size_t String::rfind(const String& substring) const{
  for(size_t i = size; i > 0; --i) {
    int ans = equal_symbols(*this, substring, i, 1);
    if (ans != -1) return ans;
  }
  return size;
}

