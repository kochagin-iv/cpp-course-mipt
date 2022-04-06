#include <iostream>
#include <vector>

template <bool B, typename T, typename F>
struct conditional {
  using type = F;
};

template <typename T, typename F>
struct conditional<true, T, F>{
  using type = T;
};

template <bool B, typename T, typename F>
using conditional_t = typename conditional<B, T, F>::type;

template <typename T>
class Deque {
private:
  T** arr;
  size_t sz;
  size_t cap;
  size_t first;
  size_t last;
  size_t len = 10;

  void realloc(const size_t &new_max_size);
  
  template<bool IsConst>
  struct common_iterator {
    size_t len = 10;
    T** ptr;
    T* cur_ptr;
    size_t idx;
    
    common_iterator(const T** ptr, const T* cur_ptr, size_t idx): ptr(ptr), cur_ptr(cur_ptr), idx(idx){};
    
    common_iterator(T** ptr, T* cur_ptr, size_t idx): ptr(ptr), cur_ptr(cur_ptr), idx(idx){};
    
    common_iterator(const T** ptr): ptr(ptr), cur_ptr(*ptr), idx(0){};
    
    common_iterator(const common_iterator &it) {
      ptr = it.ptr;
      cur_ptr = it.cur_ptr;
      idx = it.idx;
    }
    common_iterator& operator=(const common_iterator& it) {
      ptr = it.ptr;
      cur_ptr = it.cur_ptr;
      idx = it.idx;
      return *this;
    }
    conditional_t<IsConst, const T&, T&> operator*() {
      return *cur_ptr;
    }
    conditional_t<IsConst, const T*, T*> operator->() {
      return cur_ptr;
    }
    
    common_iterator& operator++() {
      ++idx;
      ++cur_ptr;
      if (idx != 0 && idx % len == 0) {
        ++ptr;
        idx = 0;
        cur_ptr = (*ptr);
      }
      return *this;
    }
    common_iterator& operator--() {
      if (idx == 0) {
        --ptr;
        idx = len - 1;
        cur_ptr = (*ptr) + idx;
        return *this;
      }
      --idx;
      --cur_ptr;
      return *this;
    }
    common_iterator& operator++(int) {
      T* copy_ptr(cur_ptr);
      ++this;
      return copy_ptr;
    }
    common_iterator& operator--(int) {
      T* copy_ptr(cur_ptr);
      --this;
      return copy_ptr;
    }
    common_iterator& operator+=(size_t go) {
      idx += go % len;
      if (idx >= len) {
        ++ptr;
        idx -= len;
      }
      ptr += go / len;
      cur_ptr = (*ptr) + idx;
      return *this;
    }
    common_iterator& operator-=(size_t go) {
      idx -= go % len;
      if (idx < 0) {
        --ptr;
        idx += len;
      }
      ptr -= go / len;
      cur_ptr = (*ptr) + idx;
      return *this;
    }
    common_iterator operator+(size_t go) {
      common_iterator copy(ptr, (*ptr) + idx, idx);
      copy += go;
      return copy;
    }
    common_iterator operator-(size_t go) {
      common_iterator copy(ptr, (*ptr) + idx, idx);
      copy -= go;
      return copy;
    }
    long operator-(const common_iterator& it) {
      return (ptr - it.ptr) * len + idx - it.idx;
    }
    bool operator==(const common_iterator& it) const{
      return cur_ptr == it.cur_ptr;
    }
    bool operator!=(const common_iterator& it) const{
      return cur_ptr != it.cur_ptr;
    }
    bool operator<(const common_iterator& it) const{
      return cur_ptr < it.cur_ptr;
    }
    bool operator>(const common_iterator& it) const{
      return cur_ptr > it.cur_ptr;
    }
    bool operator<=(const common_iterator& it) const{
      return !(cur_ptr > it.cur_ptr);
    }
    bool operator>=(const common_iterator& it) const{
      return !(cur_ptr < it.cur_ptr);
    }
  };
  template<typename Iterator>
  struct reverse_iter {
    Iterator iter;
    
    reverse_iter(const Iterator& iter): iter(iter) {};
    reverse_iter<Iterator>& operator ++() {
      --iter;
      return *this;
    }
    
    reverse_iter<Iterator>& operator --() {
      ++iter;
      return *this;
    }
  };
public:
  using iterator = common_iterator<false>;
  using const_iterator = common_iterator<true>;
  
  using reverse_iterator = reverse_iter<iterator>;
  using const_reverse_iterator = reverse_iter<const_iterator>;
  
  
  
  Deque();
  ~Deque();
  Deque(int t);
  Deque(int t, const T&);
  Deque(const Deque &deq);
  
  Deque& operator=(const Deque &b);
  const Deque& operator=(const Deque &b) const;
  
  
  size_t size() const;
  T& operator[](size_t pos);
  const T& operator[](size_t pos) const;
  
  T& at(size_t pos);
  const T& at(size_t pos) const;
  
  void push_back(const T& data);
  void pop_back();
  void push_front(const T& data);
  void pop_front();
  void erase(Deque::iterator it);
  void insert(Deque::iterator it, const T& elem);
  
  
  iterator begin() {
    return iterator(arr + (first / len), *(arr + (first / len)) + first % len, first % len);
  }
  iterator end() {
    return iterator(arr + (last / len), *(arr + (last / len)) + last % len, last % len);
  }
  const_iterator begin() const {
    return const_iterator(arr + (first / len), *(arr + (first / len)) + first % len, first % len);
  }
  const_iterator end() const {
    return const_iterator((arr + (last / len)), *(arr + (last / len)) + last % len, last % len);
  }
  const_iterator cbegin() const {
    return const_iterator(arr + (first / len), *(arr + (first / len)) + first % len, first % len);
  }
  const_iterator cend() const {
    return const_iterator(arr + (last / len), *(arr + (last / len)) + last % len, last % len);
  }
  
  reverse_iterator rbegin() {
    return end();
  }
  reverse_iterator rend() {
    return begin();
  }
  reverse_iterator crbegin() const{
    return cend();
  }
  reverse_iterator crend() const{
    return cbegin();
  }
  
};

template<class T>
size_t Deque<T>::size() const{
  return sz;
}


template<class T>
Deque<T>::Deque() {
  sz = 0;
  cap = 0;
  first = 0;
  last = 0;
}

template<class T>
Deque<T>::Deque(int t) {
  try {
    int kol_lvls = t / len + (t % len != 0);
    arr = reinterpret_cast<T**>(new int8_t[kol_lvls * sizeof(T*)]);
    for (int i = 0; i < kol_lvls; ++i) {
      *(arr + i) = reinterpret_cast<T*>(new int8_t[len * sizeof(T)]);
    }
    sz = t;
    cap = kol_lvls;
    first = 0;
    last = t;
  } catch (...) {
    sz = 0;
    cap = 0;
    first = 0;
    last = 0;
    delete [] arr;
    return;
  }
  
}

template<class T>
Deque<T>::Deque(int t, const T& num) {
  try {
    int kol_lvls = (t + 1) / len + ((t + 1) % len != 0);
    arr = reinterpret_cast<T**>(new int8_t[kol_lvls * sizeof(T*)]);
    for (int i = 0; i < kol_lvls; ++i) {
      *(arr + i) = reinterpret_cast<T*>(new int8_t[len * sizeof(T)]);
    }
    for (int i = 0; i < t; ++i) {
      new(*(arr + (i / len)) + i % len) T(num);
    }
    sz = t;
    cap = kol_lvls;
    first = 0;
    last = t;
  } catch (...) {
    sz = 0;
    cap = 0;
    first = 0;
    last = 0;
    delete [] arr;
  }
  
}
template<class T>
Deque<T>::Deque(const Deque<T> &deq) {
  try {
    arr = reinterpret_cast<T**>(new int8_t[deq.cap * sizeof(T*)]);
    for (size_t i = 0; i < deq.cap; ++i) {
      *(arr + i) = reinterpret_cast<T*>(new int8_t[len * sizeof(T)]);
      for (size_t j = 0; j < len; ++j) {
        new(*(arr + i) + j) T(deq.arr[i][j]);
      }
    }
    first = deq.first;
    last = deq.last;
    sz = deq.sz;
    cap = deq.cap;
  } catch (...) {
    sz = 0;
    cap = 0;
    first = 0;
    last = 0;
    delete [] arr;
    return;
  }
  
}

template<class T>
T& Deque<T>::operator[](size_t pos) {
  pos += first;
  return arr[pos / len][pos % len];
}

template<class T>
const T& Deque<T>::operator[](size_t pos) const {
  return arr[pos / len][pos % len];
}

template<class T>
T& Deque<T>::at(size_t pos) {
  if (pos >= sz) {
    throw std::out_of_range("ERROR::FAIL::AT");
  }
  return arr[pos / len][pos % len];
}

template <class T>
const T& Deque<T>::at(size_t pos) const {
  if (pos + first >= sz) {
    throw std::out_of_range("ERROR::FAIL::AT");
  }
  pos += first;
  return arr[pos / len][pos % len];
}


template<class T>
void Deque<T>::realloc(const size_t &new_max_size) {
  if (new_max_size <= cap) return;
  T** newarr = reinterpret_cast<T**>(new int8_t[new_max_size * sizeof(T*)]);
  for (size_t i = 0; i < cap; ++i) {
    *(newarr + i) = reinterpret_cast<T*>(new int8_t[len * sizeof(T)]);
  }
  for (size_t i = cap; i < 2 * cap; ++i) {
    *(newarr + i) = *(arr + i - cap);
  }
  for (size_t i = 2 * cap; i < 3 * cap; ++i) {
    *(newarr + i) = reinterpret_cast<T*>(new int8_t[len * sizeof(T)]);
  }
  first += cap * len;
  last += cap * len;
  cap *= 3;
  arr = newarr;
}

template<class T>
Deque<T>& Deque<T>::operator=(const Deque &b) {
  Deque copy(*this);
  try {
    if (sz != 0)
      delete [] reinterpret_cast<int8_t*>(arr);
    if (b.sz == 0) {
      *this = Deque();
      return *this;
    }
    cap = b.cap;
    sz = b.sz;
    arr = reinterpret_cast<T**>(new int8_t[b.cap * sizeof(T*)]);
    for (size_t i = 0; i < b.cap; ++i) {
      new(arr + i) T*(b.arr[i]);
    }
    first = b.first;
    last = b.last;
    return *this;
  } catch (...) {
    *this = copy;
    return *this;
  }
  
}


template<class T>
Deque<T>::~Deque() {
  delete [] arr;
}


template<class T>
void Deque<T>::push_back(const T& elem) {
  if (sz == 0) {
    *this = Deque<T>(1, elem);
    return;
  }
  if (last >= cap * len) {
    realloc(cap * 3);
  }
  
  *(*(arr + (last / len)) + last % len) = elem;
  ++last;
  ++sz;
  return;
}

template<class T>
void Deque<T>::pop_back() {
  --last;
  --sz;
  return;
}

template<class T>
void Deque<T>::push_front(const T &elem) {
  if (sz == 0) {
    *this = Deque<T>(1, elem);
    return;
  }
  if (first == 0) {
    realloc(cap * 3);
  }
  --first;
  ++sz;
  *(*(arr + (first / len)) + first % len) = elem;
  
}

template<class T>
void Deque<T>::pop_front() {
  ++first;
  --sz;
  return;
}

template<typename T>
void Deque<T>::erase(Deque::iterator it) {
  Deque copy(*this);
  try {
    for (auto i = it; i != end() - 1; ++i) {
      std::swap(*i, *(i + 1));
    }
    pop_back();
  } catch(...) {
    *this = copy;
  }
  
}

template<typename T>
void Deque<T>::insert(Deque::iterator it, const T& elem) {
  Deque copy(*this);
  try {
    push_back(elem);
    for (auto i = end() - 1; i != it; --i) {
      std::swap(*i, *(i - 1));
    }
  } catch(...) {
    *this = copy;
  }
}

