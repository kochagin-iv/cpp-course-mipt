#include <vector>
#include <iostream>

template<size_t chunkSize>
class FixedAllocator {
private:
  struct block {
    int kol_busy_in_pool = 0;
    int pool_sz = 400;
    std::vector<int8_t*> pool;
    int8_t* last_free;
    block() {
      pool.push_back(new int8_t[pool_sz * chunkSize]);
      last_free = pool[0];
    }
    void make_new_block() {
      kol_busy_in_pool = 0;
      pool.push_back(new int8_t[pool_sz * chunkSize]);
      last_free = pool.back();
    }
    ~block() {
      for (size_t i = 0; i < pool.size(); ++i) {
        delete [] pool[i];
      }
      pool.clear();
    }
  };
  std::shared_ptr<block> out_pool;
  std::vector<int8_t*> free_memory;
public:
  FixedAllocator() {
    out_pool = *new std::shared_ptr<block>(new block);
  }
  ~FixedAllocator() {
    for (size_t i = 0; i < out_pool->pool.size(); ++i) {
      delete [] out_pool->pool[i];
    }
    out_pool->pool.clear();
    free_memory.clear();
    
    
  };
  void* allocate(size_t) {
    if (free_memory.size()) {
      void* mem = free_memory.back();
      free_memory.pop_back();
      return mem;
    }
    if (out_pool->kol_busy_in_pool >= out_pool->pool_sz) {
      out_pool->make_new_block();
    }
    int8_t* last = out_pool->last_free;
    out_pool->last_free += chunkSize;
    out_pool->kol_busy_in_pool++;
    return last;
  }
  void deallocate(void* ptr, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      free_memory.push_back(static_cast<int8_t*>(ptr));
    }
  }
};

template<typename T>
class FastAllocator {
public:
  FixedAllocator<24> fixed24;
  FixedAllocator<16> fixed16;
  
  using value_type = T;
  using pointer = T*;
  using const_pointer = const T*;
  using size_type = size_t;
  
  FastAllocator() = default;
  ~FastAllocator() = default;
  template<class U>
  FastAllocator(const FastAllocator<U>& other) {
    fixed24 = other.fixed24;
    fixed16 = other.fixed16;
  }
  
  template<typename U>
  struct rebind{
    using other = FastAllocator<U>;
  };

  T* allocate(size_t n) {
    if (sizeof(T) == 24) {
      return static_cast<T*>(fixed24.allocate(n));
    }
    if (sizeof(T) == 16) {
      return static_cast<T*>(fixed16.allocate(n));
    }
    return reinterpret_cast<T*>(::operator new(n * sizeof(T)));
  }
  
  void deallocate(T* ptr, size_t n) {
    if (sizeof(T) == 16) {
      fixed16.deallocate(ptr, n);
      return;
    }
    if (sizeof(T) == 24) {
      fixed24.deallocate(ptr, n);
      return;
    }
    ::operator delete(ptr);
  }
};



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


template <typename T, typename Allocator = std::allocator<T>>
class List {
private:
  struct Node {
    Node* next{nullptr};
    Node* prev{nullptr};
    T data;
    Node() = default;
    ~Node() = default;
    Node(const Node& node) {
      this->data = node.data;
      this->next = node.next;
      this->prev = node.prev;
    }
    Node(const Node* node) {
      this->data = node->data;
      this->next = node->next;
      this->prev = node->prev;
    }
    Node(const T& value): data(value){};
    Node(Node* n, Node* p, T d): next(n), prev(p), data(d){}
  };
  Node* head{nullptr};
  Node* tail{nullptr};
  size_t sz = 0;
  template<bool IsConst>
  struct common_iterator{
    
    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type = int;
    using value_type = T;
    using pointer = conditional_t<IsConst, const Node*, Node*>;
    using reference = conditional_t<IsConst, const T&, T&>;
    
    Node* ptr;
    common_iterator(const Node* ptr): ptr(ptr){};
    common_iterator(Node* ptr): ptr(ptr){};
    common_iterator(const common_iterator<false> &it) {
      ptr = it.ptr;
    }
    common_iterator& operator=(const common_iterator& it) {
      ptr = it.ptr;
      return *this;
    }
    conditional_t<IsConst, const T&, T&> operator*() const{
      return ptr->data;
    }
    
    conditional_t<IsConst, const Node*, Node*> operator->() {
      return ptr;
    }
    common_iterator& operator++() {
      ptr = ptr->next;
      return *this;
    }
    common_iterator& operator--() {
      ptr = ptr->prev;
      return *this;
    }
    common_iterator operator++(int) {
      Node* copy_ptr(ptr);
      ++*this;
      return copy_ptr;
    }
    common_iterator operator--(int) {
      Node* copy_ptr(ptr);
      --*this;
      return copy_ptr;
    }
    bool operator==(const common_iterator& it) {
      return ptr == it.ptr;
    }
    bool operator!=(const common_iterator& it) {
      return ptr != it.ptr;
    }
    bool operator==(const common_iterator& it) const{
      return ptr == it.ptr;
    }
    bool operator!=(const common_iterator& it) const{
      return ptr != it.ptr;
    }
  };
public:
  
  using iterator = common_iterator<false>;
  using const_iterator = common_iterator<true>;
  
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  
  iterator begin() {
    return iterator(head);
  }
  iterator end() {
    return iterator(tail);
  }
  const_iterator begin() const {
    return const_iterator(head);
  }
  const_iterator end() const {
    return const_iterator(tail);
  }
  const_iterator cbegin() const {
    return const_iterator(head);
  }
  const_iterator cend() const {
    return const_iterator(tail);
  }
  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }
  reverse_iterator rend() {
    return reverse_iterator(begin());
  }
  const_reverse_iterator rbegin() const{
    return const_reverse_iterator(end());
  }
  const_reverse_iterator rend() const{
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crbegin() const{
    return const_reverse_iterator(cend());
  }
  const_reverse_iterator crend() const{
    return const_reverse_iterator(cbegin());
  }
  
  void insert(const_iterator iter, const T& data) {
    ++this->sz;
    Node* prev = iter.ptr->prev;
    Node* newNode = NodeAllocTraits::allocate(node_alloc, 1);
    newNode->data = data;
    if (prev == nullptr) {
      if (head != nullptr) {
        newNode->next = head;
        newNode->next->prev = newNode;
      } else {
        newNode->next = nullptr;
      }
      head = newNode;
      head->prev = nullptr;
    } else {
      if (prev->next == nullptr){
        prev->next = newNode;
        newNode->next = nullptr;
      } else {
        newNode->next = prev->next;
        if (newNode->next != nullptr){
          newNode->next->prev = newNode;
        }
        prev->next = newNode;
        newNode->prev = prev;
      }
    }
  }
  
  void erase(const_iterator iter) {
    --this->sz;
    Node* cur = iter.ptr;
    if (cur->prev == nullptr){
      if (cur->next == nullptr){
        head = nullptr;
      }else {
        head = cur->next;
        head->prev = nullptr;
      }
    }else {
      if (cur->next == nullptr){
        cur->prev->next = nullptr;
      }else {
        cur->prev->next = cur->next;
        cur->next->prev = cur->prev;
      }
    }
    NodeAllocTraits::deallocate(node_alloc, cur, 1);

  }
  Allocator alloc;
  typename std::allocator_traits<Allocator>::template rebind_alloc<Node> node_alloc;
  using AllocTraits = std::allocator_traits<Allocator>;
  using NodeAllocTraits = std::allocator_traits<typename Allocator::template rebind<Node>::other>;
  explicit List(const Allocator& alloc = Allocator()){
    this->alloc = alloc;
  };
  List(size_t count, const T& value, const Allocator& alloc) {
    this->sz = count;
    Node* lastNode = NodeAllocTraits::allocate(node_alloc, 1);
    lastNode->prev = nullptr;
    lastNode->data = value;
    this->head = lastNode;
    for (size_t i = 1; i < count; ++i) {
      Node* newNode = NodeAllocTraits::allocate(node_alloc, 1);
      lastNode->data = value;
      newNode->prev = lastNode;
      lastNode->next = newNode;
      lastNode = newNode;
    }
    this->tail = NodeAllocTraits::allocate(node_alloc, 1);
    this->tail->prev = lastNode;
  }
  
  List(size_t count) {
    this->sz = count;
    Node* lastNode = NodeAllocTraits::allocate(node_alloc, 1);
    lastNode->prev = nullptr;
    NodeAllocTraits::construct(node_alloc, lastNode);
    this->head = lastNode;
    for (size_t i = 1; i < count; ++i) {
      Node* newNode = NodeAllocTraits::allocate(node_alloc, 1);
      NodeAllocTraits::construct(node_alloc, newNode);
      newNode->prev = lastNode;
      lastNode->next = newNode;
      lastNode = newNode;
    }
    this->tail = NodeAllocTraits::allocate(node_alloc, 1);
    lastNode->next = this->tail;
    this->tail->prev = lastNode;
    this->tail->next = nullptr;
  }
  const Allocator& get_allocator() {
    return this->alloc;
  }
  
  
  List(const List&);
  ~List();
  List& operator=(const List&);
  size_t size() const;
  void push_back(const T& data);
  void push_front(const T& data);
  void pop_back();
  void pop_front();
  
};

template <typename T, typename Allocator>
List<T, Allocator>::~List() {
  while (this->sz != 0) {
    this->pop_back();
  }
  this->head = nullptr;
  this->tail = nullptr;
}


template <typename T, typename Allocator>
List<T, Allocator>::List(const List& list) {
  this->alloc = AllocTraits::select_on_container_copy_construction(list.alloc);
  this->head = nullptr;
  this->tail = nullptr;
  for (auto& node: list) {
    this->push_back(node);
  }
}

template <typename T, typename Allocator>
List<T, Allocator>& List<T, Allocator>::operator=(const List& list) {
  if (AllocTraits::propagate_on_container_copy_assignment::value) {
    this->alloc = list.alloc;
  }
  
  while(this->size() != 0) {
    this->pop_back();
  }
  for (auto& u: list){
    this->push_back(u);
  }
  return *this;
}

template <typename T, typename Allocator>
size_t List<T, Allocator>::size() const {
  return this->sz;
}

template <typename T, typename Allocator>
void List<T, Allocator>::push_back(const T& data) {
  ++this->sz;
  if (!this->head) {
    this->head = NodeAllocTraits::allocate(node_alloc, 1);
    this->tail = NodeAllocTraits::allocate(node_alloc, 1);
    NodeAllocTraits::construct(node_alloc, this->head, data);
    this->head->next = this->tail;
    this->head->prev = nullptr;
    this->tail->prev = this->head;
    this->tail->next = nullptr;
    return;
  }
  Node* tmp_tail = NodeAllocTraits::allocate(node_alloc, 1);
  NodeAllocTraits::construct(node_alloc, tmp_tail, data);
  tmp_tail->next = this->tail;
  this->tail->prev->next = tmp_tail;
  tmp_tail->prev = this->tail->prev;
  this->tail->prev = tmp_tail;
}

template <typename T, typename Allocator>
void List<T, Allocator>::push_front(const T& data) {
  ++this->sz;
  if (!this->head) {
    this->head = NodeAllocTraits::allocate(node_alloc, 1);
    this->tail = NodeAllocTraits::allocate(node_alloc, 1);
    NodeAllocTraits::construct(node_alloc, this->head, data);
    this->head->next = this->tail;
    this->tail->prev = this->head;
    this->tail->next = nullptr;
    this->head->prev = nullptr;
    return;
  }
  Node* tmp_head = NodeAllocTraits::allocate(node_alloc, 1);
  NodeAllocTraits::construct(node_alloc, tmp_head, data);
  tmp_head->next = this->head;
  tmp_head->prev = nullptr;
  this->head->prev = tmp_head;
  this->head = this->head->prev;
}

template <typename T, typename Allocator>
void List<T, Allocator>::pop_back() {
  if (this->sz == 1) {
    --this->sz;
    NodeAllocTraits::destroy(node_alloc, this->head);
    NodeAllocTraits::deallocate(node_alloc, this->head, 1);
    NodeAllocTraits::deallocate(node_alloc, this->tail, 1);
    this->head = nullptr;
    this->tail = nullptr;
    return;
  }
  --this->sz;
  this->tail->prev->next = nullptr;
  NodeAllocTraits::destroy(node_alloc, this->tail->prev);
  Node* tmp_tail = this->tail->prev;
  NodeAllocTraits::deallocate(node_alloc, this->tail, 1);
  this->tail = tmp_tail;
}

template <typename T, typename Allocator>
void List<T, Allocator>::pop_front() {
  if (this->sz == 1) {
    --this->sz;
    NodeAllocTraits::destroy(node_alloc, this->head);
    NodeAllocTraits::deallocate(node_alloc, this->head, 1);
    NodeAllocTraits::deallocate(node_alloc, this->tail, 1);
    this->head = nullptr;
    this->tail = nullptr;
    return;
  }
  --this->sz;
  this->head->next->prev = nullptr;
  NodeAllocTraits::destroy(node_alloc, this->head);
  Node* tmp_head = this->head->next;
  NodeAllocTraits::deallocate(node_alloc, this->head, 1);
  this->head = tmp_head;
}

