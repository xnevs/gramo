#include <vector>

template <typename T>
class multi_stack {
 private:
  using size_type = typename std::vector<T>::size_type;

  std::vector<T> data;
  std::vector<size_type> level_counts;
  size_type count;
  size_type level;
  
 public:
  multi_stack(size_type n)
      : data(n),
        level(n),
        count{0},
        level{0} {
  }
  
  bool empty() const {
    return count == 0;
  }
  
  size_type size() const {
    return count;
  }
  
  T const & top() const {
    return data[count-1];
  }
  
  void pop() {
    --count;
  }
  
  void push(T item) {
    data[count++] = item;
  }
  
  bool level_empty() const {
    if (level == 0) {
      return count == 0;
    } else {
      return count == levels_count[level-1];
    }
  }
  
  size_type level_size() const {
    if (level == 0) {
      return count;
    } else {
      return count - level_counts[level-1];
    }
  }
  
  void push_level() {
    level_counts[level++] = count;
  }
  
  void pop_level() {
    --level;
  }
  
};
