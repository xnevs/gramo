#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

#include <boost/range/iterator_range.hpp>

template <
    typename Index,
    typename Compare = std::less<Index>>
class index_heap {
 private:
  std::vector<Index> heap;
  std::vector<Index> pos;
  
  Compare compare;
  
  Index count;
  
  void swap(Index i, Index j) {
    std::swap(heap[i], heap[j]);
    pos[heap[i]] = i;
    pos[heap[j]] = j;
  }
  
  void bubble_up(Index child) {
    while (child > 0) {
      Index parent = (child - 1) / 2;
      if (compare(heap[child], heap[parent])) {
        break;
      } else {
        swap(parent, child);
      }
      child = parent;
    }
  }
  
  void trickle_down(Index parent) {
    Index child;
    while ((child = 2*parent+1) < count) {
      if (child+1 < count && compare(heap[child], heap[child+1])) {
        child += 1;
      }
      if (compare(heap[child], heap[parent])) {
        break;
      } else {
        swap(parent, child);
      }
      parent = child;
    }
  }
  
 public:
  index_heap(
      Index n,
      Compare const & compare = Compare())
      : heap(n),
        pos(n, n),
        compare{compare},
        count{n} {
    std::iota(std::begin(heap), std::end(heap), 0);
    std::iota(std::begin(pos), std::end(pos), 0);
  }
  
  void heapify() {
    for (Index i=count/2-1; i >= 0; --i) {
      trickle_down(i);
    }
  }
  
  Index const & top() const {
    return heap.front();
  }
  
  bool empty() const {
    return count == 0;
  }
  
  bool full() const {
    return count == heap.size();
  }
  
  Index size() const {
    return count;
  }
  
  void push() {
    bubble_up(count++);
  }
  
  void pop() {
    --count;
    swap(0, count);
    trickle_down(0);
  }
  
  void decrease(Index item) {
    trickle_down(pos[item]);
  }
  
  void increase(Index item) {
    bubble_up(pos[item]);
  }
  
  auto pushed() {
    return boost::make_iterator_range(std::cbegin(heap), std::cbegin(heap)+count);  
  }
  
  auto popped() {
    return boost::make_iterator_range(std::crbegin(heap), std::crend(heap)-count);
  }
};
