#ifndef UTILITY_H_
#define UTILITY_H_

template <
    typename G_,
    typename IndexOrder>
G_ transform_in_order(G_ const & g_, IndexOrder index_order) {
  using index_type = G_::value_type::value_type;

  auto n = g_.size();
  G_ r(n);
  
  std::vector<index_type> index_pos(n);
  for (index_type i=0; i<n; ++i) {
    index_pos[index_order[i]] = i;
  }
  
  for (index_type u=0; u<n; ++u) {
    auto uu = index_pos[u];
    r[uu].resize(g_[u].size());
    std::transform(
        std::begin(g_[u]),
        std::end(g_[u]),
        std::begin(r[uu]),
        [&index_pos](auto v) {
          return index_pos[v];
        });
    std::sort(std::begin(r[uu]), std::end(r[uu]));
  }
  
  return r;
}

#endif  // UTILITY_H_
