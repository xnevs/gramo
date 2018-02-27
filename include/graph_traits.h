#ifndef GRAPH_TRAITS_H_
#define GRAPH_TRAITS_H_

struct direction_category_tag {};
struct directed_tag : public direction_category_tag {};
struct undirected_tag : public direction_category_tag {};
struct bidirectional_tag : public directed_tag {};

template <typename G>
struct is_directed {
    using D = typename G::directed_category;
    static_assert(std::is_base_of<direction_category_tag, D>::value);
    static constexpr bool value = std::is_base_of<directed_tag, D>::value;
};

#endif  // GRAPH_TRAITS_H_
