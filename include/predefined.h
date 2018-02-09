#ifndef PREDEFINED_H_
#define PREDEFINED_H_

#include "adjacency_set.h"
#include "adjacency_matrix.h"
#include "adjacency_list.h"
#include "adjacency_listmat.h"
#include "ordered_adjacency_listmat.h"
#include "ordered_adjacency_list_with_not_after.h"

#include "ullmann_state.h"
#include "ullmann_oalwna_state.h"
#include "simple_state.h"
#include "ri_state.h"

#include "compatibility_matrix.h"
#include "packed_compatibility_matrix.h"
#include "bitset_compatibility_matrix.h"

#include "vertex_order.h"
#include "explore.h"

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ullmann_mono(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_matrix<typename G_::index_type> g{g_};
  adjacency_matrix<typename H_::index_type> h{h_};
  
  auto index_order_g = vertex_order_DEG(g);
  
  ullmann_state_mono<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      compatibility_matrix,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ullmann_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_matrix<typename G_::index_type> g{g_};
  adjacency_matrix<typename H_::index_type> h{h_};
  
  auto index_order_g = vertex_order_DEG(g);
  
  ullmann_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      compatibility_matrix,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ullmann_mono_RDEG_CNC(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat<typename G_::index_type> gas{g_};
  auto index_order_g = vertex_order_RDEG_CNC(gas);
  
  adjacency_matrix<typename G_::index_type> g{g_};
  adjacency_matrix<typename H_::index_type> h{h_};
  
  ullmann_state_mono<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      compatibility_matrix,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ullmann_ind_RDEG_CNC(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat<typename G_::index_type> gas{g_};
  auto index_order_g = vertex_order_RDEG_CNC(gas);
  
  adjacency_matrix<typename G_::index_type> g{g_};
  adjacency_matrix<typename H_::index_type> h{h_};
  
  ullmann_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      compatibility_matrix,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ullmann_oalwna_mono(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat<typename G_::index_type> galm{g_};
  auto index_order_g = vertex_order_RDEG_CNC(galm);
  
  ordered_adjacency_list_with_not_after<typename G_::index_type> g(galm, index_order_g);
  adjacency_listmat<typename H_::index_type> h{h_};
  
  ullmann_oalwna_state_mono<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      compatibility_matrix,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void simple_mono(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_listmat<typename G_::index_type> g{g_};
  adjacency_listmat<typename H_::index_type> h{h_};

  auto index_order_g = vertex_order_DEG(g);
  
  simple_state_mono<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void simple_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_listmat<typename G_::index_type> g{g_};
  adjacency_listmat<typename H_::index_type> h{h_};

  auto index_order_g = vertex_order_RDEG_CNC(g);
  
  simple_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void simple_ind2(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_listmat<typename G_::index_type> g{g_};
  adjacency_listmat<typename H_::index_type> h{h_};

  auto index_order_g = vertex_order_RDEG_CNC(g);
  
  simple_state_ind2<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void simple_ind3(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_listmat<typename G_::index_type> g{g_};
  adjacency_listmat<typename H_::index_type> h{h_};

  auto index_order_g = vertex_order_RDEG_CNC(g);
  
  simple_state_ind3<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ri_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_list<typename G_::index_type> g{g_};
  adjacency_listmat<typename H_::index_type> h{h_};

  auto index_order_g = vertex_order_GreatestConstraintFirst(g);
  
  ri_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

#endif  // PREDEFINED_H_
