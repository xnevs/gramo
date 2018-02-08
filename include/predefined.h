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
  adjacency_matrix g{g_};
  adjacency_matrix h{h_};
  
  auto index_order_g = vertex_order_DEG(g);
  
  ullmann_state_mono<
      adjacency_matrix,
      adjacency_matrix,
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
  adjacency_matrix g{g_};
  adjacency_matrix h{h_};
  
  auto index_order_g = vertex_order_DEG(g);
  
  ullmann_state_ind<
      adjacency_matrix,
      adjacency_matrix,
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
  
  adjacency_listmat gas{g_};
  auto index_order_g = vertex_order_RDEG_CNC(gas);
  
  adjacency_matrix g{g_};
  adjacency_matrix h{h_};
  
  ullmann_state_mono<
      adjacency_matrix,
      adjacency_matrix,
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
  
  adjacency_listmat gas{g_};
  auto index_order_g = vertex_order_RDEG_CNC(gas);
  
  adjacency_matrix g{g_};
  adjacency_matrix h{h_};
  
  ullmann_state_ind<
      adjacency_matrix,
      adjacency_matrix,
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
  
  adjacency_listmat gas{g_};
  auto index_order_g = vertex_order_RDEG_CNC(gas);
  
  ordered_adjacency_list_with_not_after g(gas, index_order_g);
  adjacency_listmat h{h_};
  
  ullmann_oalwna_state_mono<
      ordered_adjacency_list_with_not_after,
      adjacency_listmat,
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
  adjacency_listmat g{g_};
  adjacency_listmat h{h_};

  auto index_order_g = vertex_order_DEG(g);
  
  simple_state_mono<
      adjacency_listmat,
      adjacency_listmat,
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
  adjacency_listmat g{g_};
  adjacency_listmat h{h_};

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
  adjacency_listmat g{g_};
  adjacency_listmat h{h_};

  auto index_order_g = vertex_order_RDEG_CNC(g);
  
  simple_state_ind2<
      adjacency_listmat,
      adjacency_listmat,
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
  adjacency_listmat g{g_};
  adjacency_listmat h{h_};

  auto index_order_g = vertex_order_RDEG_CNC(g);
  
  simple_state_ind3<
      adjacency_listmat,
      adjacency_listmat,
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
  adjacency_list g{g_};
  adjacency_listmat h{h_};

  auto index_order_g = vertex_order_GreatestConstraintFirst(g); // TODO drugacen od RI!!
  
  ri_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

#endif  // PREDEFINED_H_
