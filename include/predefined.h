#ifndef PREDEFINED_H_
#define PREDEFINED_H_

#include "adjacency_set.h"
#include "adjacency_matrix.h"
#include "adjacency_list.h"
#include "adjacency_list_in_order.h"
#include "ordered_adjacency_list.h"
#include "ordered_adjacency_list_with_not_after.h"
#include "adjacency_listmat.h"
#include "adjacency_listmat_with_not.h"
#include "ordered_adjacency_listmat.h"
#include "ordered_adjacency_listmat_with_not_after.h"
#include "orderable_adjacency_listmat.h"

#include "ullmann_state.h"
#include "ullmann_oalwna_state.h"
#include "neighborhood_filter_state.h"
#include "ullimp_state.h"
#include "ullimp_no_after_state.h"
#include "ullimp2_state.h"
#include "ullimp3_state.h"
#include "ullimp4_state.h"
#include "simple_state.h"
#include "ri_state.h"
#include "ri2_state.h"
#include "ullimp_ri_state.h"
#include "ri_lookahead_state.h"
#include "refined_ri_state.h"
#include "ri_dynamic_parent_state.h"
#include "dynamic_state.h"
#include "dynamic_sorted_vector_state.h"
#include "dynamic_mat_state.h"

#include "compatibility_matrix.h"
#include "packed_compatibility_matrix.h"
#include "bitset_compatibility_matrix.h"
#include "reduced_compatibility_matrix.h"
#include "reduced_compatibility_matrix2.h"
#include "reduced_compatibility_matrix2_with_count.h"

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
      compatibility_matrix<typename decltype(g)::index_type, typename decltype(h)::index_type>,
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
      compatibility_matrix<typename decltype(g)::index_type, typename decltype(h)::index_type>,
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
      compatibility_matrix<typename decltype(g)::index_type, typename decltype(h)::index_type>,
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
      compatibility_matrix<typename decltype(g)::index_type, typename decltype(h)::index_type>,
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
      compatibility_matrix<typename decltype(g)::index_type, typename decltype(h)::index_type>,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void neighborhood_filter_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat<typename G_::index_type> galm{g_};
  auto index_order_g = vertex_order_RDEG_CNC(galm);
  //auto index_order_g = vertex_order_GreatestConstraintFirst(galm);
  
  ordered_adjacency_listmat_with_not_after<typename G_::index_type> g(g_, index_order_g);
  adjacency_listmat_with_not<typename H_::index_type> h{h_};
  
  neighborhood_filter_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      reduced_compatibility_matrix2<typename decltype(g)::index_type, typename decltype(h)::index_type>,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ullimp_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat<typename G_::index_type> galm{g_};
  //auto index_order_g = vertex_order_RDEG_CNC(galm);
  auto index_order_g = vertex_order_GreatestConstraintFirst(galm);
  
  ordered_adjacency_listmat_with_not_after<typename G_::index_type> g{g_, index_order_g};
  adjacency_listmat_with_not<typename H_::index_type> h{h_};
  
  ullimp_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      reduced_compatibility_matrix2<typename decltype(g)::index_type, typename decltype(h)::index_type>,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ullimp_no_after_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat_with_not<typename G_::index_type> g{g_};
  adjacency_listmat_with_not<typename H_::index_type> h{h_};
  
  auto index_order_g = vertex_order_GreatestConstraintFirst(g);
  
  ullimp_no_after_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      reduced_compatibility_matrix2<typename decltype(g)::index_type, typename decltype(h)::index_type>,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ullimp2_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat<typename G_::index_type> galm{g_};
  auto index_order_g = vertex_order_RDEG_CNC(galm);
  
  adjacency_matrix<typename G_::index_type> g{g_};
  adjacency_matrix<typename H_::index_type> h{h_};
  
  ullimp2_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      compatibility_matrix<typename decltype(g)::index_type, typename decltype(h)::index_type>,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ullimp3_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat<typename G_::index_type> galm{g_};
  //auto index_order_g = vertex_order_RDEG_CNC(galm);
  auto index_order_g = vertex_order_GreatestConstraintFirst(galm);
  
  adjacency_matrix<typename G_::index_type> g{g_};
  adjacency_matrix<typename H_::index_type> h{h_};
  
  ullimp3_state_ind<
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
void ullimp4_mono(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat<typename G_::index_type> galm{g_};
  auto index_order_g = vertex_order_RDEG_CNC(galm);
  
  ordered_adjacency_list<typename G_::index_type> g(g_, index_order_g);
  adjacency_listmat<typename H_::index_type> h{h_};
  
  ullimp4_state_mono<
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
void ullimp4_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat<typename G_::index_type> galm{g_};
  //auto index_order_g = vertex_order_GreatestConstraintFirst(galm);
  auto index_order_g = vertex_order_GreatestConstraintFirst(galm);
  
  ordered_adjacency_list<typename G_::index_type> g(g_, index_order_g);
  adjacency_listmat<typename H_::index_type> h{h_};
  
  ullimp4_state_ind<
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
void ullimp4_ind2(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat<typename G_::index_type> galm{g_};
  auto index_order_g = vertex_order_GreatestConstraintFirst(galm);
  
  ordered_adjacency_list<typename G_::index_type> g(g_, index_order_g);
  adjacency_listmat<typename H_::index_type> h{h_};
  
  ullimp4_state_ind2<
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

  auto index_order_g = vertex_order_GreatestConstraintFirst(g);
  
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

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ri_RDEG_CNC_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_listmat<typename G_::index_type> g{g_};
  adjacency_listmat<typename H_::index_type> h{h_};

  auto index_order_g = vertex_order_RDEG_CNC(g);
  
  ri_state_ind<
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
void ri_RDEG_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_list<typename G_::index_type> g{g_};
  adjacency_listmat<typename H_::index_type> h{h_};

  auto index_order_g = vertex_order_RDEG(g);
  
  ri_state_ind<
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
void ri_lookahead_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_list<typename G_::index_type> g{g_};
  adjacency_listmat<typename H_::index_type> h{h_};

  auto index_order_g = vertex_order_GreatestConstraintFirst(g);
  
  ri_lookahead_state_ind<
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
void refined_ri_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_list<typename G_::index_type> g{g_};
  adjacency_listmat<typename H_::index_type> h{h_};

  auto index_order_g = vertex_order_GreatestConstraintFirst(g);
  
  refined_ri_state_ind<
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
void ri_dynamic_parent_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
    
  adjacency_list<typename G_::index_type> gal{g_};
  auto index_order_g = vertex_order_GreatestConstraintFirst(gal);
  
  ordered_adjacency_listmat<typename G_::index_type> g{g_, index_order_g};
  adjacency_listmat<typename H_::index_type> h{h_};
  
  ri_dynamic_parent_state_ind<
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
void ri2_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_list<typename G_::index_type> gal{g_};

  auto index_order_g = vertex_order_GreatestConstraintFirst(gal);
  
  ordered_adjacency_list<typename G_::index_type> g(g_, index_order_g);
  adjacency_listmat<typename H_::index_type> h{h_};
  
  ri2_state_ind<
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
void ri2_ind2(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_list<typename G_::index_type> gal{g_};

  auto index_order_g = vertex_order_GreatestConstraintFirst(gal);
  
  ordered_adjacency_list<typename G_::index_type> g(g_, index_order_g);
  adjacency_listmat<typename H_::index_type> h{h_};
  
  ri2_state_ind2<
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
void ullimp_ri_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat<typename G_::index_type> galm{g_};
  auto index_order_g = vertex_order_RDEG_CNC(galm);
  
  ordered_adjacency_listmat_with_not_after<typename G_::index_type> g{g_, index_order_g};
  adjacency_listmat_with_not<typename H_::index_type> h{h_};
  
  ullimp_ri_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      reduced_compatibility_matrix2<typename decltype(g)::index_type, typename decltype(h)::index_type>,
      decltype(index_order_g)> S{g, h, vertex_comp, edge_comp, index_order_g};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void dynamic_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat_with_not<typename G_::index_type> g{g_};
  adjacency_listmat_with_not<typename H_::index_type> h{h_};
  
  dynamic_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate> S{g, h, vertex_comp, edge_comp};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void dynamic_sorted_vector_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  orderable_adjacency_listmat/*_with_not*/<typename G_::index_type> g{g_};
  adjacency_listmat_with_not<typename H_::index_type> h{h_};
  
  dynamic_sorted_vector_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate> S{g, h, vertex_comp, edge_comp};
  
  explore(S, callback);
}

template <
    typename G_,
    typename H_,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void dynamic_mat_ind(
    G_ const & g_,
    H_ const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  
  adjacency_listmat_with_not<typename G_::index_type> g{g_};
  adjacency_listmat_with_not<typename H_::index_type> h{h_};
  
  dynamic_mat_state_ind<
      decltype(g),
      decltype(h),
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      reduced_compatibility_matrix2_with_count<typename decltype(g)::index_type, typename decltype(h)::index_type>> S{g, h, vertex_comp, edge_comp};
  
  explore(S, callback);
}

#endif  // PREDEFINED_H_
