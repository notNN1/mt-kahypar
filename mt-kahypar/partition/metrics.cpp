/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "mt-kahypar/partition/metrics.h"

#include <cmath>
#include <algorithm>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/utils/exception.h"
#include "mt-kahypar/partition/connected_components/compute_components.h"

namespace mt_kahypar {

int BalanceMetrics::numViolations() const {
  return (this->violates_balance ? 1 : 0) + (this->violates_non_empty_blocks ? 1 : 0);
}

bool BalanceMetrics::isValidPartition() const {
  return this->numViolations() == 0;
}

bool BalanceMetrics::isBetter(const BalanceMetrics& other) const {
  return this->numViolations() < other.numViolations() ||
    (this->numViolations() == other.numViolations() && imbalance_value < other.imbalance_value);
}

bool BalanceMetrics::isEqual(const BalanceMetrics& other) const {
  return this->numViolations() == other.numViolations() && imbalance_value == other.imbalance_value;
}

bool BalanceMetrics::operator==(const BalanceMetrics& other) const {
  return imbalance_value == other.imbalance_value &&
    violates_balance == other.violates_balance &&
    violates_non_empty_blocks == other.violates_non_empty_blocks;
}

bool ConnectivityMetrics::isBetter(const ConnectivityMetrics& other) const {
  return this->extra_components_count < other.extra_components_count;
}

bool ConnectivityMetrics::isEqual(const ConnectivityMetrics& other) const {
  return this->extra_components_count == other.extra_components_count;
}

bool ConnectivityMetrics::operator==(const ConnectivityMetrics& other) const {
  return this->extra_components_count == other.extra_components_count;
}

int ConnectivityMetrics::numViolations() const {
  return this->extra_components_count > 0 ? 1 : 0 + this->inefficient_components_count > 0 ? 1 : 0;
}

bool Metrics::isBetter(const Metrics& other) const {
  return this->to_tuple() < other.to_tuple();
}

std::tuple<size_t, bool, size_t, size_t, double, HyperedgeWeight> Metrics::to_tuple() const {
  return std::tuple {
    this->imbalance.numViolations(),
    this->imbalance.violates_balance,
    this->connectivity.inefficient_components_count,
    this->connectivity.extra_components_count,
    this->imbalance.imbalance_value,
    this->quality
  };
} 

bool Metrics::isEqual(const Metrics& other) const {
  return quality == other.quality && imbalance.isEqual(other.imbalance) && connectivity.isEqual(other.connectivity);
}

int Metrics::numViolations() const {
  return this->imbalance.numViolations() + this->connectivity.numViolations();
}

void Metrics::log() const {
  LOG << "Violates balance:          " << this->imbalance.violates_balance;
  LOG << "Violates non empty blocks: " << this->imbalance.violates_non_empty_blocks;
  LOG << "Imbalance:                 " << this->imbalance.imbalance_value;
  LOG << "Quality:                   " << this->quality;
  LOG << "Extra components:          " << this->connectivity.extra_components_count;
  LOG << "Inefficient components:    " << this->connectivity.inefficient_components_count;
}

std::ostream& operator<< (std::ostream& os, const BalanceMetrics& imbalance) {
  return os << "[imb: " << imbalance.imbalance_value  << "; " << (imbalance.isValidPartition() ? "valid" : "invalid") << "]";
}


namespace metrics {

namespace {

template<typename PartitionedHypergraph, Objective objective>
struct ObjectiveFunction { };

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::cut> {
  HyperedgeWeight operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    return phg.connectivity(he) > 1 ? phg.edgeWeight(he) : 0;
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::km1> {
  HyperedgeWeight operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    return std::max(phg.connectivity(he) - 1, 0) * phg.edgeWeight(he);
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::soed> {
  HyperedgeWeight operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    const PartitionID connectivity = phg.connectivity(he);
    return connectivity > 1 ? connectivity * phg.edgeWeight(he) : 0;
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::steiner_tree> {
  HyperedgeWeight operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    ASSERT(phg.hasTargetGraph());
    const TargetGraph* target_graph = phg.targetGraph();
    const HyperedgeWeight distance = target_graph->distance(phg.shallowCopyOfConnectivitySet(he));
    return distance * phg.edgeWeight(he);
  }
};

template<Objective objective, typename PartitionedHypergraph>
HyperedgeWeight compute_objective_parallel(const PartitionedHypergraph& phg) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  tbb::enumerable_thread_specific<HyperedgeWeight> obj(0);
  phg.doParallelForAllEdges([&](const HyperedgeID he) {
    obj.local() += func(phg, he);
  });
  return obj.combine(std::plus<>()) / (PartitionedHypergraph::is_graph ? 2 : 1);
}

template<Objective objective, typename PartitionedHypergraph>
HyperedgeWeight compute_objective_sequentially(const PartitionedHypergraph& phg) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  HyperedgeWeight obj = 0;
  for (const HyperedgeID& he : phg.edges()) {
    obj += func(phg, he);
  }
  return obj / (PartitionedHypergraph::is_graph ? 2 : 1);
}

template<Objective objective, typename PartitionedHypergraph>
HyperedgeWeight contribution(const PartitionedHypergraph& phg, const HyperedgeID he) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  return func(phg, he);
}

}

template<typename PartitionedHypergraph>
HyperedgeWeight quality(const PartitionedHypergraph& hg,
                        const Context& context,
                        const bool parallel) {
  return quality(hg, context.partition.objective, parallel);
}

template<typename PartitionedHypergraph>
HyperedgeWeight quality(const PartitionedHypergraph& hg,
                        const Objective objective,
                        const bool parallel) {
  switch (objective) {
    case Objective::cut:
      return parallel ? compute_objective_parallel<Objective::cut>(hg) :
        compute_objective_sequentially<Objective::cut>(hg);
    case Objective::km1:
      return parallel ? compute_objective_parallel<Objective::km1>(hg) :
        compute_objective_sequentially<Objective::km1>(hg);
    case Objective::soed:
      return parallel ? compute_objective_parallel<Objective::soed>(hg) :
        compute_objective_sequentially<Objective::soed>(hg);
    case Objective::steiner_tree:
      return parallel ? compute_objective_parallel<Objective::steiner_tree>(hg) :
        compute_objective_sequentially<Objective::steiner_tree>(hg);
    default: throw InvalidParameterException("Unknown Objective");
  }
  return 0;
}

template<typename PartitionedHypergraph>
HyperedgeWeight contribution(const PartitionedHypergraph& hg,
                             const HyperedgeID he,
                             const Objective objective) {
  switch (objective) {
    case Objective::cut: return contribution<Objective::soed>(hg, he);
    case Objective::km1: return contribution<Objective::km1>(hg, he);
    case Objective::soed: return contribution<Objective::soed>(hg, he);
    case Objective::steiner_tree: return contribution<Objective::steiner_tree>(hg, he);
    default: throw InvalidParameterException("Unknown Objective");
  }
  return 0;
}

template<typename PartitionedHypergraph>
bool isValidPartition(const PartitionedHypergraph& phg, const Context& context) {
  BalanceMetrics imbalance_metrics = imbalance(phg, context);
  return imbalance_metrics.isValidPartition();
}

template<typename PartitionedHypergraph>
BalanceMetrics imbalance(const PartitionedHypergraph& hypergraph, const Context& context) {
  ASSERT(context.partition.perfect_balance_part_weights.size() == (size_t)context.partition.k);

  size_t num_empty_parts = 0;
  double max_balance = 0.0;
  bool violates_balance = false;
  for (PartitionID i = 0; i < context.partition.k; ++i) {
    const HypernodeWeight part_weight = hypergraph.partWeight(i);
    const double balance_i = (part_weight
            / static_cast<double>(context.partition.perfect_balance_part_weights[i]));
    max_balance = std::max(max_balance, balance_i);
    if (part_weight > context.partition.max_part_weights[i]) {
      violates_balance = true;
    }
    if (part_weight == 0) {
      num_empty_parts++;
    }
  }
  bool too_many_empty_parts = num_empty_parts > hypergraph.numRemovedHypernodes();
  return BalanceMetrics{max_balance - 1.0, violates_balance,
    !context.partition.allow_empty_blocks && too_many_empty_parts};
}

template<typename PartitionedHypergraph>
double approximationFactorForProcessMapping(const PartitionedHypergraph& hypergraph, const Context& context) {
  if ( !PartitionedHypergraph::is_graph ) {
    tbb::enumerable_thread_specific<HyperedgeWeight> approx_factor(0);
    hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
      const size_t connectivity = hypergraph.connectivity(he);
      approx_factor.local() += connectivity <= context.mapping.max_steiner_tree_size ? 1 : 2;
    });
    return static_cast<double>(approx_factor.combine(std::plus<>())) / hypergraph.initialNumEdges();
  } else {
    return 1.0;
  }
}

template<typename PartitionedHypergraph>
ConnectivityMetrics connectivity(const PartitionedHypergraph& phg, const Context& context) {

  vec<vec<mt_kahypar::connected_components::ConnectedComponent>> components;
  mt_kahypar::connected_components::compute_components_per_block(phg, context, components);

  // calculate the total number of components
  u_int32_t total_component_count = 0;
  for (const vec<mt_kahypar::connected_components::ConnectedComponent>& components_per_partition : components) {
    total_component_count += components_per_partition.size();
  }

  // calculate inefficient components
  vec<vec<connected_components::ComponentInfo>> super_components;
    connected_components::compute_super_components(phg, context, components, super_components);

    vec<vec<connected_components::ComponentInfo>> filtered_super_components;
    connected_components::find_inefficient_super_components(context, super_components, filtered_super_components);

  u_int32_t extra_components_count        = total_component_count - phg.k();
  u_int32_t inefficient_components_count  = filtered_super_components.size(); 
  return { extra_components_count, inefficient_components_count };
}

namespace {
#define OBJECTIVE_1(X) HyperedgeWeight quality(const X& hg, const Context& context, const bool parallel)
#define OBJECTIVE_2(X) HyperedgeWeight quality(const X& hg, const Objective objective, const bool parallel)
#define CONTRIBUTION(X) HyperedgeWeight contribution(const X& hg, const HyperedgeID he, const Objective objective)
#define IS_VALID_PARTITION(X) bool isValidPartition(const X& phg, const Context& context)
#define IMBALANCE(X) BalanceMetrics imbalance(const X& hypergraph, const Context& context)
#define APPROX_FACTOR(X) double approximationFactorForProcessMapping(const X& hypergraph, const Context& context)
#define CONNECTIVITY(X) ConnectivityMetrics connectivity(const X& hypergraph, const Context& context)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(OBJECTIVE_1)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(OBJECTIVE_2)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(CONTRIBUTION)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IS_VALID_PARTITION)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IMBALANCE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(APPROX_FACTOR)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(CONNECTIVITY)

} // namespace metrics

} // namespace mt_kahypar 