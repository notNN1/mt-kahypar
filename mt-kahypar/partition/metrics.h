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

#pragma once

#include <iostream>
#include <limits>

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/connected_components/compute_components.h"

namespace mt_kahypar {

struct BalanceMetrics {
  double imbalance_value = std::numeric_limits<double>::max();
  bool violates_balance = true;
  bool violates_non_empty_blocks = true;

  bool isValidPartition() const;

  bool isBetter(const BalanceMetrics& other) const;

  bool isEqual(const BalanceMetrics& other) const;

  bool operator==(const BalanceMetrics& other) const;

  int numViolations() const;
};


struct ConnectivityMetrics {
  u_int32_t extra_components_count;                      // #components - k
  u_int32_t inefficient_components_count;

  bool isBetter(const ConnectivityMetrics& other) const;

  bool isEqual(const ConnectivityMetrics& other) const;

  bool operator==(const ConnectivityMetrics& other) const;

  int numViolations() const;
};

std::ostream& operator<< (std::ostream& os, const BalanceMetrics& imbalance);

struct Metrics {
  HyperedgeWeight     quality;
  BalanceMetrics      imbalance;
  ConnectivityMetrics connectivity;

  bool isBetter(const Metrics& other) const;

  bool isEqual(const Metrics& other) const;

  void log() const;

  int numViolations() const;

  std::tuple<bool, bool, size_t, size_t, double, HyperedgeWeight> to_tuple() const; 
};

namespace metrics {

// ! Computes for the given partitioned hypergraph the corresponding objective function
template<typename PartitionedHypergraph>
HyperedgeWeight quality(const PartitionedHypergraph& hg,
                        const Context& context,
                        const bool parallel = true);
template<typename PartitionedHypergraph>
HyperedgeWeight quality(const PartitionedHypergraph& hg,
                        const Objective objective,
                        const bool parallel = true);

// ! Computes for a hyperedge the contribution to the corresponding objective function
template<typename PartitionedHypergraph>
HyperedgeWeight contribution(const PartitionedHypergraph& hg,
                             const HyperedgeID he,
                             const Objective objective);

template<typename PartitionedHypergraph>
bool isValidPartition(const PartitionedHypergraph& phg, const Context& context);

template<typename PartitionedHypergraph>
BalanceMetrics imbalance(const PartitionedHypergraph& hypergraph, const Context& context);

template<typename PartitionedHypergraph>
double approximationFactorForProcessMapping(const PartitionedHypergraph& hypergraph, const Context& context);

template<typename PartitionedHypergraph>
ConnectivityMetrics connectivity(const PartitionedHypergraph& phg,
                    const Context& context);

}  // namespace metrics
}  // namespace mt_kahypar
