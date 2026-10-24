/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
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

#include <algorithm>
#include <iomanip>
#include <limits>
#include <mutex>
#include <string>

#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace utils {

struct InitialPartitionerSummary {

  explicit InitialPartitionerSummary(const InitialPartitioningAlgorithm algo) :
    algorithm(algo),
    total_sum_quality(0),
    total_sum_imbalance(0),
    total_time(0.0),
    total_best(0),
    total_calls(0),
    total_extra_components(0),
    least_extra_components(std::numeric_limits<size_t>::max()) {}

  friend std::ostream & operator<< (std::ostream& str, const InitialPartitionerSummary& summary);

  void add(const InitialPartitionerSummary& summary) {
    ASSERT(algorithm == summary.algorithm);
    // Merge centered sums of squares before updating the totals (Chan's formula).
    if (total_calls > 0 && summary.total_calls > 0) {
      const double weight = static_cast<double>(total_calls) * summary.total_calls /
        (static_cast<double>(total_calls) + summary.total_calls);
      const double quality_delta = summary.average_quality() - average_quality();
      const double imbalance_delta = summary.average_imbalance() - average_imbalance();
      const double components_delta = summary.average_components() - average_components();
      quality_m2 += weight * quality_delta * quality_delta;
      imbalance_m2 += weight * imbalance_delta * imbalance_delta;
      components_m2 += weight * components_delta * components_delta;
    }
    quality_m2 += summary.quality_m2;
    imbalance_m2 += summary.imbalance_m2;
    components_m2 += summary.components_m2;
    total_sum_quality += summary.total_sum_quality;
    total_sum_imbalance += summary.total_sum_imbalance;
    total_time += summary.total_time;
    total_calls += summary.total_calls;
    total_extra_components += summary.total_extra_components;
    
    if (least_extra_components > summary.least_extra_components) {
      least_extra_components = summary.least_extra_components;
    }
  }

  double average_quality() const {
    return static_cast<double>(total_sum_quality) / std::max(total_calls, UL(1));
  }

  double average_imbalance() const {
    return static_cast<double>(total_sum_imbalance) / std::max(total_calls, UL(1));
  }

  double average_components() const {
    return static_cast<double>(total_extra_components) / std::max(total_calls, UL(1));
  }

  // Population variance over all recorded runs; zero for zero or one run.
  double variance_quality() const {
    return quality_m2 / std::max(total_calls, UL(1));
  }

  double variance_imbalance() const {
    return imbalance_m2 / std::max(total_calls, UL(1));
  }

  double variance_components() const {
    return components_m2 / std::max(total_calls, UL(1));
  }

  double average_running_time() const {
    return static_cast<double>(total_time) / std::max(total_calls, UL(1));
  }

  double percentage_best(const size_t total_ip_calls) const {
    return ( static_cast<double>(total_best) / std::max(total_ip_calls, UL(1)) ) * 100.0;
  }

  InitialPartitioningAlgorithm algorithm;
  double total_sum_quality;
  double total_sum_imbalance;
  double total_time;
  size_t total_best;
  size_t total_calls;
  size_t total_extra_components;
  size_t least_extra_components;
  double quality_m2 = 0.0;
  double imbalance_m2 = 0.0;
  double components_m2 = 0.0;
};

inline std::ostream & operator<< (std::ostream& str, const InitialPartitionerSummary& summary) {
  str << " avg_quality_" << summary.algorithm << "=" << summary.average_quality()
      << " var_quality_" << summary.algorithm << "=" << summary.variance_quality()
      << " avg_imbalance_" << summary.algorithm << "=" << summary.average_imbalance()
      << " var_imbalance_" << summary.algorithm << "=" << summary.variance_imbalance()
      << " avg_components_" << summary.algorithm << "=" << summary.average_components()
      << " var_components_" << summary.algorithm << "=" << summary.variance_components()
      << " total_time_" << summary.algorithm << "=" << summary.total_time
      << " total_best_" << summary.algorithm << "=" << summary.total_best;
  return str;
}

class InitialPartitioningStats {

 public:
  explicit InitialPartitioningStats() :
    _stat_mutex(),
    _num_initial_partitioner(static_cast<uint8_t>(InitialPartitioningAlgorithm::UNDEFINED)),
    _ip_summary(),
    _total_ip_calls(0),
    _total_sum_number_of_threads(0) {
    for ( uint8_t algo = 0; algo < _num_initial_partitioner; ++algo ) {
      _ip_summary.emplace_back(static_cast<InitialPartitioningAlgorithm>(algo));
    }
  }

  InitialPartitioningStats(const InitialPartitioningStats& other) :
    _stat_mutex(),
    _num_initial_partitioner(other._num_initial_partitioner),
    _ip_summary(other._ip_summary),
    _total_ip_calls(other._total_ip_calls),
    _total_sum_number_of_threads(other._total_sum_number_of_threads) { }

  InitialPartitioningStats & operator= (const InitialPartitioningStats &) = delete;

  InitialPartitioningStats(InitialPartitioningStats&& other) :
    _stat_mutex(),
    _num_initial_partitioner(std::move(other._num_initial_partitioner)),
    _ip_summary(std::move(other._ip_summary)),
    _total_ip_calls(std::move(other._total_ip_calls)),
    _total_sum_number_of_threads(std::move(other._total_sum_number_of_threads)) { }

  InitialPartitioningStats & operator= (InitialPartitioningStats &&) = delete;

  void add_initial_partitioning_result(const InitialPartitioningAlgorithm best_algorithm,
                                       const size_t number_of_threads,
                                       const parallel::scalable_vector<InitialPartitionerSummary>& summary) {
    std::lock_guard<std::mutex> lock(_stat_mutex);
    ASSERT(summary.size() == _ip_summary.size());
    uint8_t best_algorithm_index = static_cast<uint8_t>(best_algorithm);
    ++_ip_summary[best_algorithm_index].total_best;

    for ( size_t i = 0; i < _num_initial_partitioner; ++i ) {
      _ip_summary[i].add(summary[i]);
    }

    _total_sum_number_of_threads += number_of_threads;
    ++_total_ip_calls;
  }

  double average_number_of_threads_per_ip_call() const {
    return static_cast<double>(_total_sum_number_of_threads) / std::max(_total_ip_calls, UL(1));
  }

  void printInitialPartitioningStats() {
    LOG << "Initial Partitioning Algorithm Summary:";
    LOG << "Number of Initial Partitioning Calls =" << _total_ip_calls;
    LOG << "Average Number of Thread per IP Call ="
        << average_number_of_threads_per_ip_call() << "\n";
    for ( const InitialPartitionerSummary& summary : _ip_summary ) {
      if (summary.total_calls == 0) continue;
      LOG << summary.algorithm << " (runs=" << summary.total_calls
          << ", time=" << summary.total_time << "s, best=" << summary.total_best
          << ", " << summary.percentage_best(_total_ip_calls) << "%)";
      LOG << "  " << std::left << std::setw(20) << "Metric"
          << std::setw(18) << "Mean" << "Variance";
      LOG << "  " << std::left << std::setw(20) << "Quality"
          << std::setw(18) << summary.average_quality() << summary.variance_quality();
      LOG << "  " << std::left << std::setw(20) << "Imbalance"
          << std::setw(18) << summary.average_imbalance() << summary.variance_imbalance();
      LOG << "  " << std::left << std::setw(20) << "Extra components"
          << std::setw(18) << summary.average_components() << summary.variance_components();
    }
  }

  friend std::ostream & operator<< (std::ostream& str, const InitialPartitioningStats& stats);

 private:
  std::mutex _stat_mutex;
  const uint8_t _num_initial_partitioner;
  parallel::scalable_vector<InitialPartitionerSummary> _ip_summary;
  size_t _total_ip_calls;
  size_t _total_sum_number_of_threads;
};

inline std::ostream & operator<< (std::ostream& str, const InitialPartitioningStats& stats) {
  str << " average_number_of_threads_per_ip_call="
      << stats.average_number_of_threads_per_ip_call();
  for ( const InitialPartitionerSummary& summary : stats._ip_summary ) {
    str << summary;
  }
  return str;
}

}  // namespace utils
}  // namespace mt_kahypar
