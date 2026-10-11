/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2026 Simon Lang <simonlang2>@kit.edu
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

#include "mt-kahypar/partition/initial_partitioning/tarjan_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/initial_partitioning_data_container.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/connected_components/tarjan_region_coarsening.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

template<typename TypeTraits>
void TarjanInitialPartitioner<TypeTraits>::partitionImpl() {
  if ( _ip_data.should_initial_partitioner_run(InitialPartitioningAlgorithm::tarjan) ) {
    const auto& enabled_ip_algos = _context.initial_partitioning.enabled_ip_algos_tarjan.empty() ?
      _context.initial_partitioning.enabled_ip_algos : _context.initial_partitioning.enabled_ip_algos_tarjan;

    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();
    std::uniform_int_distribution<PartitionID> select_random_block(0, _context.partition.k - 1);
    
    connected_components::TarjanRegionCoarsening<PartitionedHypergraph> region_coarsening(_context, _rng);
    auto communities = region_coarsening.computeCommunities(hg);

    auto& old_hg    = hg.hypergraph();
    auto new_hg     = old_hg.contract(communities, true);

    PartitionedHypergraph new_phg(hg.k(), new_hg);
    ////


    //// do initial partitioning on the new phg for each enabled initial partitioner
    InitialPartitioningDataContainer<TypeTraits> ip_data(new_phg, _context);

    auto* ip_data_ptr = ip::to_pointer(ip_data);

    for ( uint8_t i = 0; i < static_cast<uint8_t>(InitialPartitioningAlgorithm::UNDEFINED); ++i ) {
      if ( static_cast<InitialPartitioningAlgorithm>(i) != InitialPartitioningAlgorithm::tarjan 
      && enabled_ip_algos[i]) {
        auto algorithm = static_cast<InitialPartitioningAlgorithm>(i);

        // Create one initial partitioner.
        std::unique_ptr<IInitialPartitioner> initial_partitioner =
            InitialPartitionerFactory::getInstance().createObject(
                algorithm,
                algorithm,
                ip_data_ptr,
                _context,
                _context.partition.seed + i + 100,
                10000);

        // Run it.
        initial_partitioner->partition();
      }
    }
    ////


    //// Do partitioning like the best performing partitioner on the uncontracted graph
    auto& ip_data_ref = *reinterpret_cast<InitialPartitioningDataContainer<TypeTraits>*>(ip_data_ptr);
    ip_data_ref.apply();
    auto& local_phg = ip_data_ref.local_partitioned_hypergraph();

    for (const HypernodeID& node : hg.nodes()) {
      hg.setNodePart(node, local_phg.partID(communities[node]));
    }
    ////

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(InitialPartitioningAlgorithm::tarjan, _rng, _tag, time);
  }
};


INSTANTIATE_CLASS_WITH_TYPE_TRAITS(TarjanInitialPartitioner)

} // namespace mt_kahypar
