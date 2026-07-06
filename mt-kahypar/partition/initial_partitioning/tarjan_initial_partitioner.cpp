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

#include "mt-kahypar/partition/initial_partitioning/random_initial_partitioner.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

template<typename TypeTraits>
void TarjanInitialPartitioner<TypeTraits>::partitionImpl() {
  if ( _ip_data.should_initial_partitioner_run(InitialPartitioningAlgorithm::tarjan) ) {
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();
    std::uniform_int_distribution<PartitionID> select_random_block(0, _context.partition.k - 1);

    // compute components and do tarjan for each component
    vec<connected_components::ConnectedComponent> components;
    connected_components::compute_components(hg, _context, components);

    for (const connected_components::ConnectedComponent component : components) {
        this->run_tarjan(hg, component);
    }

    // for each component compact regions

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(InitialPartitioningAlgorithm::random, _rng, _tag, time);
  }
}

template<typename TypeTraits>
void TarjanInitialPartitioner<TypeTraits>::run_tarjan(
    const PartitionedHypergraph& hypergraph, 
    const ConnetedComponent& component
) {
    
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(TarjanInitialPartitioner)

} // namespace mt_kahypar
