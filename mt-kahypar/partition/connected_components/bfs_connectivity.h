/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
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

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/connected_components/compute_components.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/bitset.h"
#include "mt-kahypar/partition/connected_components/tarjan.h"

namespace mt_kahypar {
namespace connected_components {

// This implementation is very wrong, because tarjan does not check for partitions. 
// Tarjan was implemented with the tarjan initial partitioner in mind
template<typename PartitionedHypergraph>
class BFSConnectivity {
public:
    bool moveVertex(
        const PartitionedHypergraph& phg, 
        const HypernodeID& hn
    ) {
        int32_t component_count = connected_components::total_component_count(phg, NULL);

        vec<vec<ConnectedComponent>> connected_components_per_block_colored;
        
        Bitset node_colored;
        node_colored.resize(phg.initialNumNodes());
        node_colored.set((size_t) hn);

        connected_components::compute_components_per_block_colored(phg, NULL, connected_components_per_block_colored, node_colored);

        int32_t count_induced_graph = 0;
        for (vec<ConnectedComponent> block : connected_components_per_block_colored) {
            count_induced_graph += block.size();
        }

        return count_induced_graph <= component_count;
    };
};

}  // namespace connected_components
}  // namespace mt_kahypar