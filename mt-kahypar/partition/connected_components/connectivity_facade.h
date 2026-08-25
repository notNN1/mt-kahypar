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
#include "mt-kahypar/datastructures/bitset.h"
#include "mt-kahypar/partition/connected_components/compute_components.h"
#include "mt-kahypar/partition/connected_components/connectivity_facade.h"
#include "mt-kahypar/partition/connected_components/spanning_tree_connectivity.h"
#include "mt-kahypar/partition/connected_components/bfs_connectivity.h"
#include "mt-kahypar/partition/connected_components/anker_node.h"



namespace mt_kahypar {
namespace connected_components {

using Bitset = mt_kahypar::ds::Bitset;
using ConnectedComponent = mt_kahypar::connected_components::ConnectedComponent;

template<typename PartitionedHypergraph>
class ConnectivityFacade {
private:
    BFSSpanningTreeConnectivity<PartitionedHypergraph>  stc;                            // reset when: moved without connectivity, uncoarsening 
    BFSConnectivity<PartitionedHypergraph>              bfs;

    AnkerNodes<PartitionedHypergraph>                   anker;
    
    DynamicConnectivityStrategy                         last_strategy_used;

    bool can_move_out_of_partition(
        const DynamicConnectivityStrategy& strategy,
        const PartitionedHypergraph& hypergraph,
        const HypernodeID& hn,
        const PartitionID& to
    );

    bool can_move_into_partition(
        const DynamicConnectivityStrategy& strategy,
        const PartitionedHypergraph& hypergraph,
        const HypernodeID& hn,
        const PartitionID& to
    );

public:
    void initialize(const PartitionedHypergraph& hypergraph) {
        this->stc.reset(hypergraph);
    }

    bool canMoveVertex(
        const DynamicConnectivityStrategy& strategy,
        const PartitionedHypergraph& hypergraph,
        const HypernodeID& hn,
        const PartitionID& to
    );

    void moveVertex(
        const DynamicConnectivityStrategy& strategy,
        const PartitionedHypergraph& hypergraph,
        const HypernodeID& hn,
        const PartitionID& to
    );

    void reset_connectivity(
        const PartitionedHypergraph& hypergraph,
        const DynamicConnectivityStrategy& strategy
    );

    void reset_connectivity(
        const PartitionedHypergraph& hypergraph
    );

};

}  // namespace connected_components
}  // namespace mt_kahypar