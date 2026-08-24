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

#include "mt-kahypar/partition/connected_components/connectivity_facade.h"

namespace mt_kahypar {
namespace connected_components {

template<typename PartitionedHypergraph>
bool ConnectivityFacade<PartitionedHypergraph>::can_move_into_partition(
    const DynamicConnectivityStrategy& strategy,
    const PartitionedHypergraph& hypergraph,
    const HypernodeID& hn,
    const PartitionID& to
) {
    if (strategy == DynamicConnectivityStrategy::do_nothing) {
        return true;
    }

    return anker.find_node_in_partition(hypergraph, hn, to) != kInvalidHypernode;
}

template<typename PartitionedHypergraph>
bool ConnectivityFacade<PartitionedHypergraph>::can_move_out_of_partition(
    const DynamicConnectivityStrategy& strategy,
    const PartitionedHypergraph& hypergraph,
    const HypernodeID& hn,
    const PartitionID& to
) {
    bool can_move_node = true;

    if (strategy == DynamicConnectivityStrategy::bfs) {
        can_move_node = this-bfs.moveVertex(hypergraph, hn);
    }
    else if (strategy == DynamicConnectivityStrategy::h_vertex_degree) {
        auto range = hypergraph.incidentEdges(hn);
        can_move_node = std::distance(range.begin(), range.end()) > 4;
    }
    else if (strategy == DynamicConnectivityStrategy::st) {
        can_move_node = this->stc.canMoveVertex(hypergraph, hn);
    }

    return can_move_node;
}

template<typename PartitionedHypergraph>
void ConnectivityFacade<PartitionedHypergraph>::moveVertex(
    const DynamicConnectivityStrategy& strategy,
    const PartitionedHypergraph& hypergraph,
    const HypernodeID& hn,
    const PartitionID& to
) {

    HypernodeID node_to = this->anker.find_node_in_partition(hypergraph, hn, to);

    if (strategy == DynamicConnectivityStrategy::st) {
        this->stc.moveVertex(hypergraph, hn, to, node_to);
    }

    this->last_strategy_used = strategy;
}

template<typename PartitionedHypergraph>
bool ConnectivityFacade<PartitionedHypergraph>::canMoveVertex(
    const DynamicConnectivityStrategy& strategy,
    const PartitionedHypergraph& hypergraph,
    const HypernodeID& hn,
    const PartitionID& to
) {
    if (strategy != last_strategy_used) {
        this->reset_connectivity(hypergraph, strategy);
    }

    return can_move_out_of_partition(strategy, hypergraph, hn, to) && can_move_into_partition(strategy, hypergraph, hn, to);
}

INSTANTIATE_CLASS_WITH_PARTITIONED_HG(ConnectivityFacade)

}  // namespace connected_components
}  // namespace mt_kahypar
