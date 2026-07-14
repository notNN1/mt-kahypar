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
void ConnectivityFacade<PartitionedHypergraph>::reset_connectivity_in(
    const PartitionedHypergraph& hypergraph,
    const HypernodeID& hn
) {
    this->can_move_current_node_to_partition.reset();
    this->can_move_current_node_to_partition.resize(hypergraph.k());

    
    for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        for (const PartitionID& partition : hypergraph.connectivitySet(he)) {
            this->can_move_current_node_to_partition.set((size_t) partition); 
        }
    }

    this->current_node = hn;
}

template<typename PartitionedHypergraph>
void ConnectivityFacade<PartitionedHypergraph>::reset_connectivity_out(
    const PartitionedHypergraph& hypergraph,
    const HypernodeID& hn
) {
    this->stc.reset(hypergraph, _context);
}


template<typename PartitionedHypergraph>
bool ConnectivityFacade<PartitionedHypergraph>::can_move_into_partition(
    const DynamicConnectivityStrategy& strategy,
    const PartitionedHypergraph& hypergraph,
    const HypernodeID& hn,
    const PartitionID& to
) {
    if (
        strategy == DynamicConnectivityStrategy::bfs
        || strategy == DynamicConnectivityStrategy::h_vertex_degree
        || strategy == DynamicConnectivityStrategy::st
    ) {
        if (hn != this->current_node) {
            this->initialize_can_move_current_node_to_partition(hypergraph, strategy, hn);
        }

        return this->can_move_current_node_to_partition.isSet((size_t) to);
    }
    return true; // default
}

template<typename PartitionedHypergraph>
bool ConnectivityFacade<PartitionedHypergraph>::can_move_out_of_partition(
    const DynamicConnectivityStrategy& strategy,
    const PartitionedHypergraph& hypergraph,
    const HypernodeID& hn
) {
    bool can_move_node = true;

    if (strategy == DynamicConnectivityStrategy::bfs) {
        mt_kahypar::ds::BFSConnectivity<PartitionedHypergraph> dcd = mt_kahypar::ds::BFSConnectivity<PartitionedHypergraph>();
        can_move_node = dcd.moveVertex(hypergraph, hn);
    }
    else if (strategy == DynamicConnectivityStrategy::h_vertex_degree) {
        auto range = hypergraph.incidentEdges(hn);
        can_move_node = std::distance(range.begin(), range.end()) > 4;
    }
    else if (strategy == DynamicConnectivityStrategy::st) {
        can_move_node = this->stc.canMoveVertex(hypergraph, _context, hn);
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
    if (strategy == DynamicConnectivityStrategy::st) {
        this->stc.moveVertex();
    }

    this->graph_was_changed = true;
}

}  // namespace connected_components
}  // namespace mt_kahypar
