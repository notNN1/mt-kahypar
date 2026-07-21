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

#include "mt-kahypar/partition/connected_components/spanning_tree_connectivity.h"
#include <signal.h>

namespace mt_kahypar {
namespace connected_components {

using Bitset = mt_kahypar::ds::Bitset;

template<typename PartitionedHypergraph>
void BFSSpanningTreeConnectivity<PartitionedHypergraph>::reset(
    const PartitionedHypergraph& phg
) {
    DBG << "Initialized spanning tree on hypergraph instance: " << static_cast<const void*>(&phg) << " with " << phg.initialNumNodes() << " nodes";

    this->hn_to_num_children.clear();
    this->hn_to_num_children.resize(phg.initialNumNodes(), 0);
    this->hn_to_parent.clear();
    this->hn_to_parent.resize(phg.initialNumNodes());
    this->has_connection_to_other_partition.reset();
    this->has_connection_to_other_partition.resize(phg.initialNumNodes());
    this->hn_is_locked.reset();
    this->hn_is_locked.resize(phg.initialNumNodes());

    for (const HypernodeID& node : phg.nodes()) {
        this->hn_to_parent[node] = node;
    }

    // find nodes that have connections to other partitions
    Bitset edge_colored;
    edge_colored.resize(phg.initialNumEdges());

    for (const HypernodeID& hn : phg.nodes()) {            
        PartitionID current_partition = phg.partID(hn);

        for (const HyperedgeID& he : phg.incidentEdges(hn)) {

            if (edge_colored.isSet((size_t) he)) {
                continue;
            }

            edge_colored.set((size_t) he);

            for (const HypernodeID& incident_hn : phg.pins(he)) {

                if (current_partition != phg.partID(incident_hn)) {
                    this->has_connection_to_other_partition.set((size_t) hn);
                    break;
                }
            }

            if (this->has_connection_to_other_partition.isSet((size_t) hn)) {
                // now all the nodes have connections to different partitions as well

                for (const HypernodeID& incident_hn : phg.pins(he)) {
                    this->has_connection_to_other_partition.set((size_t) incident_hn);
                }
            }
        }
    }

    // do BFS with two queues

    Bitset node_colored;
    node_colored.resize(phg.initialNumNodes());
    
    edge_colored.reset();

    std::queue<HypernodeID> node_queue;
    std::queue<HypernodeID> colored_node_queue;
    PartitionID current_partition;

    for (const HypernodeID& hn : phg.nodes()) {
        if (node_colored.isSet((size_t) hn)) {
            continue;
        }

        current_partition = phg.partID(hn);

        if (this->has_connection_to_other_partition.isSet((size_t) hn)) {
            colored_node_queue.push(hn);
        }
        else {
            node_queue.push(hn);
        }

        node_colored.set((size_t) hn);
        
        edge_colored.reset();

        while (node_queue.size() > 0 || colored_node_queue.size() > 0) {
            
            HypernodeID current;

            if (node_queue.size() > 0) {
                current = node_queue.front();
                node_queue.pop();
            }
            else {
                current = colored_node_queue.front();
                colored_node_queue.pop();
            }


            for (const HyperedgeID& he : phg.incidentEdges(current)) {

                if (edge_colored.isSet((size_t) he)) {
                    continue;
                }

                edge_colored.set((size_t) he);

                for (const HypernodeID& incident_hn : phg.pins(he)) {
                    if (node_colored.isSet((size_t) incident_hn)) {
                        continue;
                    }
                    
                    if (phg.partID(incident_hn) != current_partition) {
                        continue;
                    }

                    if (incident_hn == current) {
                        continue;
                    }

                    node_colored.set((size_t) incident_hn);

                    if (this->has_connection_to_other_partition.isSet((size_t) incident_hn)) {
                        colored_node_queue.push(incident_hn);
                    }
                    else {
                        node_queue.push(incident_hn);
                    }

                    this->hn_to_num_children[current]++;
                    this->hn_to_parent[incident_hn] = current;
                }
            }
        }   
    }
}

template<typename PartitionedHypergraph>
bool BFSSpanningTreeConnectivity<PartitionedHypergraph>::canMoveVertex(
    const PartitionedHypergraph& phg, 
    const HypernodeID& hn
) {
    //LOG << "The size is wron: " << this->hn_to_num_children.size();
    //LOG << "HN: " << hn;
    if(hn >= this->hn_to_num_children.size()) {
        LOG << "hn nonexistant in array";
        raise(SIGSEGV);
    }

    if (this->hn_to_parent[hn] == hn || this->hn_is_locked.isSet((size_t) hn)) {
        return false;
    }

    return this->hn_to_num_children[hn] == 0;
}


template<typename PartitionedHypergraph>
void BFSSpanningTreeConnectivity<PartitionedHypergraph>::moveVertex(
    const PartitionedHypergraph& phg,
    const HypernodeID& hn,
    const PartitionID& to
) {
    vec<std::string> res = is_tree_valid(phg);
    if (res.size() != 0) {

        LOG << "Beginning of moveVertex";
        for (const std::string& err : res) {
            LOG << err;
        }

        raise(SIGTRAP); 
    }
    
    if (!canMoveVertex(phg, hn)) {
        LOG << "Move node without canMoveVertex being true";
        raise(SIGTRAP); 
    }

    if (phg.partID(hn) == to) {

        LOG << "Node is already in the selected partition";

        raise(SIGTRAP); 
    }

    this->hn_is_locked.set((size_t) hn);

    // find candidate to attach to
    for (const HyperedgeID& he : phg.incidentEdges(hn)) {
        for (const HypernodeID& incident_hn : phg.pins(he)) {
            if (phg.partID(incident_hn) != to) {
                continue;
            }

            bool can_move_to_colored_node = this->has_connection_to_other_partition.isSet((size_t) incident_hn) && this->hn_to_num_children[incident_hn] > 0;
            bool can_move_to_regular_node = !this->has_connection_to_other_partition.isSet((size_t) incident_hn);
            
            if (can_move_to_colored_node || can_move_to_regular_node) {
                HypernodeID parent = this->hn_to_parent[hn];
                this->hn_to_num_children[parent]--;

                this->hn_to_parent[hn] = incident_hn;
                this->hn_to_num_children[incident_hn]++;

                return;
            }
        }
    }

    // attach to the next best node
    for (const HyperedgeID& he : phg.incidentEdges(hn)) {
        for (const HypernodeID& incident_hn : phg.pins(he)) {
            if (phg.partID(incident_hn) != to) {
                continue;
            }

            HypernodeID parent = this->hn_to_parent[hn];
            this->hn_to_num_children[parent]--;

            this->hn_to_parent[hn] = incident_hn;
            this->hn_to_num_children[incident_hn]++;

            return;
        }
    }

    LOG << "There has not been found a node to attach to for node " << hn;
    while(true);
}

INSTANTIATE_CLASS_WITH_PARTITIONED_HG(BFSSpanningTreeConnectivity)

}  // namespace connected_components
}  // namespace mt_kahypar