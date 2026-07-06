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
#include <list>

namespace mt_kahypar {
namespace ds {

using Bitset = mt_kahypar::ds::Bitset;
using ConnectedComponent = mt_kahypar::connected_components::ConnectedComponent;
using ComponentID = uint32_t;

template<typename PartitionedHypergraph>
class BFSConnectivity {
public:
    bool moveVertex(
        const PartitionedHypergraph& phg, 
        const HypernodeID& hn
    );
};

template<typename PartitionedHypergraph>
class BFSSpanningTreeConnectivity {
private:
    vec<HypernodeID> hn_to_parent;
    vec<size_t> hn_to_num_children;

    Bitset has_connection_to_other_partition;
public:
    BFSSpanningTreeConnectivity(const PartitionedHypergraph& phg,
                                const Context& context);
    
    int size() const {
        return hn_to_num_children.size();
    } 

    vec<std::string> is_tree_valid(const PartitionedHypergraph& phg) {
        vec<std::string> result;
        if (wrong_size_of_arrays(phg)) {
            result.push_back("Wrong size of arrays");
        } 
        if (!all_nodes_go_to_one_of_the_heads(phg)) {
            result.push_back("At least one node doesn't reach the head");
        } 
        if (!connections_exist(phg)) {
            result.push_back("At least one connection does not exist");
        } 
        if (!nodes_in_correct_partition(phg)) {
            result.push_back("At least one node is not in the correct partition");
        }
        if (!nodes_have_correct_amount_of_children(phg)) {
            result.push_back("Incorrect amound of children detected");
        }

        return result;
    }

    bool nodes_have_correct_amount_of_children(const PartitionedHypergraph& phg) {
        vec<size_t> children_alternative_count;
        children_alternative_count.resize(phg.initialNumNodes());

        
        for (const HypernodeID& node : phg.nodes()) {
            if (this->hn_to_parent[node] != node) {
                children_alternative_count[this->hn_to_parent[node]]++;
            }
        }

        for (const HypernodeID& node : phg.nodes()) {
            if (children_alternative_count[node] != this->hn_to_num_children[node]) {
                return false;
            }
        }

        return true;
    }

    bool nodes_in_correct_partition(const PartitionedHypergraph& phg) {
        for (const HypernodeID& node : phg.nodes()) {
            if (phg.partID(node) != phg.partID(this->hn_to_parent[node])) {
                LOG << "Node not in correct parititon";
                LOG << "Node:   " << node << " in partition " << phg.partID(node);
                LOG << "Parent: " << this->hn_to_parent[node] << " in partition " << phg.partID(this->hn_to_parent[node]);
                return false;
            }
        }

        return true;
    }

    bool connections_exist(const PartitionedHypergraph& phg) {
        bool found = false;
        for (const HypernodeID& node : phg.nodes()) {
            for (const HyperedgeID& edge : phg.incidentEdges(node)) {
                for (const HypernodeID& incident_hn : phg.pins(edge)) {
                    if (incident_hn == this->hn_to_parent[incident_hn]) {
                        found = true;
                    }
                }
            }
        }

        return found;
    }

    bool wrong_size_of_arrays(const PartitionedHypergraph& phg) {
        return this->hn_to_parent.size() != phg.initialNumNodes() || this->hn_to_num_children.size() != phg.initialNumNodes();
    }

    bool all_nodes_go_to_one_of_the_heads(const PartitionedHypergraph& phg) {
        vec<HypernodeID> heads;

        for (const HypernodeID& node : phg.nodes()) {
            if (this->hn_to_parent[node] == node) {
                heads.push_back(node);
            }
        }

        for (const HypernodeID& node : phg.nodes()) {

            HypernodeID parent = node;

            while (this->hn_to_parent[parent] != parent) {
                parent = this->hn_to_parent[parent];
            }

            bool found = false;
            
            for (const HypernodeID& head : heads) {
                if (head == parent) {
                    found = true;
                }
            }

            if (!found) {
                return false;
            }

        }

        return true;
    }

    bool canMoveVertex(
        const PartitionedHypergraph& phg,
        const Context& context,
        const HypernodeID& hn
    );

    void moveVertex(
        const PartitionedHypergraph& phg,
        const Context& context,
        const HypernodeID& hn,
        const PartitionID& to
    );
};

template<typename PartitionedHypergraph>
class SpanningTreeConnectivity {
private:
    
    vec<HypernodeID> vertex_to_parent_compressed;
    vec<size_t> vertex_to_rank;

    struct Connection {
        typename std::list<Connection>::iterator iterator;
        HypernodeID node;

        Connection(HypernodeID node_) : node(node_) {}
    };

    std::vector<std::list<Connection>> connected_to;
    

    inline void try_connect_to_incident_without_connection(
        Bitset& has_connection_to_other_partition,
        const HypernodeID& hn,
        const HypernodeID& incident_hn
    );

    inline void try_connect_to_incident_with_connection(
        Bitset& has_connection_to_other_partition,
        const HypernodeID& hn,
        const HypernodeID& incident_hn
    );

    inline void connect_nodes(
        const HypernodeID& hn,
        const HypernodeID& incident_hn
    );

    inline bool is_same_component(
        const HypernodeID& hn1,
        const HypernodeID& hn2
    );

    int calc_parents(const PartitionedHypergraph& phg);
public:
    SpanningTreeConnectivity(
        const PartitionedHypergraph& phg,
        const Context& context
    );

    SpanningTreeConnectivity() = default;

    bool canMoveVertex(
        const Context& context,
        const HypernodeID& hn
    );

    void moveVertex(
        const PartitionedHypergraph& phg,
        const HypernodeID& hn,
        const PartitionID& to
    );
};

template<typename PartitionedHypergraph>
class ConnectivityFacade {
private:
    Bitset          can_move_current_node_to_partition;
    HypernodeID     current_node;
    int size;

    void initialize_can_move_current_node_to_partition(
        const PartitionedHypergraph& hypergraph,
        const DynamicConnectivityStrategy& strategy,
        const HypernodeID& hn
    );

public:
    static int extra_component_count(
        const PartitionedHypergraph& hypergraph,
        const Context& _context
    );

    bool can_move_out_of_partition(
        const DynamicConnectivityStrategy& strategy,
        const PartitionedHypergraph& hypergraph,
        const Context& _context,
        const HypernodeID& hn
    );

    void move_out_of_partition(
        const DynamicConnectivityStrategy& strategy,
        const PartitionedHypergraph& hypergraph,
        const HypernodeID& hn,
        const PartitionID& to
    );

    bool can_move_into_partition(
        const DynamicConnectivityStrategy& strategy,
        const PartitionedHypergraph& hypergraph,
        const Context& _context,
        const HypernodeID& hn,
        const PartitionID& to
    );
};

}  // namespace ds
}  // namespace mt_kahypar
