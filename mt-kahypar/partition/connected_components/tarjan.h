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
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/connected_components/compute_components.h"
#include "mt-kahypar/datastructures/bitset.h"

#include <algorithm>
#include <set>
#include <utility>
#include <signal.h>

namespace mt_kahypar {
namespace connected_components {

using Bitset = mt_kahypar::ds::Bitset;
using Time   = size_t;

template<typename PartitionedHypergraph>
struct Tarjan {
private:

    vec<size_t>     disc;
    vec<size_t>     low;
    
    Time            time;

    Bitset          is_ap;
    vec<vec<HypernodeID>> biconnected_components;
    vec<HypernodeID> vertex_stack;
    bool            initialized = false;

public:
    bool is_articulation_point(
        const PartitionedHypergraph& phg,
        const HypernodeID& node
    ) {
        if (!initialized) {
            initialize(phg);
        }

        return is_ap.isSet((size_t) node);
    }

    bool is_articulation_point(
        const PartitionedHypergraph& phg,
        const HypernodeID& node,
        const Bitset& excluded
    ) {
        initialize(phg, &excluded);
        return is_ap.isSet((size_t) node);
    }

    // Components of the circular expansion. Articulation vertices may belong to
    // multiple components; bridges form pairs and isolated vertices singletons.
    vec<vec<HypernodeID>> get_biconnected_components(const PartitionedHypergraph& phg) {
        if (!initialized) {
            initialize(phg);
        }
        return biconnected_components;
    }

    void initialize(const PartitionedHypergraph& phg, const Bitset* excluded = nullptr) {
        reset_fields(phg);
        vec<vec<HypernodeID>> hn_to_new_incident_nodes;
        
        if (excluded) {
            hn_to_new_incident_nodes.resize(phg.initialNumNodes());
            vec<HypernodeID> remaining_pins;

            for (const HyperedgeID he : phg.edges()) {

                remaining_pins.clear();

                for (const HypernodeID pin : phg.pins(he)) {
                    if (!excluded->isSet((size_t) pin)) {
                        remaining_pins.push_back(pin);
                    }
                }

                // Reconnect surviving pins
                const size_t size = remaining_pins.size();

                for (size_t i = 0; i < size; ++i) {

                    const HypernodeID current = remaining_pins[i];
                    const HypernodeID prev = remaining_pins[(i + size - 1) % size];
                    const HypernodeID next = remaining_pins[(i + 1) % size];

                    if (prev != current) {
                        hn_to_new_incident_nodes[current].push_back(prev);
                    }

                    if (next != current && next != prev) {
                        hn_to_new_incident_nodes[current].push_back(next);
                    }
                }
            }
        } else {
            hn_to_new_incident_nodes = phg.get_circular_edge_expansion();
        }
        
        for (const HypernodeID& node : phg.nodes()) {
            if (excluded && excluded->isSet((size_t) node)) {
                continue;
            }
            if (disc[node] == 0) {
                const size_t children = dfs_recursive_wrapper(phg, node, hn_to_new_incident_nodes);
                if (children > 1) {
                    is_ap.set((size_t) node);
                } else if (children == 0) {
                    biconnected_components.push_back({node});
                }
            }
        }
        initialized = true;
    }
private:
    void extract_biconnected_component(const HypernodeID from, const HypernodeID to) {
        vec<HypernodeID> component;

        // Pop the child's remaining subtree, including the child itself.
        // Keep the parent on the stack: it may belong to another component.
        while (!vertex_stack.empty()) {
            const HypernodeID node = vertex_stack.back();
            vertex_stack.pop_back();
            component.push_back(node);
            if (node == to) {
                break;
            }
        }
        component.push_back(from);
        biconnected_components.push_back(std::move(component));
    }

    size_t dfs_recursive_wrapper(
        const PartitionedHypergraph& phg, 
        const HypernodeID& node,
        const vec<vec<HypernodeID>>& hn_to_new_incident_nodes
    ) {
        const size_t children = dfs_recursive(phg, node, node, hn_to_new_incident_nodes);
        // Every child component has been extracted; only the DFS root remains.
        vertex_stack.pop_back();
        return children;
    }

    size_t dfs_recursive(
        const PartitionedHypergraph& phg, 
        const HypernodeID& node, 
        const HypernodeID& parent,
        const vec<vec<HypernodeID>>& hn_to_new_incident_nodes
    ) {
        size_t children = 0;
        low[node] = disc[node] = time;
        time++;
        vertex_stack.push_back(node);

        bool skipped_parent = false;
        for (const HypernodeID incident_hn : hn_to_new_incident_nodes[node]) {
            // Skip the incoming edge once; parallel edges are valid back edges.
            if (incident_hn == parent && !skipped_parent) {
                skipped_parent = true;
                continue;
            }
            if (disc[incident_hn] == 0) { // if V has not been discovered before
                children++;
                dfs_recursive(phg, incident_hn, node, hn_to_new_incident_nodes); // recursive DFS call

                if (disc[node] <= low[incident_hn]) {
                    if (node != parent) {
                        is_ap.set((size_t) node);
                    }
                    extract_biconnected_component(node, incident_hn);
                }

                low[node] = std::min(low[node], low[incident_hn]); // low[v] might be an ancestor of u
            } else if (disc[incident_hn] < disc[node]) { // back edge to an ancestor
                low[node] = std::min(low[node], disc[incident_hn]); // finds the ancestor with the least discovery time
            }
        }
        
        return children;
    }


    size_t dfs(
        const PartitionedHypergraph& phg, 
        const HypernodeID& node
    ) {
        size_t children = 0;

        vec<std::pair<HypernodeID, HypernodeID>> queue;
        queue.push_back({node, node});

        vec<std::pair<HypernodeID, HypernodeID>> bt_queue;
        bt_queue.push_back({node, node});
        
        HypernodeID current;
        HypernodeID parent;

        Bitset node_colored;
        node_colored.resize(phg.initialNumNodes());

        bool backtracking       = false;


        vec<vec<vec<HypernodeID>>> he_hn_to_new_edge = phg.circular_edge_expansion();

        while(queue.size() > 0) {
            
            current     = queue.back().first;
            parent      = queue.back().second;

            queue.pop_back();

            if (!backtracking) {
                low[current]    = time;
                disc[current]   = time;
                time++;

                backtracking = true;
                
                for (const HyperedgeID& he : phg.incidentEdges(current)) {

                    for (const HypernodeID& incident_hn : he_hn_to_new_edge[he][current]) {

                        if (incident_hn == parent) {
                            continue;
                        }

                        if (disc[incident_hn] > 0) {
                            low[current] = low[current] < disc[incident_hn] ? low[current] : disc[incident_hn];
                            continue;
                        }

                        if (current == node) {
                            children++;
                        }


                        node_colored.set((size_t) incident_hn);

                        queue.push_back({incident_hn, current});
                        bt_queue.push_back({incident_hn, current});
                        
                        backtracking = false;
                    }
                }
            }
            else {
                HypernodeID current_bt = bt_queue.back().first;
                HypernodeID parent_bt  = bt_queue.back().second;

                while(current != current_bt) {
                    if (bt_queue.empty()) {
                        LOG << "empty";
                        raise(SIGSEGV);
                    }

                    bt_queue.pop_back();
                    current_bt  = bt_queue.back().first;
                    parent_bt   = bt_queue.back().second;

                    if (disc[parent_bt] <= low[current_bt]) {
                        is_ap.set((size_t) parent_bt);
                    }

                    low[parent_bt] = low[parent_bt] < low[current_bt] ? low[parent_bt] : low[current_bt];
                }

                backtracking = false;
            }
            
        }

        return children;
    }

    void reset_fields(const PartitionedHypergraph& phg) {
        initialized = false;
        biconnected_components.clear();
        vertex_stack.clear();
        disc.clear();
        disc.resize(phg.initialNumNodes());

        low.clear();
        low.resize(phg.initialNumNodes());

        is_ap.reset();
        is_ap.resize(phg.initialNumNodes());

        time = 1;
    }


};

}  // namespace connected_components
}  // namespace mt_kahypar
