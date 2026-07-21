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

#include "mt-kahypar/partition/initial_partitioning/st_initial_partitioner.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/randomize.h"
#include <signal.h>

namespace mt_kahypar {

template<typename TypeTraits>
void STInitialPartitioner<TypeTraits>::partitionImpl() {
    if ( _ip_data.should_initial_partitioner_run(InitialPartitioningAlgorithm::st) ) {
        HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
        PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();

        // compute components and calculate st for each component
        vec<connected_components::ConnectedComponent> components;
        connected_components::compute_components<typename TypeTraits::PartitionedHypergraph>(hg, _context, components);

        vec<vec<HypernodeID>> hn_to_children;
        hn_to_children.resize(hg.initialNumNodes());

        vec<HypernodeID> hn_to_parent;
        hn_to_parent.resize(hg.initialNumNodes());

        for (const HypernodeID& hn : hg.nodes()) {
            hn_to_parent[hn] = hn;
        }

        Bitset covered;
        covered.resize(hg.initialNumNodes());

        vec<size_t> subtree_size(hg.initialNumNodes());
        size_t total_size = 0;
        for (const HypernodeID& hn : hg.nodes()) {
            subtree_size[hn] = hg.nodeWeight(hn);
            total_size       += hg.nodeWeight(hn);
        }

        size_t target_size = (total_size / _context.partition.k);


        size_t current_size = 0;
        size_t current_k = 0;
        size_t component_remaining_size;
        size_t real_component_size;


        for (const connected_components::ConnectedComponent& component : components) {
            size_t upper_bound = (size_t)((double)target_size * (1.0 + this->_context.partition.epsilon));
            size_t lower_bound = (size_t)((double)target_size * (1.0 - this->_context.partition.epsilon)) + 1;

            LOG << "Target size: " << target_size;
            LOG << "Upper bound: " << upper_bound;
            LOG << "Lower bound: " << lower_bound;

            real_component_size         = 0;
            component_remaining_size    = 0;
            for (const HypernodeID& node : component.nodes) {
                component_remaining_size    += hg.nodeWeight(node);
                real_component_size         += hg.nodeWeight(node);
            }


            while (component_remaining_size > 0) {

                if (current_size > lower_bound && current_size <= upper_bound) {
                    
                    if (current_k + 1 < _context.partition.k) {
                        current_size = 0;
                        current_k++;
                    }
                    else {
                        // put remaining nodes into current partition
                        for (const HypernodeID& hn : component.nodes) {
                            if (covered.isSet((size_t) hn)) {
                                continue;
                            }
                            hg.setNodePart(hn, current_k);
                        }

                        component_remaining_size = 0;
                    }
                }
                if (current_size + real_component_size > upper_bound) {
                    size_t target = upper_bound - current_size;
                    
                    int max = 3;
                    int count = 0;

                    std::pair<HypernodeID, size_t> best_cut;

                    HypernodeID best_node   = -1;
                    size_t absolute_size    = 0;

                    while ((current_size + absolute_size < lower_bound && count < max) || count == 0) {
                        calculate_spanning_tree(hg, component, hn_to_parent, hn_to_children, subtree_size, covered);

                        best_cut = find_best_node_to_split(component, subtree_size, covered, target);

                        best_node           = best_cut.first;
                        absolute_size       = best_cut.second;
                        
                        count++;
                    }

                    LOG << "Split with absolute size " << absolute_size << " into partition " << current_k;
                    component_remaining_size    -= absolute_size;
                    current_size                += absolute_size;                    

                    assign_subtree_of_hn(hg, hn_to_children, subtree_size, current_k, covered, best_node);
                }
                else {
                    current_size                += real_component_size;
                    component_remaining_size    = 0;

                    for (const HypernodeID& hn : component.nodes) {
                        covered.set((size_t) hn);
                        hg.setNodePart(hn, current_k);
                    }
                }
            }
        }

        vec<size_t> part_size;
        part_size.resize(2);
        for (const HypernodeID& node : hg.nodes()) {
            if (hg.partID(node) == 0) {
                part_size[0] += hg.nodeWeight(node);
            }
            else {
                part_size[1] += hg.nodeWeight(node);
            }
        }

        LOG << "Size part 0: " << part_size[0];
        LOG << "Size part 1: " << part_size[1];
        
        HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
        double time = std::chrono::duration<double>(end - start).count();
        _ip_data.commit(InitialPartitioningAlgorithm::st, _rng, _tag, time);
    }
}

template<typename TypeTraits>
void STInitialPartitioner<TypeTraits>::calculate_spanning_tree(
    const PartitionedHypergraph& hg,
    const ConnectedComponent& component,
    vec<HypernodeID>& hn_to_parent,
    vec<vec<HypernodeID>>& hn_to_children,
    vec<size_t>& subtree_size,
    Bitset& covered
) {
    assert(component.nodes.size() > 0);

    // setup datastructures for BFS
    std::deque<HypernodeID> queue;
    queue.push_back(component.nodes[0]);

    std::deque<HypernodeID> calculation_queue;

    // find first node, that is not covered
    for (const HypernodeID& node : component.nodes) {
        if (!covered.isSet((size_t) node)) {
            calculation_queue.push_back(node);
            break;
        }
    }

    Bitset node_colored;
    node_colored.resize(hg.initialNumNodes());
    node_colored.set((size_t) component.nodes[0]);

    Bitset edge_colored;
    edge_colored.resize(hg.initialNumEdges());

    
    for (const HypernodeID& node : component.nodes) {
        subtree_size[node] = hg.nodeWeight(node);
        hn_to_children[node].clear();
        hn_to_parent[node] = node;
    }

    while (queue.size() > 0) {
        HypernodeID current_node    = queue.back();
        queue.pop_back();

        for (const HyperedgeID& he : hg.incidentEdges(current_node)) {

            if (edge_colored.isSet((size_t) he)) {
                continue;
            }

            edge_colored.set((size_t) he);

            for (const HypernodeID& incident_hn : hg.pins(he)) {
                if (node_colored.isSet((size_t) incident_hn) || covered.isSet((size_t) incident_hn)) {
                    continue;
                }

                node_colored.set((size_t) incident_hn);

                queue.push_back(incident_hn);
                calculation_queue.push_back(incident_hn);

                hn_to_children[current_node].push_back(incident_hn);
                hn_to_parent[incident_hn] = current_node;
            }
        }
    }
    
    size_t best_size        = 0;
    HypernodeID best_node   = component.nodes[0];

    // calculate size of all subtrees and find
    while (calculation_queue.size() > 0) {
        HypernodeID current = calculation_queue.back();
        calculation_queue.pop_back();

        if (hn_to_parent[current] == current) {
            continue;
        }

        subtree_size[hn_to_parent[current]] += subtree_size[current];
    }
}

template<typename TypeTraits>
std::pair<HypernodeID, size_t> STInitialPartitioner<TypeTraits>::find_best_node_to_split(
    const ConnectedComponent& component,
    const vec<size_t>& subtree_size,
    const Bitset& covered,
    const size_t& target
) {
    assert(component.nodes.size() > 0);

    size_t best_size        = 0;
    HypernodeID best_node   = component.nodes[0];
    
    for (const HypernodeID& node : component.nodes) {
        if (covered.isSet((size_t) node)) {
            continue;
        }

        if (subtree_size[node] > best_size && subtree_size[node] <= target) {
            best_size = subtree_size[node];
            best_node = node;
        }
    }

    return {best_node, best_size};
}

template<typename TypeTraits>
void STInitialPartitioner<TypeTraits>::assign_subtree_of_hn(
    PartitionedHypergraph& hg,
    vec<vec<HypernodeID>>& hn_to_children,
    vec<size_t>& subtree_size,
    PartitionID partition,
    Bitset& covered,
    HypernodeID hn
) {
    std::queue<HypernodeID> queue;
    queue.push(hn);
    
    while (queue.size() > 0) {
        HypernodeID current_node = queue.front();
        queue.pop();

        covered.set((size_t) current_node);
        hg.setNodePart(current_node, partition);
        for (const HypernodeID& child : hn_to_children[current_node]) {
            if (child == current_node) {
                continue;
            }

            queue.push(child);
        }
    }
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(STInitialPartitioner)

} // namespace mt_kahypar
