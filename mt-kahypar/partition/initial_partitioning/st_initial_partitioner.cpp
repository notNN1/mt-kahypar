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
#include "mt-kahypar/partition/connected_components/tarjan.h"
#include <signal.h>
#include <algorithm>

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

        vec<size_t> covered;
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
        size_t current_split_number = 1;

        vec<HypernodeID> active_nodes;

        for (connected_components::ConnectedComponent& component : components) {
            current_split_number = 1;

            size_t upper_bound = (size_t)((double)target_size * (1.0 + 0.03));
            size_t lower_bound = (size_t)((double)target_size * (1.0 - 0.03));

            LOG << "Epsilon:     " << this->_context.partition.epsilon;
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

                if (current_size >= lower_bound && current_size <= upper_bound) {
                    
                    if (current_k + 1 < _context.partition.k) {
                        current_size = 0;
                        current_k++;
                    }
                    else {
                        // put remaining nodes into current partition
                        for (const HypernodeID& hn : component.nodes) {
                            if (covered[hn] > 0) {
                                continue;
                            }
                            hg.setNodePart(hn, current_k);
                        }

                        component_remaining_size = 0;
                    }
                }
                if (current_size + real_component_size > upper_bound) {
                    size_t target = upper_bound - current_size;
                    
                    LOG << "Target: " << target;

                    int max = 3;
                    int count = 0;

                    std::pair<HypernodeID, size_t> best_cut;

                    HypernodeID best_node   = -1;
                    size_t absolute_size    = 0;

                    while ((current_size + absolute_size < lower_bound && count < max) || count == 0) {
                        calculate_spanning_tree(hg, component, hn_to_parent, hn_to_children, subtree_size, covered, active_nodes);

                        best_cut = find_best_node_to_split(component, subtree_size, covered, target, current_split_number);

                        best_node           = best_cut.first;
                        absolute_size       = best_cut.second;
                        
                        count++;
                    }

                    LOG << "Split with absolute size " << absolute_size << " into partition " << current_k;
                    component_remaining_size    -= absolute_size;
                    current_size                += absolute_size;                    

                    assign_subtree_of_hn(hg, hn_to_children, subtree_size, current_k, covered, current_split_number, best_node);
                }
                else {
                    current_size                += real_component_size;
                    component_remaining_size    = 0;

                    for (const HypernodeID& hn : component.nodes) {
                        covered[hn] = current_split_number;
                        hg.setNodePart(hn, current_k);
                    }
                }
            }
        }

        vec<size_t> part_size;
        part_size.resize(32);
        for (const HypernodeID& node : hg.nodes()) {
            part_size[hg.partID(node)] += hg.nodeWeight(node);
        }

        int c = 0;
        for (const size_t& size : part_size) {            
            LOG << "Size part " << c << ": " << size;            
            c++;
        }

        
        HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
        double time = std::chrono::duration<double>(end - start).count();
        _ip_data.commit(InitialPartitioningAlgorithm::st, _rng, _tag, time);
    }
}

template<typename TypeTraits>
void STInitialPartitioner<TypeTraits>::calculate_spanning_tree(
    const PartitionedHypergraph& hg,
    ConnectedComponent& component,
    vec<HypernodeID>& hn_to_parent,
    vec<vec<HypernodeID>>& hn_to_children,
    vec<size_t>& subtree_size,
    vec<size_t>& covered,
    vec<HypernodeID>& active_nodes
) {
    assert(component.nodes.size() > 0);

    // find nodes that are in more than one edge
    Bitset is_in_multiple_edges;
    is_in_multiple_edges.resize(hg.initialNumNodes());
    for (const HypernodeID& node : component.nodes) {
        if (hg.incidentEdges(node).size() > 1) {
            is_in_multiple_edges.set((size_t) node);
        }
    }

    // setup datastructures for BFS
    std::deque<HypernodeID> queue;
    std::deque<HypernodeID> calculation_queue;

    Bitset node_colored;
    node_colored.resize(hg.initialNumNodes());

    // find first node, that is not covered
    
    std::shuffle(component.nodes.begin(), component.nodes.end(), _rng);

    for (const HypernodeID& node : component.nodes) {
        if (covered[node] == 0) {
            calculation_queue.push_back(node);
            queue.push_back(node);
            node_colored.set((size_t) node);
            break;
        }
    }


    Bitset edge_colored;
    edge_colored.resize(hg.initialNumEdges());

    size_t max_amount_of_branching = 2;

    
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

            size_t available_size = 1;
            for (const HypernodeID& incident_hn : hg.pins(he)) {
                if (node_colored.isSet((size_t) incident_hn) || covered[incident_hn] > 0) {
                    continue;
                }

                available_size++;
            }


            size_t branch_node_count = (available_size > max_amount_of_branching ? max_amount_of_branching : available_size);

            vec<HypernodeID> branch_nodes;
            branch_nodes.resize(branch_node_count);
            branch_nodes[0] = current_node;

            vec<std::pair<size_t, HypernodeID>> sizes;
            sizes.resize(branch_node_count);

            vec<std::pair<size_t, HypernodeID>> in_multiple_edges_sizes;
            in_multiple_edges_sizes.resize(branch_node_count);

            size_t target = available_size / branch_node_count;

            size_t found_branch_nodes = 1;

            // find branch nodes without multiple edges
            for (const HypernodeID& incident_hn : hg.pins(he)) {
                if (found_branch_nodes == branch_node_count) {
                    break;
                }

                if (node_colored.isSet((size_t) incident_hn) || covered[incident_hn] > 0) {
                    continue;
                }

                if (is_in_multiple_edges.isSet(incident_hn)) {
                    continue;
                }

                node_colored.set((size_t) incident_hn);

                branch_nodes[found_branch_nodes] = incident_hn;
                found_branch_nodes++;

                calculation_queue.push_back(incident_hn);

                hn_to_children[current_node].push_back(incident_hn);
                hn_to_parent[incident_hn] = current_node;
            }

            // find branch nodes with multiple edges
            for (const HypernodeID& incident_hn : hg.pins(he)) {
                if (found_branch_nodes == branch_node_count) {
                    break;
                }

                if (node_colored.isSet((size_t) incident_hn) || covered[incident_hn] > 0) {
                    continue;
                }

                if (!is_in_multiple_edges.isSet(incident_hn)) {
                    continue;
                }

                node_colored.set((size_t) incident_hn);

                branch_nodes[found_branch_nodes] = incident_hn;
                found_branch_nodes++;

                calculation_queue.push_back(incident_hn);

                hn_to_children[current_node].push_back(incident_hn);
                hn_to_parent[incident_hn] = current_node;
            }

            if (found_branch_nodes != branch_node_count) {
                LOG << "Branch node count: " << branch_node_count;
                LOG << "found_branch_nodes: " << found_branch_nodes;
                LOG << "available size: " << available_size;
                LOG << "rawr";
                while(true);
            }

            for (size_t i = 0; i < branch_node_count; i++) {
                sizes[i]                    = {0, branch_nodes[i]};
                in_multiple_edges_sizes[i]  = {0, branch_nodes[i]};
            }

            // distrbute nodes with multiple edges
            for (const HypernodeID& incident_hn : hg.pins(he)) {
                if (node_colored.isSet((size_t) incident_hn) || covered[incident_hn] > 0) {
                    continue;
                }

                if (!is_in_multiple_edges.isSet(incident_hn)) {
                    continue;
                }

                node_colored.set((size_t) incident_hn);

                std::sort(in_multiple_edges_sizes.begin(), in_multiple_edges_sizes.end());

                HypernodeID attachment_node = in_multiple_edges_sizes[0].second;
                in_multiple_edges_sizes[0] = {in_multiple_edges_sizes[0].first + hg.nodeWeight(incident_hn), attachment_node};


                queue.push_back(incident_hn);
                calculation_queue.push_back(incident_hn);

                hn_to_children[attachment_node].push_back(incident_hn);
                hn_to_parent[incident_hn] = attachment_node;
            }

            // distrbute nodes with multiple edges
            for (const HypernodeID& incident_hn : hg.pins(he)) {
                if (node_colored.isSet((size_t) incident_hn) || covered[incident_hn] > 0) {
                    continue;
                }

                if (!is_in_multiple_edges.isSet(incident_hn)) {
                    continue;
                }

                node_colored.set((size_t) incident_hn);

                std::sort(sizes.begin(), sizes.end());

                HypernodeID attachment_node = sizes[0].second;
                sizes[0] = {sizes[0].first + hg.nodeWeight(incident_hn), attachment_node};


                queue.push_back(incident_hn);
                calculation_queue.push_back(incident_hn);

                hn_to_children[attachment_node].push_back(incident_hn);
                hn_to_parent[incident_hn] = attachment_node;
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
    vec<size_t>& covered,
    const size_t& target,
    const size_t& current_split_number
) {
    assert(component.nodes.size() > 0);

    size_t best_size        = 0;
    HypernodeID best_node   = component.nodes[0];
    
    for (const HypernodeID& node : component.nodes) {
        if (covered[node] > 0) {
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
        vec<size_t>& covered,
        size_t& current_split_number,
        HypernodeID hn
) {
    std::queue<HypernodeID> queue;
    queue.push(hn);
    
    int count = 0;

    LOG << "partition: " << partition;

    while (queue.size() > 0) {
        HypernodeID current_node = queue.front();
        queue.pop();

        covered[current_node] = current_split_number;
        hg.setNodePart(current_node, partition);
        for (const HypernodeID& child : hn_to_children[current_node]) {
            if (child == current_node) {
                continue;
            }
            count++;

            queue.push(child);
        }
    }

    LOG << "Count: " << count;
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(STInitialPartitioner)

} // namespace mt_kahypar
