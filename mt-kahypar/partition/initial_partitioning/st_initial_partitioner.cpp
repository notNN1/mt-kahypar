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

const size_t MAX_CALCULATIONS = 10;

void check_tree(
    vec<HypernodeID>& hn_to_parent,
    const connected_components::ConnectedComponent& component
) {
    HypernodeID root = kInvalidHypernode;
    for (const HypernodeID& node : component.nodes) {
        
        HypernodeID parent = node;
        while(parent != hn_to_parent[parent]) {
            parent = hn_to_parent[parent];
        }

        if (root == kInvalidHypernode) {
            root = parent;
        }
        else if (root != parent) {
            LOG << "More than one root!";
            LOG << "Root:       " << root;
            LOG << "Other root: " << parent;
        }

    }
}

template<typename TypeTraits>
void STInitialPartitioner<TypeTraits>::partitionImpl() {
    if ( _ip_data.should_initial_partitioner_run(InitialPartitioningAlgorithm::st) ) {
        HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
        PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();

        
        ////// get components sorted by their size
        vec<connected_components::ConnectedComponent> components;
        connected_components::compute_components<typename TypeTraits::PartitionedHypergraph>(hg, _context, components);

        for (connected_components::ConnectedComponent& component : components) {
            std::shuffle(component.nodes.begin(), component.nodes.end(), _rng);
        }

        vec<std::pair<size_t, connected_components::ConnectedComponent>> components_and_size;

        for (const connected_components::ConnectedComponent& component : components) {

            size_t size = 0;
            for (const HypernodeID& node : component.nodes) {
                size += hg.nodeWeight(node);
            }

            components_and_size.push_back({size, component});
        }

        std::sort(components_and_size.begin(), components_and_size.end(),
        [](const auto& a, const auto& b) {
            return a.first > b.first;
        });


        ////// Setup important variables
        if (_context.partition.k != 2) {
            LOG << "k is not 2!!!";
            while (true);
        }

        size_t size_a = 0;
        size_t size_b = 0;

        size_t target = hg.totalWeight() / 2;

        vec<HypernodeID> nodes_to_swap;
        nodes_to_swap.reserve(hg.initialNumNodes());

        HypernodeID best_split;

        vec<HypernodeID> hn_to_parent;
        hn_to_parent.resize(hg.initialNumNodes());

        vec<vec<HypernodeID>> hn_to_children;
        hn_to_children.resize(hg.initialNumNodes());

        vec<size_t> subtree_size;
        subtree_size.resize(hg.initialNumNodes());

        vec<vec<HypernodeID>> best_hn_to_children;
        best_hn_to_children.reserve(hg.initialNumNodes());

        size_t current_origins = 0;

        ////// Split the components, only if there is not enough size
        for (const std::pair<size_t, connected_components::ConnectedComponent>& component_and_size : components_and_size) {

            best_split = kInvalidHypernode;
            nodes_to_swap.clear();
            current_origins = 0;

            size_t size                                         = component_and_size.first;
            connected_components::ConnectedComponent component  = component_and_size.second;


            if (size_a >= target) {
                for (const HypernodeID& node : component.nodes) {
                    hg.setNodePart(node, 1);
                    size_b += hg.nodeWeight(node);
                }
                continue;
            }


            if (size_a + size > target) { // split component

                size_t target_for_split = size - (target - size_a);

                for (const HypernodeID& node : component.nodes) {
                    hg.setNodePart(node, 0);
                    size_a += hg.nodeWeight(node);
                }

                //// calculate split 
                size_t split_size;
                size_t diff = 0;

                double best_split_diff = 1.0;

                do {
                    nodes_to_swap.clear();
                    split_size = 0;

                    calculate_spanning_tree(hg, component, hn_to_parent, hn_to_children, subtree_size); 

                    check_tree(hn_to_parent, component);

                    std::pair<HypernodeID, size_t> split = find_best_node_to_split(component, subtree_size, target_for_split * (1.0 + _context.partition.epsilon));

                    split_size = split.second;

                    LOG << "Split size: " << split_size;

                    diff = target_for_split >= split_size ? target_for_split - split_size : split_size - target_for_split;
                    current_origins++;

                    if (static_cast<double>(diff) / target_for_split < best_split_diff) {
                        best_split = split.first;
                        best_hn_to_children = hn_to_children;
                        best_split_diff = static_cast<double>(diff) / target_for_split;
                    }

                } while(current_origins < MAX_CALCULATIONS);

                LOG << "Target: " << target_for_split;

                //// assign nodes from best split            
                assign_subtree_of_hn(hg, best_hn_to_children, size_a, size_b, best_split);
                
            }
            else {
                for (const HypernodeID& node : component.nodes) {
                    hg.setNodePart(node, 0);
                    size_a += hg.nodeWeight(node);
                }
            }
            
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
    vec<size_t>& subtree_size
) {
    assert(component.nodes.size() > 0);

    // find nodes that are in more than one edge
    Bitset is_in_multiple_edges;

    if (_context.extra_options.st_options == STOptions::advanced_st_calculation) {
        is_in_multiple_edges.resize(hg.initialNumNodes());

        for (const HypernodeID& node : component.nodes) {
            if (hg.incidentEdges(node).size() > 1) {
                is_in_multiple_edges.set((size_t) node);
            }
        }
    }

    // setup datastructures for BFS
    std::deque<HypernodeID> queue;
    std::deque<HypernodeID> calculation_queue;

    Bitset node_colored;
    node_colored.resize(hg.initialNumNodes());

    const size_t start = std::uniform_int_distribution<size_t>(0, component.nodes.size() - 1)(_rng);
    HypernodeID starter_node = component.nodes[start];

    calculation_queue.push_back(starter_node);
    queue.push_back(starter_node);
    node_colored.set((size_t) starter_node);


    Bitset edge_colored;
    edge_colored.resize(hg.initialNumEdges());

    size_t max_amount_of_branching = 2;

    
    for (const HypernodeID& node : component.nodes) {
        subtree_size[node] = hg.nodeWeight(node);
        hn_to_children[node].clear();
        hn_to_parent[node] = node;
    }

    while (queue.size() > 0) {
        HypernodeID current_node = queue.back();
        queue.pop_back();

        for (const HyperedgeID& he : hg.incidentEdges(current_node)) {

            if (edge_colored.isSet((size_t) he)) {
                continue;
            }

            edge_colored.set((size_t) he);

            if (_context.extra_options.st_options == STOptions::normal_st) {
                for (const HypernodeID& incident_hn : hg.pins(he)) {
                    if (node_colored.isSet((size_t) incident_hn)) {
                        continue;
                    }

                    node_colored.set((size_t) incident_hn);

                    calculation_queue.push_back(incident_hn);
                    queue.push_back(incident_hn);

                    hn_to_children[current_node].push_back(incident_hn);
                    hn_to_parent[incident_hn] = current_node;
                }
            }
            else if (_context.extra_options.st_options == STOptions::advanced_st_calculation) {

                size_t available_size = 1;

                for (const HypernodeID& incident_hn : hg.pins(he)) {
                    if (node_colored.isSet((size_t) incident_hn)) {
                        continue;
                    }

                    ++available_size;
                }

                const size_t branch_node_count =
                    std::min(max_amount_of_branching, available_size);

                vec<HypernodeID> branch_nodes(branch_node_count);
                branch_nodes[0] = current_node;

                size_t found_branch_nodes = 1;

                for (const HypernodeID& incident_hn : hg.pins(he)) {
                    if (found_branch_nodes == branch_node_count) {
                        break;
                    }

                    if (node_colored.isSet((size_t) incident_hn)) {
                        continue;
                    }

                    if (is_in_multiple_edges.isSet(incident_hn)) {
                        continue;
                    }

                    node_colored.set((size_t) incident_hn);

                    branch_nodes[found_branch_nodes] = incident_hn;
                    ++found_branch_nodes;

                    calculation_queue.push_back(incident_hn);
                    queue.push_back(incident_hn);

                    hn_to_children[current_node].push_back(incident_hn);
                    hn_to_parent[incident_hn] = current_node;
                }

                for (const HypernodeID& incident_hn : hg.pins(he)) {
                    if (found_branch_nodes == branch_node_count) {
                        break;
                    }

                    if (node_colored.isSet((size_t) incident_hn)) {
                        continue;
                    }

                    if (!is_in_multiple_edges.isSet(incident_hn)) {
                        continue;
                    }

                    node_colored.set((size_t) incident_hn);

                    branch_nodes[found_branch_nodes] = incident_hn;
                    ++found_branch_nodes;

                    calculation_queue.push_back(incident_hn);
                    queue.push_back(incident_hn);

                    hn_to_children[current_node].push_back(incident_hn);
                    hn_to_parent[incident_hn] = current_node;
                }

                vec<std::pair<size_t, HypernodeID>> sizes;
                sizes.resize(branch_node_count);

                for (size_t i = 0; i < branch_node_count; ++i) {
                    sizes[i] = {
                        hg.nodeWeight(branch_nodes[i]),
                        branch_nodes[i]
                    };
                }

                for (const HypernodeID& incident_hn : hg.pins(he)) {
                    if (node_colored.isSet((size_t) incident_hn)) {
                        continue;
                    }

                    node_colored.set((size_t) incident_hn);

                    // Find the branch with the smallest current weight.
                    std::sort(sizes.begin(), sizes.end());

                    const HypernodeID attachment_node = sizes[0].second;

                    // Add this node to that branch.
                    sizes[0].first += hg.nodeWeight(incident_hn);

                    queue.push_back(incident_hn);
                    calculation_queue.push_back(incident_hn);

                    hn_to_children[attachment_node].push_back(incident_hn);
                    hn_to_parent[incident_hn] = attachment_node;
                }
            }
        }
    }

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
    const size_t& target
) {
    assert(component.nodes.size() > 0);

    size_t best_size        = 0;
    HypernodeID best_node   = kInvalidHypernode;
    
    for (const HypernodeID& node : component.nodes) {
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
    size_t& size_a,
    size_t& size_b,
    HypernodeID hn
) {
    std::queue<HypernodeID> queue;
    queue.push(hn);

    size_t total_size = 0;

    while (queue.size() > 0) {
        HypernodeID current_node = queue.front();
        queue.pop();

        hg.changeNodePart(current_node, 0, 1, DynamicConnectivityStrategy::do_nothing);

        size_a -= hg.nodeWeight(current_node);
        size_b += hg.nodeWeight(current_node);

        total_size += hg.nodeWeight(current_node);

        for (const HypernodeID& child : hn_to_children[current_node]) {
            if (child == current_node) {
                continue;
            }

            queue.push(child);
        }
    }

    LOG << "total size: " << total_size;
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(STInitialPartitioner)

} // namespace mt_kahypar
