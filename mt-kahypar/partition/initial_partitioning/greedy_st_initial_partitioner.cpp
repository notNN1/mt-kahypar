
/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/initial_partitioning/greedy_st_initial_partitioner.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/randomize.h"

#include <deque>
#include <algorithm>



namespace mt_kahypar {

template<typename TypeTraits>
void GreedySTInitialPartitioner<TypeTraits>::partitionImpl() {
  if ( _ip_data.should_initial_partitioner_run(InitialPartitioningAlgorithm::st) ) {
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();
    std::uniform_int_distribution<PartitionID> select_random_block(0, _context.partition.k - 1);

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

    ////// Split the components, only if there is not enough size
    for (const std::pair<size_t, connected_components::ConnectedComponent>& component_and_size : components_and_size) {

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
            size_t max_splits = 3;
            size_t cur_splits = 0;

            size_t split_size;
            vec<HypernodeID> nodes_to_swap;
            size_t diff = 0;

            double best_split_diff = 1.0;
            vec<HypernodeID> best_split;

            do {
                nodes_to_swap.clear();
                split_size = 0;

                calculate_split(hg, component, target_for_split, nodes_to_swap, split_size);

                diff = target_for_split >= split_size ? target_for_split - split_size : split_size - target_for_split;
                cur_splits++;

                if (static_cast<double>(diff) / target_for_split < best_split_diff) {
                    best_split = nodes_to_swap;
                    best_split_diff = static_cast<double>(diff) / target_for_split;
                }

            } while(cur_splits <= max_splits);

            ////

            LOG << "Best Split diff: " << best_split_diff;

            for (const HypernodeID& node : best_split) {
                hg.changeNodePart(node, 0, 1, false);
                size_a -= hg.nodeWeight(node);
                size_b += hg.nodeWeight(node);
            }
            
        }
        else {
            for (const HypernodeID& node : component.nodes) {
                hg.setNodePart(node, 0);
                size_a += hg.nodeWeight(node);
            }
        }
        
    }
    

    if (connected_components::total_component_count(hg, _context) != 2) {
        LOG << "WHTTFFTT: " << connected_components::total_component_count(hg, _context);
        LOG << "size_a: " << size_a;
        LOG << "size_b: " << size_b;
        LOG << "Component: " << components[0].nodes.size();
    }

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(InitialPartitioningAlgorithm::st, _rng, _tag, time);
  }
}

template<typename TypeTraits>
void GreedySTInitialPartitioner<TypeTraits>::calculate_component_spanning_tree(
    const PartitionedHypergraph& phg,
    ConnectedComponent& component,
    vec<HypernodeID>& hn_to_parent,
    vec<vec<HypernodeID>>& hn_to_children,
    vec<size_t>& subtree_size,
    vec<size_t>& hn_to_num_children,
    const HypernodeID& starter_node,
    HypernodeID& farthest_leaf_node,
    const Bitset& covered
) {
    farthest_leaf_node = kInvalidHypernode;

    for (const HypernodeID& node : phg.nodes()) {
        hn_to_children[node].clear();
    }

    Bitset node_colored;
    node_colored.resize(phg.initialNumNodes());
    node_colored.set((size_t) starter_node);

    Bitset edge_colored;
    edge_colored.resize(phg.initialNumEdges());

    std::deque<HyperedgeID> node_queue;
    node_queue.push_back(starter_node);

    std::deque<HyperedgeID> lp_node_queue;
    
    Bitset covered_nb;
    covered_nb.resize(phg.initialNumNodes());

    for (const HypernodeID& node : phg.nodes()) {
        if (covered.isSet((size_t) node)) {
            for (const HyperedgeID& he : phg.incidentEdges(node)) {

                if (edge_colored.isSet((size_t) he)) {
                    continue;
                }

                edge_colored.set((size_t) he);

                for (const HypernodeID& incident_he : phg.pins(he)) {
                    if (covered.isSet((size_t) incident_he)) {
                        continue;
                    }

                    covered_nb.set((size_t) incident_he);
                }
            }
        }
    }

    edge_colored.reset();
    std::deque<HypernodeID> calculation_queue;


    while (node_queue.size() > 0 || lp_node_queue.size() > 0) {

        HypernodeID current_node;

        if (node_queue.size() > 0) {
            current_node = node_queue.back();
            node_queue.pop_back();
        }
        else {
            current_node = lp_node_queue.back();
            lp_node_queue.pop_back();
        }

        for (const HyperedgeID& he : phg.incidentEdges(current_node)) {

            if (edge_colored.isSet((size_t) he)) {
                continue;
            }

            edge_colored.set((size_t) he);

            for (const HypernodeID& incident_hn : phg.pins(he)) {
                if (incident_hn == current_node) {
                    continue;
                }

                if (node_colored.isSet((size_t) incident_hn)) {
                    continue;
                }

                if (covered.isSet((size_t) incident_hn)) {
                    continue;
                }

                node_colored.set((size_t) incident_hn);

                hn_to_parent[incident_hn] = current_node;
                hn_to_num_children[current_node]++;

                hn_to_children[current_node].push_back(incident_hn);
                calculation_queue.push_back(incident_hn);

                if (covered_nb.isSet((size_t) incident_hn)) {
                    lp_node_queue.push_back(incident_hn);
                }
                else {
                    node_queue.push_back(incident_hn);
                }
            }
        }

        if (node_queue.size() == 0) {
            farthest_leaf_node = current_node;
        }
    }

    for (const HypernodeID& node : phg.nodes()) {
        subtree_size[node] = phg.nodeWeight(node);
    }

    for (const HypernodeID& node : calculation_queue) {
        subtree_size[hn_to_parent[node]] += subtree_size[node];
    }
}


template<typename TypeTraits>
void GreedySTInitialPartitioner<TypeTraits>::calculate_split(
    const PartitionedHypergraph& phg,
    ConnectedComponent& component,
    const size_t target,
    vec<HypernodeID>& result,
    size_t& current_split
) {
    ////
    vec<HypernodeID> hn_to_parent;
    hn_to_parent.resize(phg.initialNumNodes());

    vec<size_t> hn_to_num_children;
    hn_to_num_children.resize(phg.initialNumNodes());

    vec<vec<HypernodeID>> hn_to_children;
    hn_to_children.resize(phg.initialNumNodes());

    vec<size_t> subtree_size;
    subtree_size.resize(phg.initialNumNodes());

    Bitset covered;
    covered.resize(phg.initialNumNodes());
    ////

    HypernodeID starter_node    = kInvalidHypernode;

    size_t max_iterations       = 20;
    size_t current_iteration    = 0;

    Bitset node_colored;
    node_colored.resize(phg.initialNumNodes());

    std::deque<HyperedgeID> node_queue;

    while (current_split < target * (1.0 - _context.partition.epsilon) && current_iteration < max_iterations) {
        HypernodeID starter_node_st     = kInvalidHypernode;
        HypernodeID farthest_leaf_node  = kInvalidHypernode;

        std::shuffle(component.nodes.begin(), component.nodes.end(), _rng);

        for (const HypernodeID& node : component.nodes) {
            if (!covered.isSet((size_t) node)) {
                starter_node_st = node;
                hn_to_num_children[node] = 0;
            }
        }

        calculate_component_spanning_tree(phg, component, hn_to_parent, hn_to_children, subtree_size, hn_to_num_children, starter_node_st, farthest_leaf_node, covered);

        if (current_split == 0) {

            for (const HypernodeID& node : component.nodes) {
                if (hn_to_num_children[node] == 0 && phg.nodeWeight(node) < target) {
                    starter_node = node;
                }
            }

            node_colored.set((size_t) starter_node);
            result.push_back(starter_node);
            node_queue.push_back(starter_node);
            current_split += phg.nodeWeight(starter_node);
        }

        while (node_queue.size() > 0) {

            HypernodeID current_node = node_queue.front();
            node_queue.pop_front();

            for (const HyperedgeID& he : phg.incidentEdges(current_node)) {

                for (const HypernodeID& incident_hn : phg.pins(he)) {
                    if (incident_hn == current_node) {
                        continue;
                    }

                    if (node_colored.isSet((size_t) incident_hn)) {
                        continue;
                    }

                    if (hn_to_num_children[incident_hn] > 0 && subtree_size[incident_hn] + current_split > target) {
                        continue;
                    }

                    if (current_split + phg.nodeWeight(incident_hn) > target) {
                        continue;
                    }

                    if (subtree_size[incident_hn] + current_split <= target) {

                        std::deque<HypernodeID> asignment_queue;
                        asignment_queue.push_back(incident_hn);
                        
                        while (asignment_queue.size() > 0) {
                            HypernodeID current_assignment_node = asignment_queue.back();
                            asignment_queue.pop_back();

                            HypernodeID parent_of_incident_hn = hn_to_parent[current_assignment_node];
                            hn_to_num_children[current_assignment_node]--;
                            
                            node_colored.set((size_t) current_assignment_node);
                            node_queue.push_back(current_assignment_node);

                            current_split += phg.nodeWeight(current_assignment_node);

                            result.push_back(current_assignment_node);
                            covered.set((size_t) current_assignment_node);


                            for (const HypernodeID& ch_hn : hn_to_children[current_assignment_node]) {
                                if (covered.isSet((size_t) ch_hn)) {
                                    continue;
                                }

                                asignment_queue.push_back(ch_hn);
                            }
                        }
                    }
                    else {
                        HypernodeID parent_of_incident_hn = hn_to_parent[incident_hn];
                        hn_to_num_children[parent_of_incident_hn]--;
                        
                        node_colored.set((size_t) incident_hn);
                        node_queue.push_back(incident_hn);

                        current_split += phg.nodeWeight(incident_hn);

                        result.push_back(incident_hn);
                        covered.set((size_t) incident_hn);
                    }
                }
            }   
        }

        for (const HypernodeID& node : result) {
            node_queue.push_back(node);
        }

        current_iteration++;
    }

}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(GreedySTInitialPartitioner)

} // namespace mt_kahypar
