
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

const size_t MAX_ORIGINS  = 1;
const size_t MAX_SPLITS   = 3;

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

    vec<HypernodeID> nodes_to_swap;
    nodes_to_swap.reserve(hg.initialNumNodes());

    vec<HypernodeID> best_split;
    best_split.reserve(hg.initialNumNodes());

    size_t current_origins = 0;

    ////// Split the components, only if there is not enough size
    for (const std::pair<size_t, connected_components::ConnectedComponent>& component_and_size : components_and_size) {

        best_split.clear();
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

              calculate_split(hg, component, target_for_split * (1.0 + _context.partition.epsilon), nodes_to_swap, split_size);

              diff = target_for_split >= split_size ? target_for_split - split_size : split_size - target_for_split;
              current_origins++;

              if (static_cast<double>(diff) / target_for_split < best_split_diff) {
                best_split = nodes_to_swap;
                best_split_diff = static_cast<double>(diff) / target_for_split;
              }

          } while(current_origins <= MAX_ORIGINS);

          //// assign nodes from best split            
          for (const HypernodeID& node : best_split) {
            hg.changeNodePart(node, 0, 1, DynamicConnectivityStrategy::do_nothing);
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

    vec<vec<connected_components::ConnectedComponent>> extra_components;
    connected_components::compute_components_per_block(hg, _context, extra_components);
    for (const vec<connected_components::ConnectedComponent>& components_per_partition : extra_components) {
        if (components_per_partition.size() != 1) {
            LOG << "components in partititon=" << components_per_partition.size();
        } 
    }

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(InitialPartitioningAlgorithm::st, _rng, _tag, time);
  }
}

template<typename TypeTraits>
inline void GreedySTInitialPartitioner<TypeTraits>::calculate_component_spanning_tree(
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
        hn_to_parent[node] = node;
        hn_to_num_children[node] = 0;
    }

    Bitset node_colored;
    node_colored.resize(phg.initialNumNodes());
    node_colored.set((size_t) starter_node);

    Bitset edge_colored;
    edge_colored.resize(phg.initialNumEdges());

    vec<HyperedgeID> node_queue;
    node_queue.reserve(phg.initialNumNodes());
    node_queue.push_back(starter_node);

    vec<HyperedgeID> lp_node_queue;
    lp_node_queue.reserve(phg.initialNumNodes());
    
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
                calculation_queue.push_front(incident_hn);

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

    size_t current_iteration    = 0;

    Bitset node_colored;
    node_colored.resize(phg.initialNumNodes());

    Bitset edge_colored;
    edge_colored.resize(phg.initialNumEdges());

    Bitset already_started_from;
    already_started_from.resize(phg.initialNumNodes());

    std::deque<HyperedgeID> node_queue;

    vec<HypernodeID> asignment_queue;
    asignment_queue.reserve(phg.initialNumNodes());

    while (current_split < target && current_iteration <= MAX_SPLITS) {
        current_iteration++;

        HypernodeID starter_node_st     = kInvalidHypernode;
        HypernodeID farthest_leaf_node  = kInvalidHypernode;

        if (!component.nodes.empty()) {
          const size_t start = std::uniform_int_distribution<size_t>(
              0, component.nodes.size() - 1
          )(_rng);

          for (size_t i = 0; i < component.nodes.size(); ++i) {
              const HypernodeID node = component.nodes[(start + i) % component.nodes.size()];

              if (!covered.isSet(static_cast<size_t>(node))) {

                if (starter_node_st == kInvalidHypernode) {
                    starter_node_st = node;
                }
                else if (already_started_from.isSet((size_t) starter_node_st) && !already_started_from.isSet((size_t) node)) {
                    starter_node_st = node;
                }

                hn_to_num_children[node] = 0;
                break;
              }
          }
        }

        already_started_from.set((size_t) starter_node_st);

        //LOG << "starter_node_st: " << starter_node_st;

        calculate_component_spanning_tree(phg, component, hn_to_parent, hn_to_children, subtree_size, hn_to_num_children, starter_node_st, farthest_leaf_node, covered);

        if (result.size() == 0) {

            for (const HypernodeID& node : component.nodes) {
                if (hn_to_num_children[node] == 0 && phg.nodeWeight(node) < target) {
                    starter_node = node;
                }
            }

            if (starter_node == kInvalidHypernode) {
                continue;
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

                    if (current_split + phg.nodeWeight(incident_hn) > target) {
                        continue;
                    }

                    uint32_t true_size = phg.nodeWeight(incident_hn);

                    if (subtree_size[incident_hn] + current_split <= target) {

                        add_node_to_split(
                            incident_hn,
                            hn_to_parent,
                            hn_to_num_children,
                            node_colored,
                            covered,
                            current_split,
                            phg.nodeWeight(incident_hn),
                            node_queue,
                            result
                        );

                        asignment_queue.clear();
                        asignment_queue.push_back(incident_hn);
                        
                        while (asignment_queue.size() > 0) {
                            HypernodeID current_assignment_node = asignment_queue.back();
                            asignment_queue.pop_back();

                            for (const HypernodeID& ch_hn : hn_to_children[current_assignment_node]) {
                                if (covered.isSet((size_t) ch_hn) || node_colored.isSet((size_t) ch_hn)) {
                                    continue;
                                }

                                true_size += phg.nodeWeight(ch_hn);

                                add_node_to_split(
                                    ch_hn,
                                    hn_to_parent,
                                    hn_to_num_children,
                                    node_colored,
                                    covered,
                                    current_split,
                                    phg.nodeWeight(ch_hn),
                                    node_queue,
                                    result
                                );

                                asignment_queue.push_back(ch_hn);
                            }
                        }

                        if (true_size > subtree_size[incident_hn]) {
                            //LOG << "true size: " << true_size << " subtree_size: " << subtree_size[incident_hn];
                        }
                        
                    }
                    else if (hn_to_num_children[incident_hn] == 0) {
                        add_node_to_split(
                            incident_hn,
                            hn_to_parent,
                            hn_to_num_children,
                            node_colored,
                            covered,
                            current_split,
                            phg.nodeWeight(incident_hn),
                            node_queue,
                            result
                        );
                    } else {
                        continue;
                    }

                    HypernodeID parent = hn_to_parent[incident_hn];
                    while (parent != hn_to_parent[parent]) {
                        subtree_size[parent] -= true_size;
                        parent = hn_to_parent[parent];
                    }

                    subtree_size[parent] -= true_size;

                }
            }   
        }

        //LOG << "Iteration: " << current_iteration << " and gain: " << current_split << " with target: " << target;

        for (const HypernodeID& node : result) {
            node_queue.push_back(node);
        }

        edge_colored.reset();
    }

}

template<typename TypeTraits>
void GreedySTInitialPartitioner<TypeTraits>::add_node_to_split(
    const HypernodeID& hn,
    const vec<HypernodeID>& hn_to_parent,
    vec<size_t>& hn_to_num_children,
    Bitset& node_colored,
    Bitset& node_covered,
    size_t& current_split,
    const size_t& node_weight,
    std::deque<HypernodeID>& node_queue,
    vec<HypernodeID>& result
) {
    HypernodeID parent_of_hn = hn_to_parent[hn];
    hn_to_num_children[parent_of_hn]--;
    
    node_colored.set((size_t) hn);
    node_queue.push_back(hn);

    current_split += node_weight;

    result.push_back(hn);
    node_covered.set((size_t) hn);
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(GreedySTInitialPartitioner)

} // namespace mt_kahypar
