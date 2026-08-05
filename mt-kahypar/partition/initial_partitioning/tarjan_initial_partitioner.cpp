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

#include "mt-kahypar/partition/initial_partitioning/tarjan_initial_partitioner.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

template<typename TypeTraits>
void TarjanInitialPartitioner<TypeTraits>::partitionImpl() {
  if ( _ip_data.should_initial_partitioner_run(InitialPartitioningAlgorithm::tarjan) ) {
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();
    std::uniform_int_distribution<PartitionID> select_random_block(0, _context.partition.k - 1);
    
    // compute components and do tarjan for each component
    connected_components::Tarjan<PartitionedHypergraph> tarjan;
    tarjan.initialize(hg);

    vec<PackedComponentID> vertex_to_packed_component;
    vec<PackedComponentInfo> packed_component_info;
    compact_regions(hg, vertex_to_packed_component, packed_component_info, tarjan);

    for (const PackedComponentInfo& comp_info : packed_component_info) {

      if (comp_info.nodes.size() <= 2) {
        continue;
      }

      LOG << "Type:         " << comp_info.type;
      //LOG << "total_weight: " << comp_info.total_weight;
      LOG << "Nodes count:  " << comp_info.nodes.size();
      //LOG << "Border nodes: " << comp_info.connected_to.size();

      LOG << "";

      /*if (comp_info.type == NodeType::normal) {
        for (const HypernodeID& node : comp_info.connected_to) {
          LOG << "Border node: " << node;
        }
      } else if (comp_info.type == NodeType::articulation) {
        for (const HypernodeID& node : comp_info.nodes) {
          LOG << "Internal node: " << node;
        }
      }*/

      //LOG << "#################";
    }

    while (true);
    

    vec<PackedComponentID>  component_to_parent;
    vec<size_t>             subtree_size;

    calculate_master_spanning_tree(vertex_to_packed_component, packed_component_info, component_to_parent, subtree_size);
    /*for (const PackedComponentInfo& comp_info : packed_component_info) {
      LOG << "Component id: " << comp_info.id;
      LOG << "Own size:     " << comp_info.total_weight;
      LOG << "Subtree size: " << subtree_size[comp_info.id];
      LOG << "Parent:       " << component_to_parent[comp_info.id];
    }*/

    // for each component compact regions

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(InitialPartitioningAlgorithm::tarjan, _rng, _tag, time);
  }
};

template<typename TypeTraits>
void TarjanInitialPartitioner<TypeTraits>::compact_regions(
  const PartitionedHypergraph& hypergraph,
  vec<PackedComponentID>& vertex_to_packed_component,
  vec<PackedComponentInfo>& packed_component_info,
  connected_components::Tarjan<PartitionedHypergraph>& tarjan
) {

  vertex_to_packed_component.resize(hypergraph.initialNumNodes());

  Bitset node_colored;
  node_colored.resize(hypergraph.initialNumNodes());
  
  Bitset edge_colored;
  edge_colored.resize(hypergraph.initialNumEdges());

  std::queue<HypernodeID> node_queue;
  PackedComponentID current_component = 0;

  for (const HypernodeID& hn : hypergraph.nodes()) {

    if (node_colored.isSet((size_t) hn)) {
      continue;
    }

    NodeType current_node_type = NodeType::normal;
    
    if (tarjan.is_articulation_point(hypergraph, hn)) {
      current_node_type = NodeType::articulation;
    }

    node_queue.push(hn);
    node_colored.set((size_t) hn);

    vec<HypernodeID> nodes;
    std::set<HypernodeID> border_nodes;
    size_t total_weight = 0;
    
    edge_colored.reset();

    while (node_queue.size() > 0) {
      HypernodeID current = node_queue.front();
      node_queue.pop();

      nodes.push_back(current);
      total_weight += hypergraph.nodeWeight(current);
      vertex_to_packed_component[current] = current_component;

      for (const HyperedgeID& he : hypergraph.incidentEdges(current)) {

        if (edge_colored.isSet((size_t) he)) {
          continue;
        }

        edge_colored.set((size_t) he);

        for (const HypernodeID& incident_hn : hypergraph.pins(he)) {
          if (node_colored.isSet((size_t) incident_hn)) {
            continue;
          }

          if (current_node_type == NodeType::normal && tarjan.is_articulation_point(hypergraph, incident_hn)) {
            border_nodes.insert(incident_hn);
            continue;
          }

          if (current_node_type == NodeType::articulation && !tarjan.is_articulation_point(hypergraph, incident_hn)) {
            border_nodes.insert(incident_hn);
            continue;
          }

          node_colored.set((size_t) incident_hn);
          node_queue.push(incident_hn);
        }
      }
    }
    
    packed_component_info.push_back({(uint32_t) current_component, total_weight, current_node_type, nodes, border_nodes});
    current_component++;
  }

};

template<typename TypeTraits>
void TarjanInitialPartitioner<TypeTraits>::calculate_master_spanning_tree(
  const vec<PackedComponentID>& vertex_to_packed_component,
  const vec<PackedComponentInfo>& packed_component_info,
  vec<PackedComponentID>& component_to_parent,
  vec<size_t>& subtree_size
) {

  // find the biggest component
  PackedComponentID biggest_comp_id;
  size_t biggest = 0;
  for (const PackedComponentInfo& comp_info : packed_component_info) {
    if (comp_info.total_weight > biggest) {
      biggest           = comp_info.total_weight;
      biggest_comp_id   = comp_info.id; 
    }
  }

  // make biggest component root of the spanning tree
  component_to_parent.resize(packed_component_info.size());

  Bitset component_colored;
  component_colored.resize(packed_component_info.size());

  vec<PackedComponentID> component_queue;
  vec<PackedComponentID> calculation_queue;

  component_queue.push_back(biggest_comp_id);
  component_colored.set((size_t) biggest_comp_id);

  for (const PackedComponentInfo& comp_info : packed_component_info) {
    if (component_queue.size() == 0) {
      if (component_colored.isSet((size_t) comp_info.id)) {
        continue;
      }

      component_queue.push_back(comp_info.id);
      component_colored.set((size_t) comp_info.id);
    }    
  
    while (component_queue.size() > 0) {
      PackedComponentID current_id = component_queue.back();
      component_queue.pop_back();

      PackedComponentInfo current_component = packed_component_info[current_id];

      for (const HypernodeID& incident_hn : current_component.connected_to) {
        PackedComponentID component_id = vertex_to_packed_component[incident_hn];

        if (component_colored.isSet((size_t) component_id)) {
          continue;
        }

        component_queue.push_back(component_id);
        calculation_queue.push_back(component_id);

        component_colored.set((size_t) component_id);
        component_to_parent[component_id] = current_id; 
      } 
    }
  }

  // calculate spanning tree sizes
  subtree_size.resize(packed_component_info.size());

  for (const PackedComponentInfo& component : packed_component_info) {
    subtree_size[component.id] = component.total_weight;
  }

  for (const PackedComponentID& component_id : calculation_queue) {
    subtree_size[component_to_parent[component_id]] += subtree_size[component_id];
  }
};


INSTANTIATE_CLASS_WITH_TYPE_TRAITS(TarjanInitialPartitioner)

} // namespace mt_kahypar
