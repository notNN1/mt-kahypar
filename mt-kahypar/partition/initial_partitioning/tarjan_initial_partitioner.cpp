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
#include "mt-kahypar/partition/initial_partitioning/initial_partitioning_data_container.h"

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

    //// shuffle everything
    std::shuffle(packed_component_info.begin(), packed_component_info.end(), _rng);
    for (size_t i = 0; i < packed_component_info.size(); i++) {
      packed_component_info[i].id = i;

      for (const HypernodeID& node : packed_component_info[i].nodes) {
        vertex_to_packed_component[node] = i;
      }
    }

    for (PackedComponentInfo& pci : packed_component_info) {
      std::shuffle(pci.nodes.begin(), pci.nodes.end(), _rng);
    }    
    ////

    //// calculate spanning tree
    vec<size_t>                         subtree_size;
    vec<PackedComponentID>              component_to_parent;
    vec<PackedComponentID>              heads;

    calculate_master_spanning_tree(vertex_to_packed_component, packed_component_info, subtree_size, component_to_parent, heads);

    
    //// calculate communities and create new hypergraph
    parallel::scalable_vector<HypernodeID> communities(hg.initialNumNodes());
    calculate_communities(packed_component_info, component_to_parent, communities, subtree_size, heads);
    
    auto& old_hg    = hg.hypergraph();
    auto new_hg     = old_hg.contract(communities, true);

    PartitionedHypergraph new_phg(hg.k(), new_hg);
    ////


    //// do initial partitioning on the new phg for each enabled initial partitioner
    InitialPartitioningDataContainer<TypeTraits> ip_data(new_phg, _context);

    auto* ip_data_ptr = ip::to_pointer(ip_data);

    for ( uint8_t i = 0; i < static_cast<uint8_t>(InitialPartitioningAlgorithm::UNDEFINED); ++i ) {
      if ( static_cast<InitialPartitioningAlgorithm>(i) != InitialPartitioningAlgorithm::tarjan 
      && _context.initial_partitioning.enabled_ip_algos[i]) {
        auto algorithm = static_cast<InitialPartitioningAlgorithm>(i);

        // Create one initial partitioner.
        std::unique_ptr<IInitialPartitioner> initial_partitioner =
            InitialPartitionerFactory::getInstance().createObject(
                algorithm,
                algorithm,
                ip_data_ptr,
                _context,
                _context.partition.seed + i + 100,
                10000);

        // Run it.
        initial_partitioner->partition();
      }
    }
    ////


    //// Do partitioning like the best performing partitioner on the uncontracted graph
    auto& ip_data_ref = *reinterpret_cast<InitialPartitioningDataContainer<TypeTraits>*>(ip_data_ptr);
    ip_data_ref.apply();
    auto& local_phg = ip_data_ref.local_partitioned_hypergraph();

    for (const HypernodeID& node : hg.nodes()) {
      hg.setNodePart(node, local_phg.partID(communities[node]));
    }
    ////

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(InitialPartitioningAlgorithm::tarjan, _rng, _tag, time);
  }
};

template<typename TypeTraits>
void TarjanInitialPartitioner<TypeTraits>::calculate_communities(
  const vec<PackedComponentInfo>& packed_component_info,
  vec<PackedComponentID>& component_to_parent,
  parallel::scalable_vector<HypernodeID>& communities,
  const vec<size_t>& subtree_size,
  const vec<PackedComponentID> heads

) {
  //// get packed component to head
  vec<PackedComponentID> packed_component_to_head;
  packed_component_to_head.resize(packed_component_info.size());

  for (const PackedComponentInfo& comp_info : packed_component_info) {
    packed_component_to_head[comp_info.id] = comp_info.id;
  }

  Bitset updated;
  updated.resize(packed_component_info.size());

  vec<PackedComponentID> to_update;

  for (const PackedComponentInfo& comp_info : packed_component_info) {
    if (component_to_parent[comp_info.id] == comp_info.id) {
      packed_component_to_head[comp_info.id] = comp_info.id;
      continue;
    }

    PackedComponentID parent = comp_info.id;
    while (parent != component_to_parent[parent]) {
      if (updated.isSet((size_t) parent)) {
        parent = packed_component_to_head[parent];
        break;
      }
      to_update.push_back(parent);
      parent = component_to_parent[parent];
    }

    for (const PackedComponentID& comp : to_update) {
      updated.set((size_t) comp);
      packed_component_to_head[comp] = parent;
    }

    to_update.clear();
  }

  for (const PackedComponentInfo& comp_info : packed_component_info) {

  }

  //// find subtree size for each head with the correct size in the same tree
  vec<size_t>             best_size(packed_component_info.size(), std::numeric_limits<size_t>::max());
  vec<PackedComponentID>  best_component(packed_component_info.size(), std::numeric_limits<uint32_t>::max());

  for (const PackedComponentInfo& comp_info : packed_component_info) {
    PackedComponentID head = packed_component_to_head[comp_info.id];
    size_t target = subtree_size[head] / 2;

    if (subtree_size[comp_info.id] > target && subtree_size[comp_info.id] < best_size[head]) {
      best_size[head]       = subtree_size[comp_info.id];
      best_component[head]  = comp_info.id;
    }
  }
  ////

  //// find communities
  updated.reset();
  
  HypernodeID         community = 0;
  
  vec<HypernodeID>    packed_component_to_community;
  packed_component_to_community.resize(packed_component_info.size());

  for (const PackedComponentInfo& comp_info : packed_component_info) {
    PackedComponentID head = packed_component_to_head[comp_info.id];

    if (comp_info.id == best_component[head]) {
      continue;
    }

    PackedComponentID parent = comp_info.id;
    while (parent != component_to_parent[parent]) {
      if (updated.isSet((size_t) parent) || parent == best_component[head]) {
        break;
      }

      to_update.push_back(parent);
      parent = component_to_parent[parent];
    }

    HypernodeID parent_community = kInvalidHypernode;

    if (parent == best_component[head]) {
      parent_community = community++;
    } else if (parent == component_to_parent[parent] && !updated.isSet((size_t) parent)) {
      parent_community = community++;
      to_update.push_back(parent);
    } 
    else {
      parent_community = packed_component_to_community[parent];
    }

    for (const PackedComponentID& comp : to_update) {
      updated.set((size_t) comp);
      packed_component_to_community[comp] = parent_community;
    }

    to_update.clear();
  }
  ////

  //// assign nodes
  for (const PackedComponentInfo& comp_info : packed_component_info) {
    PackedComponentID head = packed_component_to_head[comp_info.id];

    if (comp_info.id == best_component[head]) {
      for (const HypernodeID& node : comp_info.nodes) {
        communities[node] = community++;
      }
      continue;
    }

    for (const HypernodeID& node : comp_info.nodes) {
      communities[node] = packed_component_to_community[comp_info.id];
    }
  }
  ////
}

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
          if (current_node_type == NodeType::normal && tarjan.is_articulation_point(hypergraph, incident_hn)) {
            border_nodes.insert(incident_hn);
            continue;
          }

          if (current_node_type == NodeType::articulation && !tarjan.is_articulation_point(hypergraph, incident_hn)) {
            border_nodes.insert(incident_hn);
            continue;
          }

          if (node_colored.isSet((size_t) incident_hn)) {
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
  vec<size_t>& subtree_size,
  vec<PackedComponentID>& component_to_parent,
  vec<PackedComponentID>& heads
) {

  //// setup variables
  component_to_parent.resize(packed_component_info.size());

  Bitset component_used;
  component_used.resize(packed_component_info.size());

  Bitset component_colored;
  component_colored.resize(packed_component_info.size());

  vec<PackedComponentID> component_queue;
  std::deque<PackedComponentID> calculation_queue;
  ////

  for (const PackedComponentInfo& comp_info : packed_component_info) {
    component_to_parent[comp_info.id] = comp_info.id;
  }

  vec<PackedComponentInfo> packed_component_info_copy = packed_component_info;

  std::sort(packed_component_info_copy.begin(), packed_component_info_copy.end(),
    [](const PackedComponentInfo& a, const PackedComponentInfo& b) {
      return a.nodes.size() < b.nodes.size();
    });

  for (const PackedComponentInfo& parent : packed_component_info_copy) {
    if (component_colored.isSet((size_t) parent.id)) {
      continue;
    }

    heads.push_back(parent.id);

    component_queue.push_back(parent.id);
    component_colored.set((size_t) parent.id);

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
        calculation_queue.push_front(component_id);

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
