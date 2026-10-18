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

#pragma once

#include <algorithm>
#include <deque>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/connected_components/tarjan.h"
#include "mt-kahypar/utils/connected_component_stats.h"

namespace mt_kahypar {
namespace connected_components {

// Shared region selection and branch compaction for initial partitioning and refinement.
template<typename PartitionedHypergraph>
class TarjanRegionCoarsening {
  using PackedComponentID = uint32_t;
  using Bitset = ds::Bitset;
  struct RegionInfo {
    size_t total_weight;
    utils::NodeType type;
    vec<HypernodeID> nodes;
    std::set<HypernodeID> connected_to;
  };

 public:
  TarjanRegionCoarsening(const Context& context, std::mt19937& rng) :
    _context(context), _rng(rng) { }

    parallel::scalable_vector<HypernodeID> computeCommunities(const PartitionedHypergraph& hg) {
      connected_components::Tarjan<PartitionedHypergraph> tarjan;
      vec<PackedComponentID> vertex_to_region;
      vec<RegionInfo> region_info;
      compact_regions(hg, vertex_to_region, region_info, tarjan);

      std::shuffle(region_info.begin(), region_info.end(), _rng);
      for (size_t i = 0; i < region_info.size(); ++i) {
        for (const HypernodeID node : region_info[i].nodes) {
          vertex_to_region[node] = i;
        }
        std::shuffle(region_info[i].nodes.begin(), region_info[i].nodes.end(), _rng);
      }

      vec<size_t> subtree_size;
      vec<PackedComponentID> region_to_parent;
      vec<PackedComponentID> heads;
      calculate_master_spanning_tree(vertex_to_region, region_info, subtree_size, region_to_parent, heads);
      parallel::scalable_vector<HypernodeID> communities(hg.initialNumNodes());
      calculate_communities(region_info, region_to_parent, communities, subtree_size, heads);
      return communities;
    }

    // Keep the existing partition representable and each coarse vertex connected.
    // Fixed vertices remain singletons, so their constraints cannot conflict.
    static void respectPartition(const PartitionedHypergraph& hg,
                                 parallel::scalable_vector<HypernodeID>& communities) {
      vec<HypernodeID> groups(hg.initialNumNodes(), kInvalidHypernode);
      std::map<std::pair<HypernodeID, PartitionID>, HypernodeID> group_ids;
      HypernodeID next_group = 0;
      for (const HypernodeID node : hg.nodes()) {
        if (hg.isFixed(node)) {
          groups[node] = next_group++;
        } else {
          const auto key = std::make_pair(communities[node], hg.partID(node));
          auto entry = group_ids.emplace(key, next_group);
          if (entry.second) ++next_group;
          groups[node] = entry.first->second;
        }
      }

      vec<HypernodeID> edge_group(hg.initialNumEdges(), kInvalidHypernode);
      ds::Bitset visited(hg.initialNumNodes());
      std::queue<HypernodeID> queue;
      HypernodeID next_community = 0;
      for (const HypernodeID root : hg.nodes()) {
        if (visited.isSet(root)) continue;
        const HypernodeID group = groups[root];
        visited.set(root);
        queue.push(root);
        while (!queue.empty()) {
          const HypernodeID node = queue.front();
          queue.pop();
          communities[node] = next_community;
          for (const HyperedgeID edge : hg.incidentEdges(node)) {
            if (edge_group[edge] == group) continue;
            edge_group[edge] = group;
            for (const HypernodeID pin : hg.pins(edge)) {
              if (!visited.isSet(pin) && groups[pin] == group) {
                visited.set(pin);
                queue.push(pin);
              }
            }
          }
        }
        next_community++;
      }
    }

 private:
void calculate_communities(
  const vec<RegionInfo>& region_info,
  vec<PackedComponentID>& region_to_parent,
  parallel::scalable_vector<HypernodeID>& communities,
  const vec<size_t>& subtree_size,
  const vec<PackedComponentID> heads

) {
  //// get packed region to head
  vec<PackedComponentID> region_to_head;
  region_to_head.resize(region_info.size());

  std::iota(region_to_head.begin(), region_to_head.end(), PackedComponentID{0});

  Bitset updated;
  updated.resize(region_info.size());

  vec<PackedComponentID> to_update;

  for (PackedComponentID region_id = 0; region_id < region_info.size(); ++region_id) {
    if (region_to_parent[region_id] == region_id) {
      region_to_head[region_id] = region_id;
      continue;
    }

    PackedComponentID parent = region_id;
    while (parent != region_to_parent[parent]) {
      if (updated.isSet((size_t) parent)) {
        parent = region_to_head[parent];
        break;
      }
      to_update.push_back(parent);
      parent = region_to_parent[parent];
    }

    for (const PackedComponentID& comp : to_update) {
      updated.set((size_t) comp);
      region_to_head[comp] = parent;
    }

    to_update.clear();
  }

  //// find subtree size for each head with the correct size in the same tree
  vec<size_t>             best_size(region_info.size(), std::numeric_limits<size_t>::max());
  vec<PackedComponentID>  best_region(region_info.size(), std::numeric_limits<uint32_t>::max());

  vec<size_t>             best_size_articulation(region_info.size(), std::numeric_limits<size_t>::max());
  vec<PackedComponentID>  best_region_articulation(region_info.size(), std::numeric_limits<uint32_t>::max());


  if (_context.extra_options.tarjan_min_max_region) {
    
    vec<size_t> largest_child_subtree(region_info.size(), 0);
    for (PackedComponentID region_id = 0; region_id < region_info.size(); ++region_id) {
      const PackedComponentID parent = region_to_parent[region_id];
      if (parent != region_id) {
        largest_child_subtree[parent] =
          std::max(largest_child_subtree[parent], subtree_size[region_id]);
      }
    }

    for (PackedComponentID region_id = 0; region_id < region_info.size(); ++region_id) {
      const PackedComponentID head = region_to_head[region_id];
      const size_t tree_weight = subtree_size[head];
      
      if (subtree_size[region_id] <= tree_weight / 2) {
        continue;
      }
      const size_t largest_branch = std::max(tree_weight - subtree_size[region_id], largest_child_subtree[region_id]);
      
      if (best_region[head] == std::numeric_limits<uint32_t>::max() ||
          (region_info[region_id].type == utils::NodeType::normal &&
           region_info[best_region[head]].type != utils::NodeType::normal) ||
          (region_info[region_id].type == region_info[best_region[head]].type &&
           (largest_branch < best_size[head] ||
            (largest_branch == best_size[head] && region_id < best_region[head])))) {
        best_size[head] = largest_branch;
        best_region[head] = region_id;
      }
    }
  } else {
    for (PackedComponentID region_id = 0; region_id < region_info.size(); ++region_id) {
      PackedComponentID head = region_to_head[region_id];
      size_t target = subtree_size[head] / 2;

      if (subtree_size[region_id] >= target && subtree_size[region_id] < best_size[head] && region_info[region_id].type == mt_kahypar::utils::NodeType::normal) {
        best_size[head]       = subtree_size[region_id];
        best_region[head]     = region_id;
      }
      else if (subtree_size[region_id] >= target && subtree_size[region_id] < best_size_articulation[head] && region_info[region_id].type == mt_kahypar::utils::NodeType::articulation) {
        best_size_articulation[head]    = subtree_size[region_id];
        best_region_articulation[head]  = region_id;
      }
    }

    // take region with articulation points, if there is no normal region with a reasonable subtree size
    for (PackedComponentID region_id = 0; region_id < region_info.size(); ++region_id) {
      PackedComponentID head = region_to_head[region_id];

      if (head != region_id) {
        continue;
      }

      size_t target         = subtree_size[head] / 2;
      size_t upper_target   = target * (1 + _context.partition.epsilon);

      if (best_size[head] > best_size_articulation[head] && best_size[head] > upper_target) {
        best_size[head]   = best_size_articulation[head];
        best_region[head] =  best_region_articulation[head];
      }
    }
  }

  ////

  //mt_kahypar::utils::cc_debug.s_initialize_regions_tarjan(pci_reduced, );

  //// find communities
  updated.reset();
  
  HypernodeID         community = 0;
  
  vec<HypernodeID>    region_to_community;
  region_to_community.resize(region_info.size());

  for (PackedComponentID region_id = 0; region_id < region_info.size(); ++region_id) {
    PackedComponentID head = region_to_head[region_id];

    if (region_id == best_region[head]) {
      continue;
    }

    PackedComponentID parent = region_id;
    while (parent != region_to_parent[parent]) {
      if (updated.isSet((size_t) parent) || parent == best_region[head]) {
        break;
      }

      to_update.push_back(parent);
      parent = region_to_parent[parent];
    }

    HypernodeID parent_community = kInvalidHypernode;

    if (parent == best_region[head]) {
      parent_community = community++;
    } else if (parent == region_to_parent[parent] && !updated.isSet((size_t) parent)) {
      parent_community = community++;
      to_update.push_back(parent);
    } 
    else {
      parent_community = region_to_community[parent];
    }

    for (const PackedComponentID& comp : to_update) {
      updated.set((size_t) comp);
      region_to_community[comp] = parent_community;
    }

    to_update.clear();
  }
  ////

  //// assign nodes
  for (PackedComponentID region_id = 0; region_id < region_info.size(); ++region_id) {
    PackedComponentID head = region_to_head[region_id];

    if (region_id == best_region[head]) {
      for (const HypernodeID& node : region_info[region_id].nodes) {
        communities[node] = community++;
      }
      continue;
    }

    for (const HypernodeID& node : region_info[region_id].nodes) {
      communities[node] = region_to_community[region_id];
    }
  }
  ////
}

void compact_regions(
  const PartitionedHypergraph& hypergraph,
  vec<PackedComponentID>& vertex_to_region,
  vec<RegionInfo>& region_info,
  connected_components::Tarjan<PartitionedHypergraph>& tarjan
) {

  vertex_to_region.resize(hypergraph.initialNumNodes());

  Bitset node_colored;
  node_colored.resize(hypergraph.initialNumNodes());
  
  Bitset edge_colored;
  edge_colored.resize(hypergraph.initialNumEdges());

  std::queue<HypernodeID> node_queue;

  for (const HypernodeID& hn : hypergraph.nodes()) {

    if (node_colored.isSet((size_t) hn)) {
      continue;
    }

    utils::NodeType current_node_type = utils::NodeType::normal;
    
    if (tarjan.is_articulation_point(hypergraph, hn)) {
      current_node_type = utils::NodeType::articulation;
    }

    const PackedComponentID current_region = region_info.size();
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
      vertex_to_region[current] = current_region;

      for (const HyperedgeID& he : hypergraph.incidentEdges(current)) {

        if (edge_colored.isSet((size_t) he)) {
          continue;
        }

        edge_colored.set((size_t) he);

        for (const HypernodeID& incident_hn : hypergraph.pins(he)) {
          if (current_node_type == utils::NodeType::normal && tarjan.is_articulation_point(hypergraph, incident_hn)) {
            border_nodes.insert(incident_hn);
            continue;
          }

          if (current_node_type == utils::NodeType::articulation && !tarjan.is_articulation_point(hypergraph, incident_hn)) {
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
    
    region_info.push_back({total_weight, current_node_type, nodes, border_nodes});
  }
};

void calculate_master_spanning_tree(
  const vec<PackedComponentID>& vertex_to_region,
  const vec<RegionInfo>& region_info,
  vec<size_t>& subtree_size,
  vec<PackedComponentID>& region_to_parent,
  vec<PackedComponentID>& heads
) {

  //// setup variables
  region_to_parent.resize(region_info.size());

  Bitset region_used;
  region_used.resize(region_info.size());

  Bitset region_colored;
  region_colored.resize(region_info.size());

  vec<PackedComponentID> region_queue;
  std::deque<PackedComponentID> calculation_queue;
  ////

  std::iota(region_to_parent.begin(), region_to_parent.end(), PackedComponentID{0});

  vec<PackedComponentID> region_order(region_info.size());
  std::iota(region_order.begin(), region_order.end(), PackedComponentID{0});

  std::sort(region_order.begin(), region_order.end(),
    [&](const PackedComponentID a, const PackedComponentID b) {
      return region_info[a].nodes.size() < region_info[b].nodes.size();
    });

  for (const PackedComponentID parent : region_order) {
    if (region_colored.isSet((size_t) parent)) {
      continue;
    }

    heads.push_back(parent);

    region_queue.push_back(parent);
    region_colored.set((size_t) parent);

    while (region_queue.size() > 0) {
      PackedComponentID current_id = region_queue.back();
      region_queue.pop_back();

      const RegionInfo& current_region = region_info[current_id];

      for (const HypernodeID& incident_hn : current_region.connected_to) {
        PackedComponentID region_id = vertex_to_region[incident_hn];

        if (region_colored.isSet((size_t) region_id)) {
          continue;
        }

        region_queue.push_back(region_id);
        calculation_queue.push_front(region_id);

        region_colored.set((size_t) region_id);
        region_to_parent[region_id] = current_id; 
      } 
    }    
  }

  // calculate spanning tree sizes
  subtree_size.resize(region_info.size());

  for (PackedComponentID region_id = 0; region_id < region_info.size(); ++region_id) {
    subtree_size[region_id] = region_info[region_id].total_weight;
  }

  for (const PackedComponentID& region_id : calculation_queue) {
    subtree_size[region_to_parent[region_id]] += subtree_size[region_id];
  }
};


  const Context& _context;
  std::mt19937& _rng;
};

} // namespace connected_components
} // namespace mt_kahypar
