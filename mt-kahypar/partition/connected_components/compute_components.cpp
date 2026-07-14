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

#include "mt-kahypar/partition/connected_components/compute_components.h"

#include "mt-kahypar/definitions.h"
#include <tbb/task_group.h>
#include <set>

namespace mt_kahypar {
namespace connected_components {

using Bitset = mt_kahypar::ds::Bitset;

template<typename PartitionedHypergraph>
void compute_components_per_block(const PartitionedHypergraph& phg,
                                  const Context& context,
                                  vec<vec<ConnectedComponent>>& result) {
  result.resize(phg.k());
  for (vec<ConnectedComponent>& list: result) {
    list.clear();
  }

  (void)context;

  Bitset node_colored;
  node_colored.resize(phg.initialNumNodes());
  
  Bitset edge_colored;
  edge_colored.resize(phg.initialNumEdges());

  std::queue<HypernodeID> node_queue;
  PartitionID current_partition;

  for (const HypernodeID& hn : phg.nodes()) {
    if (node_colored.isSet((size_t) hn)) {
      continue;
    }

    current_partition = phg.partID(hn);

    if (current_partition == -1) {
      continue;
    }

    node_queue.push(hn);
    node_colored.set((size_t) hn);
    
    ConnectedComponent cc = { };
    edge_colored.reset();

    while (node_queue.size() > 0) {
      HypernodeID current = node_queue.front();
      cc.nodes.push_back(current);
      node_queue.pop();

      for (const HyperedgeID& he : phg.incidentEdges(current)) {

        if (edge_colored.isSet((size_t) he)) {
          continue;
        }

        edge_colored.set((size_t) he);

        for (const HypernodeID& incident_hn : phg.pins(he)) {
          if (node_colored.isSet((size_t) incident_hn)) {
            continue;
          }
          
          if (phg.partID(incident_hn) != current_partition) {
            continue;
          }

          node_colored.set((size_t) incident_hn);
          node_queue.push(incident_hn);
        }
      }
    }
    
    result[current_partition].push_back(cc);
  }
}

template<typename PartitionedHypergraph>
void compute_components(
  const PartitionedHypergraph& phg,
  const Context& context,
  vec<ConnectedComponent>& result
) {
  (void)context;

  Bitset node_colored;
  node_colored.resize(phg.initialNumNodes());
  
  Bitset edge_colored;
  edge_colored.resize(phg.initialNumEdges());

  std::queue<HypernodeID> node_queue;
  PartitionID current_partition;

  for (const HypernodeID& hn : phg.nodes()) {
    if (node_colored.isSet((size_t) hn)) {
      continue;
    }

    node_queue.push(hn);
    node_colored.set((size_t) hn);
    
    ConnectedComponent cc = { };
    edge_colored.reset();

    while (node_queue.size() > 0) {
      HypernodeID current = node_queue.front();
      cc.nodes.push_back(current);
      node_queue.pop();

      for (const HyperedgeID& he : phg.incidentEdges(current)) {

        if (edge_colored.isSet((size_t) he)) {
          continue;
        }

        edge_colored.set((size_t) he);

        for (const HypernodeID& incident_hn : phg.pins(he)) {
          if (node_colored.isSet((size_t) incident_hn)) {
            continue;
          }

          node_colored.set((size_t) incident_hn);
          node_queue.push(incident_hn);
        }
      }
    }

    result.push_back(cc);
  }
}

template<typename PartitionedHypergraph>
void compute_super_components(const PartitionedHypergraph& phg,
                                                              const Context& context,
                                                              const vec<vec<ConnectedComponent>>& components,
                                                              vec<vec<ComponentInfo>>& result) {

  vec<size_t> hn_to_component;
  hn_to_component.resize(phg.initialNumNodes());

  vec<size_t> component_to_size;
  vec<PartitionID> component_to_paritition;

  size_t component_i = 0;

  PartitionID partition = 0;

  for (vec<ConnectedComponent> components_per_partition : components) {
    for (ConnectedComponent component : components_per_partition) {
      for (HypernodeID hn : component.nodes) {
        hn_to_component[hn] = component_i;
      }
      component_to_size.push_back(component.nodes.size());
      component_to_paritition.push_back(partition);
      component_i++;
    }
    partition++;
  }

  vec<size_t> component_to_parent;
  component_to_parent.resize(component_i);

  for (size_t i = 0; i < component_i; i++) {
    component_to_parent[i] = i;
  }

  Bitset edge_colored;
  edge_colored.resize(phg.initialNumEdges());

  for (vec<ConnectedComponent> components_per_partition : components) {
    for (ConnectedComponent component : components_per_partition) {
      for (HypernodeID hn : component.nodes) {
        
        size_t component        = hn_to_component[hn];
        size_t parent_component = component;

        while (component_to_parent[parent_component] != parent_component) parent_component = component_to_parent[parent_component];

        for (const HyperedgeID& he : phg.incidentEdges(hn)) {

          if (edge_colored.isSet((size_t) he)) {
            continue;
          }

          edge_colored.set((size_t) he);

          for (const HypernodeID& incident_hn : phg.pins(he)) {
            size_t incident_component = hn_to_component[incident_hn];

            if (incident_component != component) {
              // calculate parent
              size_t incident_parent_component  = incident_component;
              while (component_to_parent[incident_parent_component] != incident_component) incident_component = component_to_parent[incident_component];

              if (parent_component != incident_component) {
                component_to_parent[incident_component] = parent_component;
              }
            }   
          }
        }
      }
    }
  }
  

  result.resize(component_i);
  for (size_t component = 0; component < component_i; component++) {
    size_t parent_component = component;
    while (component_to_parent[parent_component] != parent_component) parent_component = component_to_parent[parent_component];

    result[parent_component].push_back({component_to_size[component], component_to_paritition[component]});
  }
}

void find_inefficient_super_components(                               
  const Context& context,
  const vec<vec<ComponentInfo>>& super_components,
  vec<vec<ComponentInfo>>& result
) {
  
  Bitset has_component_in_parititon;
  has_component_in_parititon.resize(context.partition.k);

  for (const vec<connected_components::ComponentInfo>& components_parent : super_components) {

    has_component_in_parititon.reset();
    bool is_inefficient = false;

    for (const connected_components::ComponentInfo& component : components_parent) {
      if (has_component_in_parititon.isSet((size_t) component.partition)) {
        is_inefficient = true;
        break;
      }

      has_component_in_parititon.set((size_t) component.partition);
    }

    if (is_inefficient) {
      result.push_back(components_parent);
    }
  }
}

namespace {
#define COMPUTE_COMPONENTS_PER_BLOCK(X) void compute_components_per_block(const X& phg, const Context& context, vec<vec<ConnectedComponent>>& result)
#define COMPUTE_COMPONENTS(X) void compute_components(const X& phg, const Context& context, vec<ConnectedComponent>& result)
#define COMPUTE_VC(X) void compute_super_components(const X& phg, const Context& context, const vec<vec<ConnectedComponent>>& components, vec<vec<ComponentInfo>>& result)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(COMPUTE_COMPONENTS_PER_BLOCK)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(COMPUTE_COMPONENTS)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(COMPUTE_VC)

}  // namespace connected_components
}  // namespace mt_kahypar
