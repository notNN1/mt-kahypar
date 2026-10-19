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
#include "mt-kahypar/partition/connected_components/restore_connectivity.h"
#include "mt-kahypar/definitions.h"

#include <algorithm>

namespace mt_kahypar {
namespace connected_components {

using Bitset = mt_kahypar::ds::Bitset;

template<typename PartitionedHypergraph>
void restore_connectivity(
  PartitionedHypergraph& phg,
  const Context& context,
  vec<vec<ComponentInfo>>& super_components,
  const uint32_t level,
  const vec<vec<ConnectedComponent>>& components_per_block
) {

  vec<uint32_t> components_per_partition;

  for (vec<ConnectedComponent> components : components_per_block) {
    components_per_partition.push_back(components.size());
  }


  ////// just throw components over
  for (vec<connected_components::ComponentInfo>& super_component : super_components) {

    std::sort(super_component.begin(), super_component.end(),
        [](const connected_components::ComponentInfo& a, const connected_components::ComponentInfo& b) {
            return a.weight > b.weight;
        }
    );

    for (const connected_components::ComponentInfo& component : super_component) {
      // find partition to move component to
      if (components_per_partition[component.partition] <= 1) {
        continue;
      }

      PartitionID move_to             = kInvalidPartition;
      HyperedgeWeight best_imbalance  = std::numeric_limits<int32_t>::max();

      Bitset node_colored;
      node_colored.resize(phg.initialNumNodes());

      Bitset edge_colored;
      edge_colored.resize(phg.initialNumEdges());

      for (const HypernodeID& node : component.nodes) {

        if (node_colored.isSet((size_t) node)) {
          continue;
        }

        node_colored.set((size_t) node);

        for (const HyperedgeID& edge : phg.incidentEdges(node)) {

          if (edge_colored.isSet((size_t) edge)) {
            continue;
          }

          edge_colored.set((size_t) edge);

          for (const HypernodeID& incident_hn : phg.pins(edge)) {
            if (node_colored.isSet((size_t) incident_hn)) {
              continue;
            }

            node_colored.set((size_t) incident_hn);

            if (phg.partID(incident_hn) != component.partition) {
              PartitionID part_id = phg.partID(incident_hn);

              if (phg.partWeight(part_id) + component.weight > (phg.totalWeight() / phg.k()) * (1.0 + (level * 2 + 1) * context.partition.epsilon)) {
                continue;
              }

              if (phg.partWeight(part_id) + component.weight < best_imbalance) {
                best_imbalance = phg.partWeight(part_id) + component.weight;
                move_to = part_id;
              }

            }
          }
        }
      }
      
      if (move_to == kInvalidPartition) {
        continue;
      }

      LOG << "moved component";
      for (const HypernodeID& node : component.nodes) {
        phg.changeNodePart(node, component.partition, move_to, DynamicConnectivityStrategy::do_nothing);
      }

    }
  }
}


namespace {
#define RESTORE_CONNECTIVITY(X) void restore_connectivity(X& phg, const Context& context, vec<vec<ComponentInfo>>& super_components, const uint32_t level, const vec<vec<ConnectedComponent>>& components_per_block)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(RESTORE_CONNECTIVITY)

}  // namespace connected_components
}  // namespace mt_kahypar