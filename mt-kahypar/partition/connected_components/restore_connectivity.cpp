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
  vec<vec<ComponentInfo>>& super_components
) {
  // find components
  vec<vec<ComponentInfo>> infos;
  vec<size_t> component_weight;
  infos.resize(context.partition.k);

  ////// just throw components over
  for (vec<connected_components::ComponentInfo>& super_component : super_components) {

    std::sort(super_component.begin(), super_component.end(),
        [](const connected_components::ComponentInfo& a, const connected_components::ComponentInfo& b) {
            return a.weight < b.weight;
        }
    );

    for (const connected_components::ComponentInfo& component : super_component) {
      // find partition to move component to

      PartitionID move_to = kInvalidPartition;
      for (const HypernodeID& node : component.nodes) {
        for (const HyperedgeID& edge : phg.incidentEdges(node)) {
          for (const HypernodeID& incident_hn : phg.pins(edge)) {
            if (phg.partID(incident_hn) != component.partition) {
              move_to = phg.partID(incident_hn);
            }
          }
        }
      }
      
      if (move_to == kInvalidPartition) {
        continue;
      }

      if (phg.partWeight(move_to) + component.weight < (phg.totalWeight() / phg.k()) * (1.0 + 3 * context.partition.epsilon)) {
        LOG << "hapened";
        for (const HypernodeID& node : component.nodes) {
          phg.changeNodePart(node, component.partition, move_to, DynamicConnectivityStrategy::do_nothing);
        }
      }
    }
  }
}


namespace {
#define RESTORE_CONNECTIVITY(X) void restore_connectivity(X& phg, const Context& context, vec<vec<ComponentInfo>>& super_components)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(RESTORE_CONNECTIVITY)

}  // namespace connected_components
}  // namespace mt_kahypar