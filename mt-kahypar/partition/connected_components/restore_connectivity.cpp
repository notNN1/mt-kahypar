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
  const vec<vec<ComponentInfo>>& super_components,
  vec<ComponentInfo>& result,
  const Context& context,
  const double total_weight_ratio
) {
  // find components
  vec<vec<ComponentInfo>> infos;
  vec<size_t> component_weight;
  infos.resize(context.partition.k);

  ////// just throw components over
  for (const vec<connected_components::ComponentInfo>& super_component : super_components) {

    vec<size_t> components_per_partition = find_possible_partitions(super_component, phg, context);

    for (const connected_components::ComponentInfo& component : super_component) {
      
      
    }
  }
}

/*if (component.parititon_weight < total_weight_ratio * phg.totalWeight()) {
    for (const HypernodeID& node : component.nodes) {
        phg.setNodePart(node, selected_part);
    }
}*/

template<typename PartitionedHypergraph>
vec<size_t> find_possible_partitions(
  const vec<connected_components::ComponentInfo>& super_component,
  const PartitionedHypergraph& phg,
  const Context& context
) {

  vec<size_t> components_per_partition;
  components_per_partition.resize(phg.k());

  for (const connected_components::ComponentInfo& component : super_component) {
    components_per_partition[component.partition]++;
  }

  return components_per_partition;
}



namespace {
#define FIND_PARTITIONS(X) vec<size_t> find_possible_partitions(const vec<ComponentInfo>& super_component, const X& phg, const Context& context)
#define RESTORE_CONNECTIVITY(X) void restore_connectivity(X& phg, const vec<vec<ComponentInfo>>& super_components, vec<ComponentInfo>& result, const Context& context, const double total_weight_ratio)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(FIND_PARTITIONS)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(RESTORE_CONNECTIVITY)

}  // namespace connected_components
}  // namespace mt_kahypar