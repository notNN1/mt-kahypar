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

#include "mt-kahypar/partition/connected_components/anker_node.h"

namespace mt_kahypar {
namespace connected_components {


template<typename PartitionedHypergraph>
void AnkerNodes<PartitionedHypergraph>::initialize(
    const PartitionedHypergraph& phg,
    const vec<uint32_t>& node_priority
) {
    vec<vec<HypernodeID>> node_to_partition(phg.initialNumNodes(), vec<HypernodeID>(phg.k(), kInvalidHypernode));

    vec<HypernodeID> nodes_for_partition(phg.k(), kInvalidHypernode);
    vec<uint32_t> priority(phg.k(), 0);

    uint32_t count = 0;

    for (const HyperedgeID& he : phg.edges()) {
        
        
        for (const HypernodeID& incident_hn : phg.pins(he)) {
            PartitionID part = phg.partID(incident_hn);

            if (nodes_for_partition[part] != kInvalidHypernode && node_priority[nodes_for_partition[part]] >= node_priority[incident_hn]) {
                continue;
            }

            count++;
            nodes_for_partition[part] = incident_hn;
        }

        if (count >= 2) {
            for (uint32_t i = 0; i < phg.k(); i++) {
                if (nodes_for_partition[i] == kInvalidHypernode) {
                    continue;
                }

                for (const HypernodeID& incident_hn : phg.pins(he)) {
                    if (node_to_partition[incident_hn][i] != kInvalidHypernode) {
                        continue;
                    }

                    node_to_partition[incident_hn][i] = nodes_for_partition[i];
                }
            }
        }

        nodes_for_partition.clear();
    }

    this->node_to_partition = node_to_partition;
}

template<typename PartitionedHypergraph>
HypernodeID AnkerNodes<PartitionedHypergraph>::find_node_in_partition(
    const PartitionedHypergraph& phg,
    const HypernodeID& hn,
    const PartitionID& partition
) {
    HypernodeID anker_node = this->node_to_partition[hn][partition];
    if (partition == phg.partID(anker_node)) {
        return anker_node;
    }

    // else find new anker node
    for (const HyperedgeID& he : phg.incidentEdges(hn)) {
        for (const HypernodeID& incident_hn : phg.pins(he)) {
            if (phg.partID(incident_hn) == partition) {
                this->node_to_partition[hn][partition] = incident_hn;
                return incident_hn;
            }
        }
    }
}

INSTANTIATE_CLASS_WITH_PARTITIONED_HG(ConnectivityFacade)

}  // namespace connected_components
}  // namespace mt_kahypar

