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

#pragma once

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/connected_components/compute_components.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/bitset.h"
#include "mt-kahypar/datastructures/dynamic_connectivity_datastructures.h"

namespace mt_kahypar {
namespace connected_components {


template<typename PartitionedHypergraph>
class BFSConnectivity {
public:
    bool moveVertex(
        const PartitionedHypergraph& phg, 
        const HypernodeID& hn
    ) {

        // compute components with hn set, so we can't move over hn
        int components = 0;

        Bitset node_colored;
        node_colored.resize(phg.initialNumNodes());
        node_colored.set((size_t) hn);

        Bitset edge_colored;
        edge_colored.resize(phg.initialNumEdges());
        
        std::queue<HypernodeID> node_queue;
        PartitionID current_partition = phg.partID(hn);

        for (const HypernodeID& hn_ : phg.nodes()) {
            if (node_colored.isSet((size_t) hn_) || phg.partID(hn_) != current_partition) {
                continue;
            }

            node_queue.push(hn_);
            node_colored.set((size_t) hn_);

            edge_colored.reset();

            while (node_queue.size() > 0) {
                HypernodeID current = node_queue.front();
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

                        node_queue.push(incident_hn);
                        node_colored.set((size_t) incident_hn);
                    }
                }
            }
            
            components++;
            if (components > 1) {
                return false;
            }
        }
        return true;
    };
};

}  // namespace connected_components
}  // namespace mt_kahypar