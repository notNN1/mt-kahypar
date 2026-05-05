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

#include "mt-kahypar/datastructures/dynamic_connectivity_datastructures.h"
#include "mt-kahypar/partition/connected_components/compute_components.h"

#include "mt-kahypar/definitions.h"
#include <tbb/task_group.h>

namespace mt_kahypar {
namespace ds {

    using Bitset = mt_kahypar::ds::Bitset;
    using ConnectedComponent = mt_kahypar::connected_components::ConnectedComponent;

    template<typename PartitionedHypergraph>    
    bool BFSConnectivity<PartitionedHypergraph>::moveVertex(
        const PartitionedHypergraph& phg, 
        const Context& context,
        HypernodeID hn, 
        PartitionID from
    ) {
        (void) context;

        // compute components with hn set, so we can't move over hn
        int components = 0;

        Bitset node_colored;
        node_colored.resize(phg.initialNumNodes());
        node_colored.set((size_t) hn);
        
        std::queue<HypernodeID> node_queue;
        PartitionID current_partition = from;

        for (const HypernodeID& hn_ : phg.nodes()) {
            if (node_colored.isSet((size_t) hn_) || phg.partID(hn_) != current_partition) {
                continue;
            }

            node_queue.push(hn_);
            node_colored.set((size_t) hn_);

            while (node_queue.size() > 0) {
                HypernodeID current = node_queue.front();
                node_queue.pop();

                for (const HyperedgeID& he : phg.incidentEdges(current)) {
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
    


namespace {
    // This macro defines how to refer to your class for a specific type X
    #define BFS_CONNECTIVITY(X) BFSConnectivity<X>
}

// This tells the compiler: "Build BFSConnectivity for all Hypergraph types"
INSTANTIATE_CLASS_WITH_PARTITIONED_HG(BFSConnectivity)

}  // namespace ds
}  // namespace mt_kahypar
