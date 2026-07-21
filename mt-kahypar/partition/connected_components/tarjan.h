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


#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/connected_components/compute_components.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/bitset.h"
#include "mt-kahypar/datastructures/dynamic_connectivity_datastructures.h"
#include <signal.h>

namespace mt_kahypar {
namespace connected_components {

using Bitset = mt_kahypar::ds::Bitset;
using Time   = size_t;

template<typename PartitionedHypergraph>
struct Tarjan {
private:

    vec<size_t>     disc;
    vec<size_t>     low;
    
    Time            time;

    Bitset          is_ap;
    bool            initialized = false;

public:
    bool is_articulation_point(
        const PartitionedHypergraph& phg,
        const HypernodeID& node
    ) {
        if (!initialized) {
            initialize(phg);
        }

        return is_ap.isSet((size_t) node);
    }

    void initialize(const PartitionedHypergraph& phg) {
        reset_fields(phg);
        
        for (const HypernodeID& node : phg.nodes()) {
            if (disc[node] == 0) {
                if (dfs(phg, node) > 1) {
                    is_ap.set((size_t) node);
                }
            }
        }

    }
private:
    size_t dfs(
        const PartitionedHypergraph& phg, 
        const HypernodeID& node
    ) {
        size_t children = 0;

        vec<std::pair<HypernodeID, HypernodeID>> queue;
        queue.push_back({node, node});

        vec<std::pair<HypernodeID, HypernodeID>> bt_queue;
        bt_queue.push_back({node, node});
        
        HypernodeID current;
        HypernodeID parent;

        bool backtracking       = false;

        while(queue.size() > 0) {
            
            current     = queue.back().first;
            parent      = queue.back().second;

            queue.pop_back();

            if (!backtracking) {
                low[current]    = time;
                disc[current]   = time;
                time++;

                backtracking = true;
                for (const HyperedgeID& he : phg.incidentEdges(current)) {
                    for (const HypernodeID& incident_hn : phg.pins(he)) {
                        if (incident_hn == parent || incident_hn == current) {
                            continue;
                        }

                        if (disc[incident_hn] != 0) {
                            low[current] = low[current] < disc[incident_hn] ? low[current] : disc[incident_hn];
                            continue;
                        }

                        if (current == node) {
                            children++;
                        }

                        queue.push_back({incident_hn, current});
                        bt_queue.push_back({incident_hn, current});
                        backtracking = false;
                    }
                }
            }
            else {
                HypernodeID current_bt = bt_queue.back().first;
                HypernodeID parent_bt  = bt_queue.back().second;

                while(current != current_bt) {
                    if (bt_queue.empty()) {
                        LOG << "empty";
                        while(true);
                    }

                    bt_queue.pop_back();
                    current_bt  = bt_queue.back().first;
                    parent_bt   = bt_queue.back().second;

                    if (disc[parent_bt] <= low[current_bt]) {
                        is_ap.set((size_t) parent_bt);
                    }

                    low[parent_bt] = low[parent_bt] < low[current_bt] ? low[parent_bt] : low[current_bt];
                }

                backtracking = false;
            }
            
        }

        return children;
    }

    void reset_fields(const PartitionedHypergraph& phg) {
        disc.clear();
        disc.resize(phg.initialNumNodes());

        low.clear();
        low.resize(phg.initialNumNodes());

        is_ap.reset();
        is_ap.resize(phg.initialNumNodes());

        time = 1;
    }
};

}  // namespace connected_components
}  // namespace mt_kahypar