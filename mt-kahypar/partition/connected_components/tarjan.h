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
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"

#include <algorithm>
#include <signal.h>

namespace mt_kahypar {
namespace connected_components {

using Bitset = mt_kahypar::ds::Bitset;
using Time   = size_t;

using Hypergraph = ds::StaticHypergraph;
using Factory = Hypergraph::Factory;

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
            this->initialized = true;
        }

        return is_ap.isSet((size_t) node);
    }

    void initialize(const PartitionedHypergraph& phg) {
        reset_fields(phg);

        vec<vec<HypernodeID>> hn_to_new_incident_nodes = phg.circular_edge_expansion();
        
        for (const HypernodeID& node : phg.nodes()) {
            if (disc[node] == 0) {
                if (dfs_recursive_wrapper(phg, node, hn_to_new_incident_nodes) > 1) {
                    is_ap.set((size_t) node);
                }
            }
        }

    }
private:
    size_t dfs_recursive_wrapper(
        const PartitionedHypergraph& phg, 
        const HypernodeID& node,
        const vec<vec<HypernodeID>>& hn_to_new_incident_nodes
    ) {
        return dfs_recursive(phg, node, node, hn_to_new_incident_nodes);
    }

    size_t dfs_recursive(
        const PartitionedHypergraph& phg, 
        const HypernodeID& node, 
        const HypernodeID& parent,
        const vec<vec<HypernodeID>>& hn_to_new_incident_nodes
    ) {
        size_t children = 0;
        low[node] = disc[node] = time;
        time++;

        for (const HypernodeID& incident_hn : hn_to_new_incident_nodes[node]) {
            if (incident_hn == parent) continue; // we don't want to go back through the same path.
                            // if we go back is because we found another way back
            if (disc[incident_hn] == 0) { // if V has not been discovered before
                children++;
                dfs_recursive(phg, incident_hn, node, hn_to_new_incident_nodes); // recursive DFS call

                if (node != parent && disc[node] <= low[incident_hn]) // condition #1
                    is_ap.set((size_t) node);

                low[node] = std::min(low[node], low[incident_hn]); // low[v] might be an ancestor of u
            } else // if v was already discovered means that we found an ancestor
                low[node] = std::min(low[node], disc[incident_hn]); // finds the ancestor with the least discovery time
        }
        
        return children;
    }


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

        Bitset node_colored;
        node_colored.resize(phg.initialNumNodes());

        bool backtracking       = false;


        vec<vec<vec<HypernodeID>>> he_hn_to_new_edge = phg.circular_edge_expansion();

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

                    for (const HypernodeID& incident_hn : he_hn_to_new_edge[he][current]) {

                        if (incident_hn == parent) {
                            continue;
                        }

                        if (disc[incident_hn] > 0) {
                            low[current] = low[current] < disc[incident_hn] ? low[current] : disc[incident_hn];
                            continue;
                        }

                        if (current == node) {
                            children++;
                        }


                        node_colored.set((size_t) incident_hn);

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

public: 
    void test_tarjan() {
        LOG << "Begin tarjan test";
        test_tarjan_path_middle();
        test_tarjan_leaf();
        test_tarjan_cycle();
        //test_circular_expansion();
        //test_contraction();
        LOG << "End tarjan test";
    }

    void test_circular_expansion() {

        auto hg = Factory::construct(
            5,
            4,
            {
                {0,1},
                {1,2},
                {2,3},
                {3,4}
            },
            nullptr,
            nullptr,
            true
        );


        PartitionedHypergraph phg(
            2,     // two partitions
            hg
        );

        phg.setNodePart(0, 0);
        phg.setNodePart(1, 0);
        phg.setNodePart(2, 0);
        phg.setNodePart(3, 1);
        phg.setNodePart(4, 1);

        vec<vec<vec<HypernodeID>>> he_hn_to_new_edge = phg.circular_edge_expansion();

        for (const HyperedgeID& he : phg.edges()) {
            for (const HypernodeID& hn : phg.pins(he)) {
                LOG << "HN: " << hn;
                for (const HypernodeID& incident_hn : he_hn_to_new_edge[he][hn]) {
                    LOG << "incident_hn: " << incident_hn;
                }
            }
            LOG << "############";
        }
    }

    void test_contraction() {
        auto hg = Factory::construct(
            5,
            4,
            {
                {0,1},
                {1,2},
                {2,3},
                {3,4}
            },
            nullptr,
            nullptr,
            true
        );

        PartitionedHypergraph phg(2, hg);

        mt_kahypar::ds::StaticHypergraph& hg_p = phg.hypergraph();

        parallel::scalable_vector<HypernodeID> communities(5);
        communities[0] = 0;
        communities[1] = 0;
        communities[2] = 0;
        communities[3] = 1;
        communities[4] = 1;

        mt_kahypar::ds::StaticHypergraph new_hg = hg_p.contract(communities, true);
        PartitionedHypergraph phg_n(2, new_hg);

        LOG << "phg incident net: " << new_hg.incidentEdges(0).size();
        LOG << "initial num nodes: " << new_hg.initialNumNodes();
    }

    void test_tarjan_path_middle() {
        auto hg = Factory::construct(
            5,
            4,
            {
                {0,1},
                {1,2},
                {2,3},
                {3,4}
            },
            nullptr,
            nullptr,
            true
        );

        PartitionedHypergraph phg(2, hg);

        phg.setNodePart(0,0);
        phg.setNodePart(1,0);
        phg.setNodePart(2,0);
        phg.setNodePart(3,1);
        phg.setNodePart(4,1);

        // Removing 2 disconnects {0,1} from {3,4}
        Tarjan t;

        if (!t.is_articulation_point(phg, 2)) {
            LOG << "Should be articulaton point!";
        }
    }

    void test_tarjan_leaf() {
        auto hg = Factory::construct(
            5,
            4,
            {
                {0,1},
                {1,2},
                {2,3},
                {3,4}
            },
            nullptr,
            nullptr,
            true
        );

        PartitionedHypergraph phg(2, hg);

        phg.setNodePart(0,0);
        phg.setNodePart(1,0);
        phg.setNodePart(2,0);
        phg.setNodePart(3,0);
        phg.setNodePart(4,0);

        Tarjan t;

        if (t.is_articulation_point(phg, 0)) {
            LOG << "A leaf should never be an articulation point";
        }
    }

    void test_tarjan_cycle() {
        auto hg = Factory::construct(
            4,
            6,
            {
                {0,1},
                {1,2},
                {2,0},
                {0,3},
                {1,3},
                {2,3},
            },
            nullptr,
            nullptr,
            true
        );

        PartitionedHypergraph phg(2, hg);

        phg.setNodePart(0,0);
        phg.setNodePart(1,0);
        phg.setNodePart(2,0);
        phg.setNodePart(3,1);

        if (!phg.canMoveVertex(DynamicConnectivityStrategy::bfs, 0, 1)) {
            LOG << "1There should be no articulation points in the cycle";
        }
        if (!phg.canMoveVertex(DynamicConnectivityStrategy::bfs, 1, 1)) {
            LOG << "2There should be no articulation points in the cycle";
        }
        if (!phg.canMoveVertex(DynamicConnectivityStrategy::bfs, 2, 1)) {
            LOG << "3There should be no articulation points in the cycle";
        }
    }
};

}  // namespace connected_components
}  // namespace mt_kahypar