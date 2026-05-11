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
        node_colored.resize(phg.initialNumnodes());
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

        
    void rotate_zick(
        vec<Node>& nodes,
        int32_t x,
        int32_t p,
        int32_t B
    ) {
        nodes[x].parent         = -1;
        nodes[x].right          = p;
        nodes[x].splay_tree     = nodes[p].splay_tree;

        nodes[p].parent         = x;
        nodes[p].left           = B;
        nodes[p].splay_tree     = -1;
    };

        
    void rotate_zack(
        vec<Node>& nodes,
        int32_t x,
        int32_t p,
        int32_t B
    ) {
        nodes[x].parent         = -1;
        nodes[x].left           = p;
        nodes[x].splay_tree     = nodes[p].splay_tree;

        nodes[p].parent         = x;
        nodes[p].right          = B;
        nodes[p].splay_tree     = -1;
    };

        
    void rotate_zick_zick(
        vec<Node>& nodes,
        int32_t x,
        int32_t p,
        int32_t g,
        int32_t gg,
        int32_t B,
        int32_t C
    ) {
        nodes[x].right          = p;
        nodes[x].parent         = gg;
        nodes[x].splay_tree     = nodes[g].splay_tree;

        nodes[p].right    = g;
        nodes[p].parent   = x;
        nodes[p].left     = B;
        
        nodes[g].left           = C;
        nodes[g].parent         = p;
        nodes[g].splay_tree     = -1;
    };

        
    void rotate_zack_zack(
        vec<Node>& nodes,
        int32_t x,
        int32_t p,
        int32_t g,
        int32_t gg,
        int32_t B,
        int32_t C
    ) {
        nodes[x].left           = p;
        nodes[x].parent         = gg;
        nodes[x].splay_tree     = nodes[g].splay_tree;

        nodes[p].left     = g;
        nodes[p].parent   = x;
        nodes[p].right    = B;

        nodes[g].right          = C;
        nodes[g].parent         = p;
        nodes[g].splay_tree     = -1;
    };
    

        
    void rotate_zick_zack(
        vec<Node>& nodes,
        int32_t x,
        int32_t p,
        int32_t g,
        int32_t gg,
        int32_t B,
        int32_t C
    ) {
        nodes[x].parent         = gg;
        nodes[x].left           = p;
        nodes[x].right          = g;
        nodes[x].splay_tree     = nodes[g].splay_tree;

        nodes[p].parent   = x;
        nodes[p].right    = B;

        nodes[g].parent         = x;
        nodes[g].left           = C;
        nodes[g].splay_tree     = -1;
    };

        
    void rotate_zack_zick(
        vec<Node>& nodes,
        int32_t x,
        int32_t p,
        int32_t g,
        int32_t gg,
        int32_t B,
        int32_t C
    ) {
        nodes[x].parent         = gg;
        nodes[x].right          = p;
        nodes[x].left           = g;
        nodes[x].splay_tree     = nodes[g].splay_tree;

        nodes[p].parent   = x;
        nodes[p].left     = B;

        nodes[g].parent         = x;
        nodes[g].right          = C;
        nodes[g].splay_tree     = -1;
    };

        
    static void splay(
        vec<Node>& nodes,
        Node x
    ) {
        int32_t p    = x.parent;

        if (p == -1) {
            return;
        }

        int32_t g    = nodes[p].parent;

        // zick / zack
        if (g == -1) {

            if (nodes[p].left == x.self) {
                rotate_zick(nodes, x.self, p, x.right);
            } 
            else {
                rotate_zack(nodes, x.self, p, x.left);
            }

            return;
        }

        int32_t gg   = nodes[g].parent;

        if (nodes[p].left == x.self) {
            if (nodes[g].left == p) {
                rotate_zick_zick(
                    nodes,
                    x.self,
                    p,
                    g,
                    gg,
                    x.right,
                    nodes[p].right
                );
            }
            else {
                rotate_zick_zack(
                    nodes,
                    x.self,
                    p,
                    g,
                    gg,
                    x.right,
                    x.left
                );
            }
        }
        else {
            if (nodes[g].left == p) {
                rotate_zack_zick(
                    nodes,
                    x.self,
                    p,
                    g,
                    gg,
                    x.left,
                    x.right
                );
            }
            else {
               rotate_zack_zack(
                nodes,
                    x.self,
                    p,
                    g,
                    gg,
                    x.left,
                    nodes[p].left
                );
            }
        }

        splay(nodes, x);
    };


    void LinkCutTree::expose(
        Node& u
    ) {
        splay(this->nodes, u);
        SplayTree current = paths[u.splay_tree];
        int32_t current_index = u.splay_tree;

        Node& lower_node = u;
        Node& upper_node = u;

        
        while (current.connected_to != -1) {
            upper_node = nodes[current.connected_to];
            splay(this->nodes, upper_node);

            // remove current splay tree
            this->paths.erase(this->paths.begin() + current_index);

            // connect to upper splay tree and add new splay tree for cut splay tree

            lower_node.splay_tree = -1;
            lower_node.parent = upper_node.self;

            if (upper_node.right != -1) {
                this->paths.push_back(SplayTree {upper_node.self, upper_node.right});
                this->nodes[upper_node.right].splay_tree = this->paths.size() - 1;
                this->nodes[upper_node.right].parent = -1;
            }

            upper_node.right = lower_node.self;

            // find the next splay tree and repeat
            current = paths[upper_node.splay_tree];
            current_index = upper_node.splay_tree;

            lower_node = upper_node;
        }
    }


namespace {
    // This macro defines how to refer to your class for a specific type X
    #define BFS_CONNECTIVITY(X) BFSConnectivity<X>
}

// This tells the compiler: "Build BFSConnectivity for all Hypergraph types"
INSTANTIATE_CLASS_WITH_PARTITIONED_HG(BFSConnectivity)

}  // namespace ds
}  // namespace mt_kahypar
