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

namespace mt_kahypar {
namespace ds {

template<typename PartitionedHypergraph>
class BFSConnectivity {
public:
    bool moveVertex(
        const PartitionedHypergraph& phg, 
        const Context& context,
        HypernodeID hn, 
        PartitionID from
    );
};

struct Node {
    int32_t left;
    int32_t right;
    int32_t parent;
    int32_t self;

    size_t splay_tree;

    bool rev;

    int value;
    int sum;
};

struct SplayTree {
    int32_t connected_to;
    int32_t root_node;
};

class LinkCutTree {
private:
    vec<SplayTree> paths;
    vec<Node> nodes;
public:
    void expose(
        Node& u
    );

    void link(
        Node u,
        Node v
    );

    void cut(
        Node u,
        Node v
    );

    void findRoot(
        Node u
    );

    void pathQuery(
        Node u,
        Node v
    );

    void pathUpdate(
        Node u,
        Node v
    );
};

/*template<typename PartitionedHypergraph>
class HolmeRotenbergThorup {
public:
    bool moveVertex(
        const PartitionedHypergraph& phg, 
        const Context& context,
        HypernodeID hn, 
        PartitionID from
    );

    HolmeRotenbergThorup(
        const PartitionedHypergraph& phg, 
        const Context& context
    );

private:
    struct Node {
        uint32_t id,
        bool is_edge
    }

    void link(
        Node v, 
        Node w
    );

    void cut(
        Node v, 
        Node w
    );

    bool connected(
        Node v, 
        Node w
    );

    void cover(
        Node v, 
        Node w,
        int32_t cover_level
    );

    void uncover(
        Node v, 
        Node w,
        int32_t cover_level
    );

    int32_t coverLevel(
        Node v
    );

    int32_t coverLevel(
        Node v,
        Node w
    );

    int32_t minCoveredEdge(
        Node v
    );
    
    int32_t minCoveredEdge(
        Node v,
        Node w
    );

    void addLabel(
        Node v,
        std::string user_label,
        int32_t cover_level
    );

    void removeLabel(
        std::string user_label
    );

    void findFirstLabel(
        Node v,
        Node w,
        int32_t cover_level
    );

    FindSize(
        Node v,
        Node w,
        int32_t cover_level
    )

    Node meet(
        Node u,
        Node v,
        Node w
    );


}*/

}  // namespace ds
}  // namespace mt_kahypar
