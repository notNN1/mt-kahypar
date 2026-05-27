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
#include <list>

namespace mt_kahypar {
namespace ds {

using Bitset = mt_kahypar::ds::Bitset;
using ConnectedComponent = mt_kahypar::connected_components::ConnectedComponent;
using ComponentID = uint32_t;

template<typename PartitionedHypergraph>
class BFSConnectivity {
public:
    bool moveVertex(
        const PartitionedHypergraph& phg, 
        const Context& context,
        HypernodeID hn
    );
};

template<typename PartitionedHypergraph>
class SpanningTreeConnectivity {
private:
    
    vec<HypernodeID> vertex_to_parent_compressed;
    vec<size_t> vertex_to_rank;

    struct Connection {
        typename std::list<Connection>::iterator iterator;
        HypernodeID node;

        Connection(HypernodeID node_) : node(node_) {}
    };

    std::vector<std::list<Connection>> connected_to;
    

    inline void try_connect_to_incident_without_connection(
        Bitset& has_connection_to_other_partition,
        const HypernodeID& hn,
        const HypernodeID& incident_hn
    );

    inline void try_connect_to_incident_with_connection(
        Bitset& has_connection_to_other_partition,
        const HypernodeID& hn,
        const HypernodeID& incident_hn
    );

    inline void connect_nodes(
        const HypernodeID& hn,
        const HypernodeID& incident_hn
    );

    inline bool is_same_component(
        const HypernodeID& hn1,
        const HypernodeID& hn2
    );

    int calc_parents(const PartitionedHypergraph& phg);
public:
    SpanningTreeConnectivity(
        const PartitionedHypergraph& phg,
        const Context& context
    );

    SpanningTreeConnectivity() = default;

    bool canMoveVertex(
        const Context& context,
        HypernodeID hn
    );

    void moveVertex(
        const PartitionedHypergraph& phg,
        const Context& context,
        HypernodeID hn,
        PartitionID to
    );

};

}  // namespace ds
}  // namespace mt_kahypar
