/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2026 Simon Lang <simonlang2>@kit.edu
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

#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/initial_partitioning_data_container.h"
#include "mt-kahypar/datastructures/bitset.h"
#include "mt-kahypar/partition/connected_components/compute_components.h"
#include <set>

namespace mt_kahypar {

template<typename TypeTraits>
class TarjanInitialPartitioner : public IInitialPartitioner {

    static constexpr bool debug = false;

    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    using Bitset                = typename mt_kahypar::ds::Bitset;
    using ConnectedComponent    = typename mt_kahypar::connected_components::ConnectedComponent;
    using PackedComponentID     = uint32_t;   

public:
    TarjanInitialPartitioner(
        const InitialPartitioningAlgorithm,
        ip_data_container_t* ip_data,
        const Context& context,
        const int seed,
        const int tag
    ) :
    _ip_data(ip::to_reference<TypeTraits>(ip_data)),
    _context(context),
    _rng(seed),
    _tag(tag) { }

 private:
    // starts the partitioning
    void partitionImpl() final;

    // for tarjan
    enum NodeType {
        normal          = 0,
        articulation    = 1,
        interface       = 2
    };

    vec<NodeType> hn_to_node_type;                  // maps each node to its type

    struct PackedComponentInfo {
        PackedComponentID id;
        size_t total_weight;
        NodeType type;                              // should only be normal or articulation
        vec<HypernodeID> nodes;                     // contains all nodes in the component -- does not include the articulation point
        std::set<HypernodeID> connected_to;              // only interface nodes
    };

    // Do BFS while not moving over nodes of other types and collect them into a packed component
    // Does not compact interface nodes and instead gives each one a separate PackedComponent
    void compact_regions(
        const PartitionedHypergraph& hypergraph,
        vec<PackedComponentID>& vertex_to_packed_component,
        vec<PackedComponentInfo>& packed_component_info,
        connected_components::Tarjan<PartitionedHypergraph>& tarjan
    );

    // calculates a spanning tree over the packed components
    void calculate_master_spanning_tree(
        const vec<PackedComponentID>& vertex_to_packed_component,
        const vec<PackedComponentInfo>& packed_component_info,
        vec<PackedComponentID>& component_to_parent,
        vec<size_t>& subtree_size
    );

    void find_biggest_leaf() {
        
    };

    void find_farthest_component_from_component(
        const vec<PackedComponentInfo>& packed_component_info,
        const PackedComponentInfo& starter_component,
        PackedComponentInfo& end_component
    );

    void build_spanning_tree_from_node_in_direction() {

    };


    InitialPartitioningDataContainer<TypeTraits>& _ip_data;
    const Context& _context;
    std::mt19937 _rng;
    const int _tag;
};

} // namespace mt_kahypar
