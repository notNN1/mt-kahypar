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
#include "mt-kahypar/utils/connected_component_stats.h"
#include <set>

namespace mt_kahypar {
    
using PackedComponentID     = uint32_t;   

struct RegionInfo {
    PackedComponentID id;
    size_t total_weight;
    utils::NodeType type;
    vec<HypernodeID> nodes;                     // contains all nodes in the component -- does not include the articulation point
    std::set<HypernodeID> connected_to;              // only interface nodes
};

template<typename TypeTraits>
class TarjanInitialPartitioner : public IInitialPartitioner {

    static constexpr bool debug = false;

    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    using Bitset                = typename mt_kahypar::ds::Bitset;
    using ConnectedComponent    = typename mt_kahypar::connected_components::ConnectedComponent;

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

    vec<utils::NodeType> hn_to_node_type;                  // maps each node to its type

    // Do BFS while not moving over nodes of other types and collect them into a packed component
    // Does not compact interface nodes and instead gives each one a separate PackedComponent
    void compact_regions(
        const PartitionedHypergraph& hypergraph,
        vec<PackedComponentID>& vertex_to_region,
        vec<RegionInfo>& region_info,
        connected_components::Tarjan<PartitionedHypergraph>& tarjan
    );

    // calculates a spanning tree over the packed components
    void calculate_master_spanning_tree(
        const vec<PackedComponentID>& vertex_to_region,
        const vec<RegionInfo>& region_info,
        vec<size_t>& subtree_size,
        vec<PackedComponentID>& component_to_parent,
        vec<PackedComponentID>& heads
    );

    void calculate_communities(
        const vec<RegionInfo>& region_info,
        vec<PackedComponentID>& component_to_parent,
        parallel::scalable_vector<HypernodeID>& communities,
        const vec<size_t>& subtree_size,
        const vec<PackedComponentID> heads
    );

    void find_farthest_component_from_component(
        const vec<RegionInfo>& region_info,
        const RegionInfo& starter_component,
        RegionInfo& end_component
    );

    void build_spanning_tree_from_node_in_direction() {

    };


    InitialPartitioningDataContainer<TypeTraits>& _ip_data;
    const Context& _context;
    std::mt19937 _rng;
    const int _tag;
};

} // namespace mt_kahypar
