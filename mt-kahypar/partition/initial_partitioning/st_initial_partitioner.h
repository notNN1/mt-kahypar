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
#include <deque>

namespace mt_kahypar {

template<typename TypeTraits>
class STInitialPartitioner : public IInitialPartitioner {

    static constexpr bool debug = false;

    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    using Bitset                = typename mt_kahypar::ds::Bitset;
    using ConnectedComponent    = typename mt_kahypar::connected_components::ConnectedComponent;

public:
    STInitialPartitioner(
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
    void partitionImpl() final;

    void calculate_spanning_tree(
        const PartitionedHypergraph& hg,
        ConnectedComponent& component,
        vec<HypernodeID>& hn_to_parent,
        vec<vec<HypernodeID>>& hn_to_children,
        vec<size_t>& subtree_size,
        vec<size_t>& covered
    ); 

    std::pair<HypernodeID, size_t> find_best_node_to_split(
        const ConnectedComponent& component,
        const vec<size_t>& subtree_size,
        vec<size_t>& covered,
        const size_t& target
    );

    void assign_subtree_of_hn(
        PartitionedHypergraph& hg,
        vec<vec<HypernodeID>>& hn_to_children,
        vec<size_t>& subtree_size,
        PartitionID partition,
        vec<size_t>& covered,
        size_t& current_split_number,
        HypernodeID hn
    );

    InitialPartitioningDataContainer<TypeTraits>& _ip_data;
    const Context& _context;
    std::mt19937 _rng;
    const int _tag;
};

} // namespace mt_kahypar
