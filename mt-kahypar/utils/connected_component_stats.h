#pragma once

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/initial_partitioning/tarjan_initial_partitioner.h"


namespace mt_kahypar::utils {

struct DebugState {
    vec<vec<size_t>> components_per_partition_per_level;

    vec<std::tuple<size_t, size_t, size_t>> eNodes_to_blockedENodes_to_allNodes_per_rebuild;

    vec<mt_kahypar::PackedComponentInfo> packed_components_tarjan;
};

template<typename PartitionedHypergraph>
double get_average_node_degree(const PartitionedHypergraph& phg);

template<typename PartitionedHypergraph>
size_t get_maximum_node_degree(const PartitionedHypergraph& phg);

template<typename PartitionedHypergraph>
size_t get_minimum_node_degree(const PartitionedHypergraph& phg);


extern DebugState cc_debug;

} // mtkahypar::utils