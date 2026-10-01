
#include "mt-kahypar/utils/connected_component_stats.h"

namespace mt_kahypar::utils {

DebugState cc_debug;

template<typename PartitionedHypergraph>
double get_average_node_degree(const PartitionedHypergraph& phg) {
    size_t total_degree = 0;

    for (const HypernodeID& node : phg.nodes()) {
        total_degree += phg.incidentEdges(node).size();
    }

    return total_degree / phg.nodes().size();
}

template<typename PartitionedHypergraph>
size_t get_maximum_node_degree(const PartitionedHypergraph& phg) {
    size_t max_degree = 0;

    for (const HypernodeID& node : phg.nodes()) {
        if (max_degree < phg.incidentEdges(node).size()) {
            max_degree = phg.incidentEdges(node).size();
        }
    }

    return max_degree;
}

template<typename PartitionedHypergraph>
size_t get_minimum_node_degree(const PartitionedHypergraph& phg) {
    size_t min_degree = get_maximum_node_degree(phg);

    for (const HypernodeID& node : phg.nodes()) {
        if (min_degree > phg.incidentEdges(node).size()) {
            min_degree = phg.incidentEdges(node).size();
        }
    }

    return min_degree;
}


}