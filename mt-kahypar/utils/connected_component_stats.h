#pragma once

#include "mt-kahypar/definitions.h"


namespace mt_kahypar::utils {

enum NodeType {
    normal          = 0,
    articulation    = 1,
    interface       = 2
};

struct PackedComponentInfo {
    size_t total_weight;
    NodeType type;                              
};

struct IPState {
    vec<PackedComponentInfo> packed_components_tarjan;  // singleton intialization
};


struct StatsPerReset {
    // refinement st rebuilds
    size_t st_blocked_enodes;
    size_t st_enodes;

    // anker nodes
    size_t recomputations;
};


struct StatePerLevel {
    size_t total_nodes;

    // connectivity
    vec<size_t> components_per_partition;

    vec<StatsPerReset> stats_per_reset;
};

class DebugState {
private:
    IPState                 ip_state;
    vec<StatePerLevel>      state_per_level;
public:

    DebugState() {
        this->add_new_level();
        this->add_new_reset();
    }

    // refinement
    void add_new_level() {
        state_per_level.push_back({ });
    }

    void set_total_nodes(size_t total_nodes) {
        state_per_level.back().total_nodes = total_nodes;
    }

    void add_new_reset() {
        state_per_level.back().stats_per_reset.push_back({ });
    }

    void s_initialize_components_per_partition(vec<size_t> components_per_partition) {
        if (this->state_per_level.back().components_per_partition.size() != 0) {
            return;
        }

        this->state_per_level.back().components_per_partition = components_per_partition;
    };

    void add_st_rebuild_stats(size_t st_blocked_enodes, size_t st_enodes) {
        this->state_per_level.back().stats_per_reset.back().st_blocked_enodes = st_blocked_enodes;
        this->state_per_level.back().stats_per_reset.back().st_enodes = st_enodes;
    }

    void increse_anker_rebuilds() {
        this->state_per_level.back().stats_per_reset.back().recomputations++;
    }

    // ip
    void s_initialize_components_tarjan(vec<PackedComponentInfo> components) {
        if (this->ip_state.packed_components_tarjan.size() != 0) {
            return;
        }

        this->ip_state.packed_components_tarjan = components;
    };
};




extern DebugState cc_debug;

} // mtkahypar::utils