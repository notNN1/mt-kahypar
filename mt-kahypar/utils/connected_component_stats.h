#pragma once

#include <mutex>
#include <utility>

#include "mt-kahypar/definitions.h"


namespace mt_kahypar::utils {

enum NodeType {
    normal          = 0,
    articulation    = 1,
    interface       = 2
};

struct RegionInfo {
    size_t total_weight;
    NodeType type;                              
};

struct BFSStats {
    size_t assigned_vertices = 0;
    // Includes the initial spanning-tree construction, as does st_bfs_max_recalculations.
    size_t recalculations = 0;
    vec<size_t> assigned_vertices_per_recalculation;
};

struct BFSAverageStats {
    double assigned_vertices = 0.0;
    double recalculations = 0.0;
};

struct IPState {
    vec<RegionInfo> regions_tarjan;  // singleton intialization
    RegionInfo chosen_region;
    vec<BFSStats> bfs_runs;
};


struct StatsPerReset {
    // refinement st rebuilds
    size_t st_blocked_enodes;
    size_t st_enodes;

    // anker nodes
    size_t recomputations;
    size_t node_moves;
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
    mutable std::mutex      bfs_stats_mutex;
public:

    DebugState() {
        add_new_level();
        add_new_reset();
    }

    // refinement
    void add_new_level() {
        state_per_level.push_back({ });
        add_new_reset();
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

    void increse_node_changes() {
        this->state_per_level.back().stats_per_reset.back().node_moves++;
    }

    // ip
    void add_bfs_stats(BFSStats stats) {
        stats.recalculations = stats.assigned_vertices_per_recalculation.size();
        std::lock_guard<std::mutex> lock(bfs_stats_mutex);
        ip_state.bfs_runs.push_back(std::move(stats));
    }

    BFSAverageStats average_bfs_stats() const {
        std::lock_guard<std::mutex> lock(bfs_stats_mutex);
        BFSAverageStats averages;
        if (ip_state.bfs_runs.empty()) {
            return averages;
        }

        for (const BFSStats& stats : ip_state.bfs_runs) {
            averages.assigned_vertices += static_cast<double>(stats.assigned_vertices);
            averages.recalculations += static_cast<double>(stats.recalculations);
        }
        const double runs = static_cast<double>(ip_state.bfs_runs.size());
        averages.assigned_vertices /= runs;
        averages.recalculations /= runs;
        return averages;
    }

    void s_initialize_regions_tarjan(vec<RegionInfo> regions, RegionInfo chosen_region) {
        if (this->ip_state.regions_tarjan.size() != 0) {
            return;
        }

        this->ip_state.regions_tarjan = regions;
        this->ip_state.chosen_region = chosen_region;
    };
};




extern DebugState cc_debug;

} // mtkahypar::utils
