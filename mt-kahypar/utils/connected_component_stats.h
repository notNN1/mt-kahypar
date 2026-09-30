#pragma once

#include "mt-kahypar/definitions.h"


namespace mt_kahypar::utils {

struct DebugState {
    vec<vec<size_t>> components_per_partition_per_level;

    vec<std::tuple<size_t, size_t, size_t>> eNodes_to_blockedENodes_to_allNodes_per_rebuild;
};


extern DebugState cc_debug;

} // mtkahypar::utils