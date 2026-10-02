#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/partition/refinement/rebalancing/advanced_rebalancer.h"

#pragma once

namespace mt_kahypar {
namespace connected_components {

using Bitset = mt_kahypar::ds::Bitset;
using Time   = size_t;

using Hypergraph = ds::StaticHypergraph;
using Factory = Hypergraph::Factory;
using Rebalancer = AdvancedRebalancer<GraphAndGainTypes<StaticHypergraphTypeTraits, Km1GainTypes>>;
using GainCache = Km1GainTypes::GainCache;
using PartitionedHypergraph = ds::PartitionedHypergraph<mt_kahypar::ds::StaticHypergraph, mt_kahypar::ds::ConnectivityInfo>;

class Test {
public:
    void test_rebalancer();
};


} // mt_kahypar
} // connected_components