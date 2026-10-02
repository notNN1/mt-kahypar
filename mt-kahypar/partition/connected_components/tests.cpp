#include "mt-kahypar/partition/connected_components/tests.h"
#include "mt-kahypar/partition/refinement/rebalancing/advanced_rebalancer.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/randomize.h"

#include <limits>
#include <memory>


namespace mt_kahypar {
namespace connected_components {

void Test::test_rebalancer() {
    using Rebalancer =
        AdvancedRebalancer<
            GraphAndGainTypes<StaticHypergraphTypeTraits, Km1GainTypes>>;

    GainCache gain_cache;
    Context context;

    context.partition.mode = Mode::direct;
    context.partition.epsilon = 0.01;
    context.partition.k = 2;
    context.partition.preset_type = PresetType::default_preset;
    context.partition.instance_type = InstanceType::hypergraph;
    context.partition.partition_type = PartitionedHypergraph::TYPE;
    context.partition.objective = Objective::km1;
    context.partition.gain_policy = GainPolicy::km1;
    context.partition.allow_empty_blocks = false;
    context.refinement.dynamic_connectivity.advanced_rebalancer_dynamic_connectivity_strategy = DynamicConnectivityStrategy::st;

    auto hg = Factory::construct(
        8, 7,
        {{0,1}, {1,2}, {2,3}, {3,4}, {4,5}, {5,6}, {6, 7}},
        nullptr, nullptr, true);

    context.setupPartWeights(hg.totalWeight());

    PartitionedHypergraph phg(2, hg);

    phg.setNodePart(0, 0);
    phg.setNodePart(1, 0);
    phg.setNodePart(2, 0);
    phg.setNodePart(3, 0);
    phg.setNodePart(4, 0);
    phg.setNodePart(5, 0);
    phg.setNodePart(6, 0);
    phg.setNodePart(7, 1);

    auto rebalancer =
        std::make_unique<Rebalancer>(
            hg.initialNumNodes(), context, gain_cache);

    mt_kahypar_partitioned_hypergraph_t phg_handle =
        utils::partitioned_hg_cast(phg);

    rebalancer->initialize(phg_handle);

    Metrics metrics;
    metrics.quality = metrics::quality(phg, context);
    metrics.imbalance = metrics::imbalance(phg, context);

    rebalancer->refine(
        phg_handle,
        {},
        metrics,
        std::numeric_limits<double>::max());


    for (const HypernodeID node : phg.nodes()) {
        LOG << "node: " << node
            << " partition: " << phg.partID(node);
    }

}


} // namespace connected_components
} // namespace mt_kahypar