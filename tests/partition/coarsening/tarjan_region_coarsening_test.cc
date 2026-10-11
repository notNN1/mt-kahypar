#include "gtest/gtest.h"

#include "tests/definitions.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/partition/connected_components/tarjan_region_coarsening.h"
#include "mt-kahypar/partition/metrics.h"

namespace mt_kahypar {
namespace {

using Hypergraph = ds::StaticHypergraph;
using PartitionedHypergraph = StaticHypergraphTypeTraits::PartitionedHypergraph;
using RegionCoarsening = connected_components::TarjanRegionCoarsening<PartitionedHypergraph>;

TEST(TarjanRegionCoarsening, PreservesWeightedPartitionAndCutWhenContractingAndProjecting) {
  const HyperedgeWeight edge_weights[] = {2, 5, 7};
  const HypernodeWeight node_weights[] = {10, 20, 30, 40};
  auto hg = ds::StaticHypergraphFactory::construct(4, 3, {{0, 1}, {1, 2}, {2, 3}},
                                                 edge_weights, node_weights, true);
  PartitionedHypergraph phg(2, hg);
  for (const HypernodeID node : phg.nodes()) phg.setNodePart(node, node / 2);
  parallel::scalable_vector<HypernodeID> communities(4, 0);
  RegionCoarsening::respectPartition(phg, communities);
  const auto original_cut = metrics::quality(phg, Objective::cut);
  auto coarse_hg = hg.contract(communities, true);
  ASSERT_EQ(coarse_hg.initialNumNodes(), 2);
  PartitionedHypergraph coarse_phg(2, coarse_hg);
  coarse_phg.setNodePart(communities[0], 0);
  coarse_phg.setNodePart(communities[2], 1);
  EXPECT_EQ(coarse_phg.partWeight(0), 30);
  EXPECT_EQ(coarse_phg.partWeight(1), 70);
  EXPECT_EQ(metrics::quality(coarse_phg, Objective::cut), original_cut);

  // A coarse region move must lift to all of its original vertices.
  coarse_phg.resetPartition();
  coarse_phg.setNodePart(communities[0], 1);
  coarse_phg.setNodePart(communities[2], 0);
  phg.resetPartition();
  for (const HypernodeID node : phg.nodes()) phg.setNodePart(node, coarse_phg.partID(communities[node]));
  EXPECT_EQ(phg.partWeight(0), 70);
  EXPECT_EQ(phg.partWeight(1), 30);
  EXPECT_EQ(metrics::quality(phg, Objective::cut), metrics::quality(coarse_phg, Objective::cut));
}

TEST(TarjanRegionCoarsening, DoesNotHideDisconnectedPiecesOfTheSameBlock) {
  auto hg = ds::StaticHypergraphFactory::construct(4, 4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
  PartitionedHypergraph phg(2, hg);
  for (const HypernodeID node : phg.nodes()) phg.setNodePart(node, node % 2);
  parallel::scalable_vector<HypernodeID> communities(4, 0);
  RegionCoarsening::respectPartition(phg, communities);
  EXPECT_NE(communities[0], communities[2]);
  EXPECT_NE(communities[1], communities[3]);
  auto coarse_hg = hg.contract(communities, true);
  EXPECT_EQ(coarse_hg.initialNumNodes(), 4);
}

TEST(TarjanRegionCoarsening, KeepsFixedVerticesAsSingletons) {
  auto hg = ds::StaticHypergraphFactory::construct(3, 2, {{0, 1}, {1, 2}});
  ds::FixedVertexSupport<Hypergraph> fixed(3, 2);
  fixed.setHypergraph(&hg);
  fixed.fixToBlock(0, 0);
  hg.addFixedVertexSupport(std::move(fixed));
  PartitionedHypergraph phg(2, hg);
  for (const HypernodeID node : phg.nodes()) phg.setNodePart(node, 0);
  parallel::scalable_vector<HypernodeID> communities(3, 0);
  RegionCoarsening::respectPartition(phg, communities);
  EXPECT_NE(communities[0], communities[1]);
  EXPECT_EQ(communities[1], communities[2]);
  auto coarse_hg = hg.contract(communities, true);
  EXPECT_TRUE(coarse_hg.isFixed(communities[0]));
  EXPECT_EQ(coarse_hg.fixedVertexBlock(communities[0]), 0);
  EXPECT_EQ(coarse_hg.nodeWeight(communities[0]), 1);
}

} // namespace
} // namespace mt_kahypar
