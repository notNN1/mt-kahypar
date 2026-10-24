#include "gtest/gtest.h"

#include "tests/definitions.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/partition/connected_components/tarjan.h"

namespace mt_kahypar {
namespace {

using PartitionedHypergraph = StaticHypergraphTypeTraits::PartitionedHypergraph;

TEST(TarjanReducedGraph, UpdatesArticulationPointsAfterEachRemoval) {
  auto hg = ds::StaticHypergraphFactory::construct(4, 4,
    {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, nullptr, nullptr, true);
  PartitionedHypergraph phg(2, hg);
  ds::Bitset excluded;
  excluded.resize(4);
  connected_components::Tarjan<PartitionedHypergraph> tarjan;

  EXPECT_FALSE(tarjan.is_articulation_point(phg, 2, excluded));
  excluded.set(0);
  EXPECT_TRUE(tarjan.is_articulation_point(phg, 2, excluded));
  EXPECT_FALSE(tarjan.is_articulation_point(phg, 1, excluded));
  excluded.set(1);
  EXPECT_FALSE(tarjan.is_articulation_point(phg, 2, excluded));
}

TEST(TarjanReducedGraph, ReconnectsRemainingPinsOfAHyperedge) {
  auto hg = ds::StaticHypergraphFactory::construct(4, 1,
    {{0, 1, 2, 3}}, nullptr, nullptr, true);
  PartitionedHypergraph phg(2, hg);
  ds::Bitset excluded;
  excluded.resize(4);
  excluded.set(0);
  connected_components::Tarjan<PartitionedHypergraph> tarjan;

  // Filtering the original cycle would make vertex 2 an articulation point.
  // The remaining three pins must instead form a new cycle.
  EXPECT_FALSE(tarjan.is_articulation_point(phg, 2, excluded));
  excluded.set(2);
  EXPECT_FALSE(tarjan.is_articulation_point(phg, 1, excluded));
  excluded.set(3);
  EXPECT_FALSE(tarjan.is_articulation_point(phg, 1, excluded));
}

}  // namespace
}  // namespace mt_kahypar
