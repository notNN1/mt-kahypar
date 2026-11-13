#include <sstream>

#include "gtest/gtest.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/partition/coarsening/coarsening_commons.h"
#include "mt-kahypar/utils/clustering_coefficient.h"

namespace mt_kahypar {
namespace {

TEST(CoarsestGraphStats, CliqueClusteringIncludesIsolatedVertices) {
  auto hg = ds::StaticHypergraphFactory::construct(5, 3,
    {{0, 1, 2}, {2, 3}, {0, 1, 2}}, nullptr, nullptr, true);
  // A triangle with a pendant vertex and an isolate: (1 + 1 + 1/3) / 5.
  // The repeated hyperedge must not change the simple clique expansion.
  EXPECT_DOUBLE_EQ(utils::meanCliqueExpansionClusteringCoefficient(hg), 7.0 / 15.0);
}

TEST(CoarsestGraphStats, CliqueClusteringIsZeroWithoutTriangles) {
  auto hg = ds::StaticHypergraphFactory::construct(4, 3,
    {{0, 1}, {1, 2}, {2, 3}}, nullptr, nullptr, true);
  EXPECT_DOUBLE_EQ(utils::meanCliqueExpansionClusteringCoefficient(hg), 0.0);
}

TEST(CoarsestGraphStats, RecordsContractedTopologyAndKeepsFirstStage) {
  auto hg = ds::StaticHypergraphFactory::construct(4, 3,
    {{0, 1}, {1, 2}, {2, 3}}, nullptr, nullptr, true);
  Context context;
  context.type = ContextType::main;
  context.partition.k = 2;
  context.partition.sp_process_output = true;
  context.coarsening.contraction_limit = 4;
  context.utility_id = utils::Utilities::instance().registerNewUtilityObjects();
  UncoarseningData<StaticHypergraphTypeTraits> data(false, hg, context);
  // Contract the middle pair. The coarsest graph is a three-vertex path.
  data.performMultilevelContraction({0, 1, 1, 2}, false,
    std::chrono::high_resolution_clock::now());
  data.finalizeCoarsening();

  std::ostringstream first;
  first << utils::Utilities::instance().getStats(context.utility_id);
  const std::string result = first.str();
  EXPECT_NE(result.find(" coarsening_num_nodes=3"), std::string::npos);
  EXPECT_NE(result.find(" coarsening_num_hyperedges=2"), std::string::npos);
  EXPECT_NE(result.find(" coarsening_mean_degree=1.33333"), std::string::npos);
  EXPECT_NE(result.find(" coarsening_mean_hyperedge_size=2"), std::string::npos);
  EXPECT_NE(result.find(" coarsening_min_hyperedge_size=2"), std::string::npos);
  EXPECT_NE(result.find(" coarsening_max_hyperedge_size=2"), std::string::npos);
  EXPECT_NE(result.find(" coarsening_num_articulation_points=1"), std::string::npos);
  EXPECT_NE(result.find(" coarsening_mean_clique_expansion_clustering_coefficient=0"), std::string::npos);

  // A later coarsening stage using the same utility object must keep the
  // entire first snapshot, even if its topology is different.
  auto triangle = ds::StaticHypergraphFactory::construct(3, 1,
    {{0, 1, 2}}, nullptr, nullptr, true);
  UncoarseningData<StaticHypergraphTypeTraits> later(false, triangle, context);
  later.finalizeCoarsening();
  std::ostringstream second;
  second << utils::Utilities::instance().getStats(context.utility_id);
  EXPECT_EQ(first.str(), second.str());
}

}  // namespace
}  // namespace mt_kahypar
