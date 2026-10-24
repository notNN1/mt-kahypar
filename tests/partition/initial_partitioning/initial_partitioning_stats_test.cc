#include "gtest/gtest.h"

#include <sstream>
#include "mt-kahypar/utils/initial_partitioning_stats.h"

namespace mt_kahypar {
namespace {

TEST(InitialPartitioningStats, MergesPopulationVariancesAndSerializesThem) {
  const auto algorithm = static_cast<InitialPartitioningAlgorithm>(0);
  utils::InitialPartitionerSummary all(algorithm), left(algorithm), right(algorithm);
  for (size_t i = 0; i < 4; ++i) {
    utils::InitialPartitionerSummary one(algorithm);
    one.total_calls = 1;
    one.total_sum_quality = 1000000000.0 + 2 * i;
    one.total_sum_imbalance = 0.1 * i;
    one.total_extra_components = 2 * i;
    all.add(one);
    (i < 2 ? left : right).add(one);
  }
  left.add(right);
  EXPECT_DOUBLE_EQ(all.variance_quality(), 5.0);
  EXPECT_NEAR(all.variance_imbalance(), 0.0125, 1e-12);
  EXPECT_DOUBLE_EQ(all.variance_components(), 5.0);
  EXPECT_DOUBLE_EQ(left.variance_quality(), all.variance_quality());
  EXPECT_NEAR(left.variance_imbalance(), all.variance_imbalance(), 1e-12);
  EXPECT_DOUBLE_EQ(left.variance_components(), all.variance_components());

  std::ostringstream output;
  output << all;
  EXPECT_NE(output.str().find(" var_quality_"), std::string::npos);
  EXPECT_NE(output.str().find(" var_imbalance_"), std::string::npos);
  EXPECT_NE(output.str().find(" var_components_"), std::string::npos);
}

TEST(InitialPartitioningStats, EmptyAndSingleRunHaveZeroVariance) {
  const auto algorithm = static_cast<InitialPartitioningAlgorithm>(0);
  utils::InitialPartitionerSummary summary(algorithm), empty(algorithm);
  EXPECT_DOUBLE_EQ(summary.variance_quality(), 0.0);
  EXPECT_DOUBLE_EQ(summary.variance_imbalance(), 0.0);
  EXPECT_DOUBLE_EQ(summary.variance_components(), 0.0);
  EXPECT_DOUBLE_EQ(summary.percentage_best(0), 0.0);
  summary.total_calls = 1;
  summary.total_sum_quality = 42;
  summary.total_sum_imbalance = 0.02;
  summary.total_extra_components = 3;
  summary.add(empty);
  EXPECT_DOUBLE_EQ(summary.variance_quality(), 0.0);
  EXPECT_DOUBLE_EQ(summary.variance_imbalance(), 0.0);
  EXPECT_DOUBLE_EQ(summary.variance_components(), 0.0);
}

}  // namespace
}  // namespace mt_kahypar
