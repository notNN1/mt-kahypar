#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace mt_kahypar::utils {

// Mean local clustering coefficient of the simple, unweighted clique expansion.
// Each hyperedge induces a clique. Parallel edges and self-neighbors are ignored;
// vertices with fewer than two neighbors contribute zero to the mean.
template<typename Hypergraph>
double meanCliqueExpansionClusteringCoefficient(const Hypergraph& hg) {
  std::vector<std::vector<size_t>> neighbors(hg.initialNumNodes());
  size_t num_nodes = 0;
  for (const auto u : hg.nodes()) {
    ++num_nodes;
    auto& adjacent = neighbors[u];
    for (const auto e : hg.incidentEdges(u)) {
      for (const auto v : hg.pins(e)) {
        if (u != v) {
          adjacent.push_back(v);
        }
      }
    }
    std::sort(adjacent.begin(), adjacent.end());
    adjacent.erase(std::unique(adjacent.begin(), adjacent.end()), adjacent.end());
  }

  // Enumerate each triangle u < v < w once and credit all three vertices.
  std::vector<uint64_t> triangles(hg.initialNumNodes(), 0);
  for (const auto u : hg.nodes()) {
    const auto& adjacent = neighbors[u];
    for (auto v = std::upper_bound(adjacent.begin(), adjacent.end(), u);
         v != adjacent.end(); ++v) {
      auto a = v + 1;
      const auto& other = neighbors[*v];
      auto b = std::upper_bound(other.begin(), other.end(), *v);
      while (a != adjacent.end() && b != other.end()) {
        if (*a < *b) {
          ++a;
        } else if (*b < *a) {
          ++b;
        } else {
          ++triangles[u];
          ++triangles[*v];
          ++triangles[*a];
          ++a;
          ++b;
        }
      }
    }
  }

  double sum = 0.0;
  for (const auto u : hg.nodes()) {
    const size_t degree = neighbors[u].size();
    if (degree > 1) {
      sum += 2.0 * triangles[u] / (static_cast<double>(degree) * (degree - 1));
    }
  }
  return num_nodes > 0 ? sum / num_nodes : 0.0;
}

}  // namespace mt_kahypar::utils
