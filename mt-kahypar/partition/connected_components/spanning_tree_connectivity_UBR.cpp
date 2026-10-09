
/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/


#include "mt-kahypar/partition/connected_components/spanning_tree_connectivity_UBR.h"
#include "mt-kahypar/utils/connected_component_stats.h"
#include <queue>
#include <cassert>
#include "mt-kahypar/definitions.h"
#include <signal.h>

namespace mt_kahypar {
namespace connected_components {

    template<typename PartitionedHypergraph>   
    SpanningTreeConnectivity<PartitionedHypergraph>::SpanningTreeConnectivity(
        const PartitionedHypergraph& phg,
        const Context& context
    ) {

        Bitset edge_already_seen;
        edge_already_seen.resize(phg.initialNumEdges());

        // find nodes that have connections to other partitions
        Bitset has_connection_to_other_partition;
        has_connection_to_other_partition.resize(phg.initialNumNodes());

        for (const HypernodeID& hn : phg.nodes()) {            
        PartitionID current_partition = phg.partID(hn);

            for (const HyperedgeID& he : phg.incidentEdges(hn)) {

                if (edge_already_seen.isSet((size_t) he)) {
                    continue;
                }

                edge_already_seen.set((size_t) he);

                for (const HypernodeID& incident_hn : phg.pins(he)) {

                    if (current_partition != phg.partID(incident_hn)) {
                        has_connection_to_other_partition.set((size_t) hn);
                        break;
                    }
                }

                if (has_connection_to_other_partition.isSet((size_t) hn)) {
                    // now all the nodes have connections to different partitions as well

                    for (const HypernodeID& incident_hn : phg.pins(he)) {
                        has_connection_to_other_partition.set((size_t) incident_hn);
                    }
                }
            }
        }

        // There are nodes, that cannot be avoided when building the spanning tree, which are also nodes, that connect to other partitions
        // Then a new run for the connected component tries to use these nodes first, to build the spanning tree

        vec<vec<ConnectedComponent>> connected_components;
        connected_components::compute_components_per_block(
            phg,
            context,
            connected_components
        );


        this->connected_to.resize(phg.initialNumNodes());
        this->vertex_to_parent_compressed.resize(phg.initialNumNodes());
        this->vertex_to_rank.resize(phg.initialNumNodes());

        for (HypernodeID hypernode : phg.nodes()) {
            this->vertex_to_parent_compressed[hypernode] = hypernode;
        }

        vec<Bitset> colored_edges_per_partition;
        colored_edges_per_partition.resize(phg.k());

        for (size_t i = 0; i < phg.k(); i++) {
            colored_edges_per_partition[i].resize(phg.initialNumEdges());
        }

        PartitionID current_partition;


        for (const HypernodeID& hn : phg.nodes()) {

            if (has_connection_to_other_partition.isSet((size_t) hn)) {
                continue;
            }

            current_partition = phg.partID(hn);

            for (       const HyperedgeID& he          : phg.incidentEdges(hn)          ) {

                if (current_partition != kInvalidPartition) {
                    if (colored_edges_per_partition[current_partition].isSet((size_t) he)) {
                        continue;
                    }

                    colored_edges_per_partition[current_partition].set((size_t) he);
                }

                for (   const HypernodeID& incident_hn : phg.pins(he)                   ) {

                    if (!has_connection_to_other_partition.isSet((size_t) incident_hn) 
                        ||(has_connection_to_other_partition.isSet((size_t) incident_hn) && this->connected_to[incident_hn].size() == 0)
                    ) {
                        if (!is_same_component(hn, incident_hn)) {
                            connect_nodes(hn, incident_hn);
                        }
                    }

                }
            }
        }

        for (size_t k = 0; k < phg.k(); k++) {
            colored_edges_per_partition[k].reset();
        }

        // now only connect nodes with 'has_connection_to_other_partition'
        for (const HypernodeID& hn : phg.nodes()) {

            if (!has_connection_to_other_partition.isSet((size_t) hn)) {
                continue;
            }

            current_partition = phg.partID(hn);


            for (       const HyperedgeID& he          : phg.incidentEdges(hn)          ) {

                if (current_partition != kInvalidPartition) {
                    if (colored_edges_per_partition[current_partition].isSet((size_t) he)) {
                        continue;
                    }

                    colored_edges_per_partition[current_partition].set((size_t) he);
                }
                
                for (   const HypernodeID& incident_hn : phg.pins(he)                   ) {

                    
                    if (phg.partID(incident_hn) != current_partition) {
                        continue;
                    }
       
                    if (!is_same_component(hn, incident_hn)) {
                        connect_nodes(hn, incident_hn);
                    }
                }
            }
        }

        // calculate lost nodes
        size_t external_nodes           = 0;
        size_t blocked_external_nodes   = 0;

        for (const HypernodeID& hn : phg.nodes()) {
            if (has_connection_to_other_partition.isSet((size_t) hn)) {
                external_nodes++;
                if (this->connected_to[hn].size() > 1) {
                    blocked_external_nodes++;
                }
            }
        }

        mt_kahypar::utils::cc_debug.add_st_rebuild_stats(blocked_external_nodes, external_nodes);
    };


    template<typename PartitionedHypergraph>  
    int SpanningTreeConnectivity<PartitionedHypergraph>::calc_parents(const PartitionedHypergraph& phg) {
      int count = 0;
        for (const HypernodeID& hn : phg.nodes()) {
            if (this->vertex_to_parent_compressed[hn] == hn) {
                count++;
            }
        }

        return count;  
    }

    // is_same_component has to be called before, so hn and incident_hn are directly under the parent in the vertex_to_parent tree
    template<typename PartitionedHypergraph>    
    inline void  SpanningTreeConnectivity<PartitionedHypergraph>::connect_nodes(
        const HypernodeID& hn,
        const HypernodeID& incident_hn
    ) {
        HypernodeID pn          = this->vertex_to_parent_compressed[hn];
        HypernodeID incident_pn = this->vertex_to_parent_compressed[incident_hn];

        if (this->vertex_to_rank[pn] > this->vertex_to_rank[incident_pn]) {
            this->vertex_to_parent_compressed[incident_pn] = this->vertex_to_parent_compressed[pn];
        }
        else if (this->vertex_to_rank[pn] == this->vertex_to_rank[incident_pn]) {
            this->vertex_to_rank[pn]++;
            this->vertex_to_parent_compressed[incident_pn] = this->vertex_to_parent_compressed[pn];
        }
        else {
            this->vertex_to_parent_compressed[pn] = this->vertex_to_parent_compressed[incident_pn];
        }

        this->connected_to[hn].emplace_back(incident_hn);
        auto it1 = std::prev(this->connected_to[hn].end());

        this->connected_to[incident_hn].emplace_back(hn);
        auto it2 = std::prev(this->connected_to[incident_hn].end());

        it1->iterator = it2;
        it2->iterator = it1;
    }

    template<typename PartitionedHypergraph>    
    inline bool SpanningTreeConnectivity<PartitionedHypergraph>::is_same_component(
        const HypernodeID& hn1,
        const HypernodeID& hn2
    ) {
        HypernodeID parent1 = hn1;
        HypernodeID parent2 = hn2;

        vec<HypernodeID> path1;
        vec<HypernodeID> path2;


        while (this->vertex_to_parent_compressed[parent1] != parent1) {
            path1.push_back(parent1);
            parent1 = this->vertex_to_parent_compressed[parent1];
        }

        for (HypernodeID hp1 : path1) {
            this->vertex_to_parent_compressed[hp1] = parent1;
        }

        while (this->vertex_to_parent_compressed[parent2] != parent2) {
            path2.push_back(parent2);
            parent2 = this->vertex_to_parent_compressed[parent2];
        }

        for (HypernodeID hp2 : path2) {
            this->vertex_to_parent_compressed[hp2] = parent2;
        }

        return parent1 == parent2;        
    };

    template<typename PartitionedHypergraph>    
    bool SpanningTreeConnectivity<PartitionedHypergraph>::canMoveVertex(
        const Context& context,
        const HypernodeID& hn
    ) {
        if (this->vertex_to_parent_compressed[hn] == hn && this->connected_to[hn].size() > 0) {
            return false;   
        }
        return this->connected_to[hn].size() <= 1;
    };

    template<typename PartitionedHypergraph>    
    HypernodeID SpanningTreeConnectivity<PartitionedHypergraph>::moveVertex(
        const PartitionedHypergraph& phg,
        const HypernodeID& hn,
        const PartitionID& to,
        const HypernodeID& node_to
    ) {
        assert(this->connected_to[hn].size() == 1);
        assert(this->vertex_to_parent_compressed[hn] != hn);

        if (node_to == kInvalidHypernode) {
            LOG << "There has not been found a node to attach to for node " << hn;
            raise(SIGSEGV);
        }

        auto it1                = this->connected_to[hn].begin();
        auto it2                = it1->iterator;
        HypernodeID incident_hn = it1->node;

        // erase parent from child
        this->connected_to[incident_hn].erase(it2);
        this->connected_to[hn].erase(it1);

        // remove hn from component tree
        this->vertex_to_parent_compressed[hn] = hn;

        /*for (       const HyperedgeID& he           : phg.incidentEdges(hn) ) {
            for (   const HypernodeID& incident_hn  : phg.pins(he)          ) {
                if (phg.partID(incident_hn) == to && !is_same_component(hn, incident_hn)) {
                    connect_nodes(hn, incident_hn);
                    return;
                }
            }
        }*/

        connect_nodes(hn, node_to);

        if (this->connected_to[incident_hn].size() == 1) {
            return incident_hn;
        }

        return kInvalidHypernode;

    };

INSTANTIATE_CLASS_WITH_PARTITIONED_HG(SpanningTreeConnectivity)

}  // namespace connected_components
}  // namespace mt_kahypar
