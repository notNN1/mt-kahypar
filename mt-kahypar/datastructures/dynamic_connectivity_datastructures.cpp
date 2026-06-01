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

#include "mt-kahypar/datastructures/dynamic_connectivity_datastructures.h"
#include <tbb/task_group.h>
#include <queue>
#include <cassert>
#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace ds {

    template<typename PartitionedHypergraph>    
    bool BFSConnectivity<PartitionedHypergraph>::moveVertex(
        const PartitionedHypergraph& phg, 
        const Context& context,
        HypernodeID hn
    ) {
        (void) context;

        // compute components with hn set, so we can't move over hn
        int components = 0;

        Bitset node_colored;
        node_colored.resize(phg.initialNumNodes());
        node_colored.set((size_t) hn);

        Bitset edge_colored;
        edge_colored.resize(phg.initialNumEdges());
        
        std::queue<HypernodeID> node_queue;
        PartitionID current_partition = phg.partID(hn);

        for (const HypernodeID& hn_ : phg.nodes()) {
            if (node_colored.isSet((size_t) hn_) || phg.partID(hn_) != current_partition) {
                continue;
            }

            node_queue.push(hn_);
            node_colored.set((size_t) hn_);

            edge_colored.reset();

            while (node_queue.size() > 0) {
                HypernodeID current = node_queue.front();
                node_queue.pop();

                for (const HyperedgeID& he : phg.incidentEdges(current)) {

                    if (edge_colored.isSet((size_t) he)) {
                        continue;
                    }

                    edge_colored.set((size_t) he);

                    for (const HypernodeID& incident_hn : phg.pins(he)) {
                        if (node_colored.isSet((size_t) incident_hn)) {
                            continue;
                        }
                        
                        if (phg.partID(incident_hn) != current_partition) {
                            continue;
                        }

                        node_queue.push(incident_hn);
                        node_colored.set((size_t) incident_hn);
                    }
                }
            }
            
            components++;
            if (components > 1) {
                return false;
            }
        }
        return true;
    };


    template<typename PartitionedHypergraph>   
    SpanningTreeConnectivity<PartitionedHypergraph>::SpanningTreeConnectivity(
        const PartitionedHypergraph& phg,
        const Context& context
    ) {

        // find nodes that have connections to other partitions
        Bitset has_connection_to_other_partition;
        has_connection_to_other_partition.resize(phg.initialNumNodes());

        for (const HypernodeID& hn : phg.nodes()) {            
            PartitionID current_partition = phg.partID(hn);

            for (const HyperedgeID& he : phg.incidentEdges(hn)) {
                for (const HypernodeID& incident_hn : phg.pins(he)) {

                    if (current_partition != phg.partID(incident_hn)) {
                        has_connection_to_other_partition.set((size_t) hn);
                        break;
                    }
                }
            }
        }

        // There are nodes, that cannot be avoided when building the spanning tree, which are also nodes, that connect to other partitions
        // Save these in nodes_unavoidable
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

        Bitset edge_already_seen;
        edge_already_seen.resize(phg.initialNumEdges());

        for (HypernodeID hypernode : phg.nodes()) {
            this->vertex_to_parent_compressed[hypernode] = hypernode;
        }


        for (const HypernodeID& hn : phg.nodes()) {

            if (has_connection_to_other_partition.isSet((size_t) hn)) {
                continue;
            }

            for (       const HyperedgeID& he          : phg.incidentEdges(hn)          ) {

                if (edge_already_seen.isSet((size_t) he)) {
                    continue;
                }

                edge_already_seen.set((size_t) he);

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

        edge_already_seen.reset();

        PartitionID current_partition;
        // now only connect nodes with 'has_connection_to_other_partition'
        for (const HypernodeID& hn : phg.nodes()) {

            if (!has_connection_to_other_partition.isSet((size_t) hn)) {
                continue;
            }

            current_partition = phg.partID(hn);


            for (       const HyperedgeID& he          : phg.incidentEdges(hn)          ) {

                if (edge_already_seen.isSet((size_t) he)) {
                    //continue;
                }

                edge_already_seen.set((size_t) he);

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

    template<typename PartitionedHypergraph>    
    inline void  SpanningTreeConnectivity<PartitionedHypergraph>::try_connect_to_incident_without_connection(
        Bitset& has_connection_to_other_partition,
        const HypernodeID& hn,
        const HypernodeID& incident_hn
    ) {
        if (
            this->connected_to[incident_hn].size() == 0
        ) {
            if (!is_same_component(hn, incident_hn)) {
                connect_nodes(hn, incident_hn);
            }
        }
    };

    // connect to nodes without has_connection_to_other_partition
    template<typename PartitionedHypergraph>    
    inline void  SpanningTreeConnectivity<PartitionedHypergraph>::try_connect_to_incident_with_connection(
        Bitset& has_connection_to_other_partition,
        const HypernodeID& hn,
        const HypernodeID& incident_hn
    ) {
    
        if (!is_same_component(hn, incident_hn)) {
            connect_nodes(hn, incident_hn);
        }
        
    };

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
        HypernodeID hn
    ) {
        if (this->vertex_to_parent_compressed[hn] == hn) {
            return false;   
        }
        return this->connected_to[hn].size() == 1;
    };

    template<typename PartitionedHypergraph>    
    void SpanningTreeConnectivity<PartitionedHypergraph>::moveVertex(
        const PartitionedHypergraph& phg,
        const Context& context,
        HypernodeID hn,
        PartitionID to
    ) {
        assert(this->connected_to[hn].size() == 1);
        assert(this->vertex_to_parent_compressed[hn] != hn);


        auto it1                = this->connected_to[hn].begin();
        auto it2                = it1->iterator;
        HypernodeID incident_hn = it1->node;

        // erase parent from child
        this->connected_to[incident_hn].erase(it2);
        this->connected_to[hn].erase(it1);


        // remove hn from component tree
        this->vertex_to_parent_compressed[hn] = hn;

        for (       const HyperedgeID& he           : phg.incidentEdges(hn) ) {
            for (   const HypernodeID& incident_hn  : phg.pins(he)          ) {
                if (phg.partID(incident_hn) == to && !is_same_component(hn, incident_hn)) {
                    connect_nodes(hn, incident_hn);
                }
            }
        }

    };


INSTANTIATE_CLASS_WITH_PARTITIONED_HG(BFSConnectivity)
INSTANTIATE_CLASS_WITH_PARTITIONED_HG(SpanningTreeConnectivity)

}  // namespace ds
}  // namespace mt_kahypar
