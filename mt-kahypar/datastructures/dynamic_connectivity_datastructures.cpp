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
#include <signal.h>

namespace mt_kahypar {
namespace ds {

    template<typename PartitionedHypergraph>    
    bool BFSConnectivity<PartitionedHypergraph>::moveVertex(
        const PartitionedHypergraph& phg, 
        const HypernodeID& hn
    ) {

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

        // calculate lost nodes
        size_t count = 0;
        size_t count1 = 0;
        for (const HypernodeID& hn : phg.nodes()) {
            if (has_connection_to_other_partition.isSet((size_t) hn)) {
                count1++;
                if (this->connected_to[hn].size() > 1) {
                    count++;
                }
            }
        }

        LOG << "Number of nodes that are not leaves and connected to other partitions: " << count;
        LOG << "Number of nodes that are: " << count1 - count;
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
        const HypernodeID& hn
    ) {
        if (this->vertex_to_parent_compressed[hn] == hn) {
            return false;   
        }
        return this->connected_to[hn].size() == 1;
    };

    template<typename PartitionedHypergraph>    
    void SpanningTreeConnectivity<PartitionedHypergraph>::moveVertex(
        const PartitionedHypergraph& phg,
        const HypernodeID& hn,
        const PartitionID& to
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

    template<typename PartitionedHypergraph>
    BFSSpanningTreeConnectivity<PartitionedHypergraph>::BFSSpanningTreeConnectivity(const PartitionedHypergraph& phg, const Context& context) {
        reset(phg, context);
    }

    template<typename PartitionedHypergraph>
    void BFSSpanningTreeConnectivity<PartitionedHypergraph>::reset(
        const PartitionedHypergraph& phg, 
        const Context& context
    ) {
        LOG << "initialized on Hypergraph instance: " << static_cast<const void*>(&phg) << " with " << phg.initialNumNodes() << " nodes";
        (void)context;

        this->hn_to_num_children.resize(phg.initialNumNodes());
        this->hn_to_parent.resize(phg.initialNumNodes());
        this->has_connection_to_other_partition.resize(phg.initialNumNodes());

        for (const HypernodeID& node : phg.nodes()) {
            this->hn_to_parent[node] = node;
        }

        // find nodes that have connections to other partitions
        Bitset edge_colored;
        edge_colored.resize(phg.initialNumEdges());

        for (const HypernodeID& hn : phg.nodes()) {            
            PartitionID current_partition = phg.partID(hn);

            for (const HyperedgeID& he : phg.incidentEdges(hn)) {

                if (edge_colored.isSet((size_t) he)) {
                    continue;
                }

                edge_colored.set((size_t) he);

                for (const HypernodeID& incident_hn : phg.pins(he)) {

                    if (current_partition != phg.partID(incident_hn)) {
                        this->has_connection_to_other_partition.set((size_t) hn);
                        break;
                    }
                }

                if (this->has_connection_to_other_partition.isSet((size_t) hn)) {
                    // now all the nodes have connections to different partitions as well

                    for (const HypernodeID& incident_hn : phg.pins(he)) {
                        this->has_connection_to_other_partition.set((size_t) incident_hn);
                    }
                }
            }
        }

        // do BFS with two queues

        Bitset node_colored;
        node_colored.resize(phg.initialNumNodes());
        
        edge_colored.reset();

        std::queue<HypernodeID> node_queue;
        std::queue<HypernodeID> colored_node_queue;
        PartitionID current_partition;

        for (const HypernodeID& hn : phg.nodes()) {
            if (node_colored.isSet((size_t) hn)) {
                continue;
            }

            current_partition = phg.partID(hn);

            if (this->has_connection_to_other_partition.isSet((size_t) hn)) {
                colored_node_queue.push(hn);
            }
            else {
                node_queue.push(hn);
            }

            node_colored.set((size_t) hn);
            
            edge_colored.reset();

            while (node_queue.size() > 0 || colored_node_queue.size() > 0) {
                
                HypernodeID current;

                if (node_queue.size() > 0) {
                    current = node_queue.front();
                    node_queue.pop();
                }
                else {
                    current = colored_node_queue.front();
                    colored_node_queue.pop();
                }


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

                        if (incident_hn == current) {
                           continue;
                        }

                        node_colored.set((size_t) incident_hn);

                        if (this->has_connection_to_other_partition.isSet((size_t) incident_hn)) {
                            colored_node_queue.push(incident_hn);
                        }
                        else {
                            node_queue.push(incident_hn);
                        }

                        this->hn_to_num_children[current]++;
                        this->hn_to_parent[incident_hn] = current;
                    }
                }
            }   
        }
        int count1 = 0;
        for (const HypernodeID& node : phg.nodes()) {
            if (this->hn_to_parent[node] == node) {
                count1++;
            }
        }

        int count2 = ConnectivityFacade<PartitionedHypergraph>::extra_component_count(phg, context);

        if (count1 != count2) {
            while (true)
            {
               LOG << "Too many roots: " << count1 << "roots and " << count2 << " components";
            }
            
        }
        is_tree_valid(phg);
    }

    template<typename PartitionedHypergraph>
    bool BFSSpanningTreeConnectivity<PartitionedHypergraph>::canMoveVertex(const PartitionedHypergraph& phg, const Context& context, const HypernodeID& hn) {
        (void) context;
        //LOG << "The size is wron: " << this->hn_to_num_children.size();
        //LOG << "HN: " << hn;
        assert(hn < this->hn_to_num_children.size());

        if (this->hn_to_parent[hn] == hn) {
            return false;
        }

        return this->hn_to_num_children[hn] == 0;
    }


    template<typename PartitionedHypergraph>
    void BFSSpanningTreeConnectivity<PartitionedHypergraph>::moveVertex(
        const PartitionedHypergraph& phg,
        const Context& context,
        const HypernodeID& hn,
        const PartitionID& to
    ) {
        vec<std::string> res = is_tree_valid(phg);
        if (res.size() != 0) {

            LOG << "Beginning of moveVertex";
            for (const std::string& err : res) {
                LOG << err;
            }

            raise(SIGTRAP); 
        }
        
        if (!canMoveVertex(phg, context, hn)) {
            LOG << "Move node without canMoveVertex being true";
            raise(SIGTRAP); 
        }

        if (phg.partID(hn) == to) {

            LOG << "Node is already in the selected partition";

            raise(SIGTRAP); 
        }


        // find candidate to attach to 
        for (const HyperedgeID& he : phg.incidentEdges(hn)) {
            for (const HypernodeID& incident_hn : phg.pins(he)) {
                if (phg.partID(incident_hn) != to) {
                    continue;
                }

                bool can_move_to_colored_node = this->has_connection_to_other_partition.isSet((size_t) incident_hn) && this->hn_to_num_children[incident_hn] > 0;
                bool can_move_to_regular_node = !this->has_connection_to_other_partition.isSet((size_t) incident_hn);
                
                if (can_move_to_colored_node || can_move_to_regular_node) {
                    HypernodeID parent = this->hn_to_parent[hn];
                    this->hn_to_num_children[parent]--;

                    this->hn_to_parent[hn] = incident_hn;
                    this->hn_to_num_children[incident_hn]++;

                    LOG << "special attachment";

                    return;
                }
            }
        }

        // attach to the next best node
        for (const HyperedgeID& he : phg.incidentEdges(hn)) {
            for (const HypernodeID& incident_hn : phg.pins(he)) {
                if (phg.partID(incident_hn) != to) {
                    continue;
                }

                HypernodeID parent = this->hn_to_parent[hn];
                this->hn_to_num_children[parent]--;

                this->hn_to_parent[hn] = incident_hn;
                this->hn_to_num_children[incident_hn]++;

                LOG << "default attachment";
                return;
            }
        }

        LOG << "There has not been found a node to attach to for node " << hn;
        while(true);
    }

    template<typename PartitionedHypergraph>
    void ConnectivityFacade<PartitionedHypergraph>::initialize_can_move_current_node_to_partition(
        const PartitionedHypergraph& hypergraph,
        const DynamicConnectivityStrategy& strategy,
        const HypernodeID& hn
    ) {
        this->can_move_current_node_to_partition.reset();
        this->can_move_current_node_to_partition.resize(hypergraph.k());

        if (
        strategy == DynamicConnectivityStrategy::bfs
        || strategy == DynamicConnectivityStrategy::h_vertex_degree
        || strategy == DynamicConnectivityStrategy::st
        ) {
            for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
                for (const PartitionID& partition : hypergraph.connectivitySet(he)) {
                    this->can_move_current_node_to_partition.set((size_t) partition); 
                }
            }

            this->current_node = hn;
        }
    }

    template<typename PartitionedHypergraph>
    bool ConnectivityFacade<PartitionedHypergraph>::can_move_into_partition(
        const DynamicConnectivityStrategy& strategy,
        const PartitionedHypergraph& hypergraph,
        const Context& _context,
        const HypernodeID& hn,
        const PartitionID& to
    ) {
        if (
            strategy == DynamicConnectivityStrategy::bfs
            || strategy == DynamicConnectivityStrategy::h_vertex_degree
            || strategy == DynamicConnectivityStrategy::st
        ) {
            if (hn != this->current_node) {
                this->initialize_can_move_current_node_to_partition(hypergraph, strategy, hn);
            }

            return this->can_move_current_node_to_partition.isSet((size_t) to);
        }
        return true; // default
    }

    template<typename PartitionedHypergraph>
    bool ConnectivityFacade<PartitionedHypergraph>::can_move_out_of_partition(
        const DynamicConnectivityStrategy& strategy,
        const PartitionedHypergraph& hypergraph,
        const Context& _context,
        const HypernodeID& hn
    ) {
        bool can_move_node = true;

        if (strategy == DynamicConnectivityStrategy::bfs) {
            mt_kahypar::ds::BFSConnectivity<PartitionedHypergraph> dcd = mt_kahypar::ds::BFSConnectivity<PartitionedHypergraph>();
            can_move_node = dcd.moveVertex(hypergraph, hn);
        }
        else if (strategy == DynamicConnectivityStrategy::h_vertex_degree) {
            auto range = hypergraph.incidentEdges(hn);
            can_move_node = std::distance(range.begin(), range.end()) > 4;
        }
        else if (strategy == DynamicConnectivityStrategy::st) {
            can_move_node = hypergraph.canMoveVertex(strategy, _context, hn);
        }

        return can_move_node;
    }

    template<typename PartitionedHypergraph>
    void ConnectivityFacade<PartitionedHypergraph>::move_out_of_partition(
        const DynamicConnectivityStrategy& strategy,
        const PartitionedHypergraph& hypergraph,
        const HypernodeID& hn,
        const PartitionID& to
    ) {
        if (strategy == DynamicConnectivityStrategy::st) {
            hypergraph.moveVertex(strategy, hn, to);
        }

        this->current_node = -1;
    }



INSTANTIATE_CLASS_WITH_PARTITIONED_HG(BFSConnectivity)
INSTANTIATE_CLASS_WITH_PARTITIONED_HG(SpanningTreeConnectivity)
INSTANTIATE_CLASS_WITH_PARTITIONED_HG(BFSSpanningTreeConnectivity)
INSTANTIATE_CLASS_WITH_PARTITIONED_HG(ConnectivityFacade)

}  // namespace ds
}  // namespace mt_kahypar
