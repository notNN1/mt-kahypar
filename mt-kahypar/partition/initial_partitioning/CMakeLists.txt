set(InitialPartitioningSources
        pool_initial_partitioner.cpp
        random_initial_partitioner.cpp
        bfs_initial_partitioner.cpp
        label_propagation_initial_partitioner.cpp
        st_initial_partitioner.cpp
        tarjan_initial_partitioner.cpp
        greedy_st_initial_partitioner.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${InitialPartitioningSources})
