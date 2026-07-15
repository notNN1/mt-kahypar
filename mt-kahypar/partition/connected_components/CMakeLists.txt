set(ConnectedComponentsSources
        compute_components.cpp
        connectivity_facade.cpp
        spanning_tree_connectivity.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${ConnectedComponentsSources})

target_sources(MtKaHyPar-ToolsSources INTERFACE ${ConnectedComponentsSources})
