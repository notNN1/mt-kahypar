set(ConnectedComponentsSources
        compute_components.cpp
        connectivity_facade.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${ConnectedComponentsSources})

target_sources(MtKaHyPar-ToolsSources INTERFACE ${ConnectedComponentsSources})
