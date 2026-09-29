set(UtilSources
      timer.cpp
      delete.cpp
      memory_tree.cpp
      connected_component_stats.cpp
    )

target_sources(MtKaHyPar-Sources INTERFACE ${UtilSources})
target_sources(MtKaHyPar-ToolsSources INTERFACE ${UtilSources})
