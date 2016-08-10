project(Cyntenator)

cmake_minimum_required(VERSION 2.6.3 FATAL_ERROR)

set(sources cyntenator.cpp localign.cpp genome.cpp flow.cpp species_tree.cpp)
set(sources_h localign.h genome.h flow.h species_tree.h)

### source group definition
source_group("Header Files" FILES ${sources_h})

add_executable(Cyntenator ${sources})