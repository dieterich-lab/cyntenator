# Makefile for cyntenator
CC      = g++            # compiler
CPP	= cyntenator.cpp localign.cpp genome.cpp flow.cpp species_tree.cpp
HEAD	= localign.h genome.h flow.h species_tree.h

cyntenator: $(CPP) $(HEAD)
	$(CC) $(CPP) -o cyntenator
