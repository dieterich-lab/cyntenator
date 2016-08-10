#ifndef GUARD_flow_h
#define GUARD_flow_h

#include<vector>
#include<string>
#include<list>
#include "localign.h"
#include "genome.h"
#include "species_tree.h"

using namespace std;

class program_flow{

  scoring_scheme scorer;
  string tree;
  string modus;
  int filter;   // if 0 only unique alignments
  int max_coverage; 
  int min_length;  
  bool print_all;
  string output;

  map<string, alignment_set> Alis;
  //map

  string align( string subtree , int i);
  void print_help();
  int get_parameters( int argc,  char** argv );

 public:
  program_flow(int argc,  char** argv);
  

};








string printHeader(  alignment_set A );
//void printHeader( string s );
string replace( string s,string search_string,  string with );




#endif
