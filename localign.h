#ifndef GUARD_localign_h
#define GUARD_localign_h
#include <vector>
#include <map>
#include <string>
#include "genome.h"
#include "species_tree.h"
#include "string_hash.h"

using namespace std;



class scoring_scheme{
 private:
  double miss;
  double gap;
  double match;
  double threshold;
  double base_missmatch;
  double base_gap;
  double base_threshold;
  double base_match;
 
  string mode;
  map<string, string> species_of;   // this is used only in phylo modus
  binary_tree species_tree;
  
  //shmap< shmap<double>::map >::map scoring_table;  // Hash_map version
  map<string, map<string, int> > scoring_table;
 
  list<string> get_all_homologs( string query );
  void read_single_species( string species , string file);

 public:

  // constructors
  scoring_scheme();
  scoring_scheme( string modus, const char* file);
  scoring_scheme( string modus, const char* file, string alignment_tree, string spec_tree );

  void read_species( string );
  double scoreID( const column& a, const column& b);
  double scoreBLAST( const column& a, const column& b );
  double score( const column& a, const column& b);
  double scoreGenes( const location & loc_a, const location& loc_b);

  // basic get and set functions  
  string get_mode(){ return mode;}
  void set_miss( double x ){  base_missmatch = x;}
  double get_miss(){ return miss; }
  void set_gap( double x){ base_gap = x;}
  double get_gap(){ return gap; }
  void set_threshold( double x){ base_threshold = x;}
  double get_threshold(){ return threshold; }

  // species tree for phylogenetic progressive alignment
  void set_species_tree( string s ){ species_tree=binary_tree( s );}
  void phylogentic_parameter_adaptation( const alignment& X, const alignment& Y);
  
  // ** THIS IS THE LINE WHERE I NEED THE SPECIES OF GENE**
  double distance ( const string&  gid , const string& gjd ){
        double dist = species_tree.get_subtree_size( gid, gjd );
	//double dist = species_tree.get_subtree_size( species_of[gid], species_of[gjd] );
    //cerr << "dist: " << gid << "\t" << gjd << "\n";
    return dist;

  }
  
  list<string> get_all_homologs( alignment A );
  
};



class local_aligner{
 
 private:
 alignment A;
 alignment B;
 vector<vector<char> > T;
 vector<vector<double> > S;
 bool adapt_scores;   // whether to adapt parameters in every step


 public: 
 local_aligner(){   adapt_scores=1; }
 void parameter_fixed(){ adapt_scores = 0; }
 list<column> trace_back(int i, int j , list<column>& newAli, const double& min , const double& max);
 double last_minimum( int i, int j , const double& final_score, double min);
 void smith_waterman( const alignment& A, const alignment& B, alignment_set& alignments, scoring_scheme& scorer);

};


class peak{
  int i;
  int j;
  double score;
 public:
  peak( const int i, const int j, const double score);
  int get_i()     const { return i;}
  int get_j()     const { return j;}
  double get_score() const { return score; }



};



vector< vector<char> > trac_back_matric( int n, int m );
vector< vector<double> > scoring_matrix( int n, int m );
bool comparePeaks( const peak& a, const peak& b);

#endif

