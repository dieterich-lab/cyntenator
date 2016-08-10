#ifndef GUARD_species_tree_h
#define GUARD_species_tree_h
#include<string>
#include<map>

using namespace std;


class node{
  string name;
  node* father;
  node* lc;
  node* rc;
  double distance_to_dad;
 public:
  node(){ name =""; lc = 0; rc = 0; father = 0; distance_to_dad = 1;}
  void set_name( string s ){ name = s; };
  void set_father( node* v){ father = v;   }
  void set_left_child( node* v){ lc = v;   }
  void set_right_child( node* v){ rc = v;  }
  void set_distance(double d){ distance_to_dad =  d;}
  double  get_distance(){ return distance_to_dad; }
  node* get_left_child(){ return lc;}
  node* get_right_child(){ return rc;}
  node* get_father(){ return father;}
  string get_name(){return name;}
 
};

class binary_tree{

  node* root;                                  // pointer to the root node
  map<string,node*> nodes;                     // Mapping of species and inner node names to nodes
  map<string, map<string,double> > distances;  // Matrix of Sum of Edge distance between two nodes
  node* get_node( string s);                   // Fetch s from map "nodes" or add a new node 
  node* set_father( string s, node* v, node* w, double d);   //
  double sum_of_edges;
  
 public:
  binary_tree(){}
  binary_tree( string s );
  node* getRoot(){ return root;}
  node* lca( node* a, node* b);
  node* lca( vector<string> alignment_species );
  double distance_to_root( node* u);
  vector<node*> get_ancestors( vector<node*>& vec, node* u);
  double sum_of_subtree_edges( node* u );
  double get_tree_size( ){ return sum_of_subtree_edges( root); }
  double get_subtree_size( string a, string b){ return distances[a][b];}  // relative subtree size (a,b,(lca(a,b)))

  vector<node*> get_leafs();
};




int getLast( list<int>& positions );
vector<string> get_alignment_order( string tree);
string* split_subtree( string subtree );
double get_branch_length( string s);
string get_species_name( string s);
#endif
