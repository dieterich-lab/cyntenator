
#include<string>
#include<iostream>
#include<list>
#include<vector>
#include "species_tree.h"
#include "genome.h"
#include "flow.h"
#include <stdlib.h>
using namespace std;


/*************************************************************************
 *************************************************************************
 ***                                                                   ***
 ***               1) Binary_Tree                                      ***
 ***                                                                   ***
 *************************************************************************
 *************************************************************************/

/*******************************************************************
 *
 * Constructor for a species tree given in the following format
 *
 * "(((maus:1 ratte:1.2):2 mensch:2.5):3 (Kuh:1.1 Pferd:1.5):2.2)"
 *
 *                       |
 *                       |
 *                   saeugetiere
 *                     /    \
 *                    /      \
 *                   3        2.2
 *                  /           \
 *             mensch            \
 *            nagetier            \
 *              /  \               \ 
 *             2    \               \
 *            /      \               \
 *         rodent     2.5          Paarhufer 
 *         /  \         \            /  \
 *        1   1.2        \        1.1    1.5 
 *       /       \        \       /        \
 *    Maus      Ratte   Mensch  Kuh       Pferd
 *
 *******************************************************************/ 
binary_tree::binary_tree(string tree ){


  vector<string> order = get_alignment_order( tree );   // split Spring into single alignments 
  sum_of_edges = 0;
  
                                           //*************************
  for( int i = 0; i != order.size(); i++){ // Insert Nodes and Edges
                                           //*************************
    string subalign = order[i];
     string anr = intToString( i+1 );
     string repl = "alignment" + anr;       
     string* alis = split_subtree( order[i] );  //split single alignments (human:1 alignment1):1
        
     node* v = get_node( alis[0] );   // u is the current node and father of
     node* w = get_node( alis[1] );   // v and w
     node* u = set_father( repl, v, w, get_branch_length( subalign ));  
     root = u;
     
     for( int j = i; j != order.size(); j++){  // replace subalign with repl in the rest of the list
      order[j] = replace( order[j], subalign, repl );
      }
     
  }
  
  sum_of_edges =  sum_of_subtree_edges(root->get_left_child()) +  sum_of_subtree_edges(root->get_right_child());
  
  vector<node*> leafs = get_leafs();
  for( int i = 0; i != leafs.size(); i++){
    for( int j = i+1; j < leafs.size(); j++){
      node* a = leafs[i];
      node* b = leafs[j];
      node* lacoan = lca( leafs[i] , leafs[j]);
      double dist_ab = distance_to_root( a ) +  distance_to_root( b )  - 2*distance_to_root( lacoan);
      // cerr << a->get_name() << "\t" << b->get_name() <<  " " << dist_ab/sum_of_edges  <<"\n";
      distances[ a->get_name() ][ b->get_name() ] = dist_ab / sum_of_edges;
      distances[ b->get_name() ][ a->get_name() ] = dist_ab / sum_of_edges;
      
    }
  }
}


/**
 * Function: sum_of_subtree_edges
 * Description: Sum up over all edges in the subtree induced by node u
 *  
 **/

double binary_tree::sum_of_subtree_edges( node* u ){

  if( u->get_left_child() == 0 || u->get_right_child() == 0)
    return u->get_distance();
  else{
    return  u->get_distance() + sum_of_subtree_edges( u->get_left_child()  ) + sum_of_subtree_edges( u->get_right_child());
  }
  
    
}


double binary_tree::distance_to_root( node* u ){
  if( u == root) return 0;
  else{
    return u->get_distance() + distance_to_root( u->get_father() );
  }

}

vector<node*> binary_tree::get_leafs(){
  vector<node*> leafs;
  for( map<string,node*>::const_iterator it = nodes.begin(); it != nodes.end(); it++ ){
    node* u = it->second;
    if( u->get_left_child() == 0)
      leafs.push_back( u );
  }

  return leafs;
}

/**
 * Determine the last common ancestor for a and b
*/

node* binary_tree::lca( node* a, node* b){
  vector<node*> ancestors_a;
  ancestors_a = get_ancestors( ancestors_a, a);
   vector<node*> ancestors_b;
  ancestors_b = get_ancestors( ancestors_b, b);
  for( int i = 0 ; i != ancestors_a.size(); i++)
    for( int j = 0 ; j != ancestors_b.size(); j++){
      if( ancestors_a[i] == ancestors_b[j] ){
	return ancestors_a[i];
      }
    }
  return 0;
}


/*
 * Search the pair of leafs with the largest
 * last common ancestor subtree size.
 *
 */

node*  binary_tree::lca( vector<string> alignment_species ){
  //cerr <<  alignment_species[0] << "\t"<<  nodes.size() <<"\n";
  
  /*for( map<string,node*>::iterator it =nodes.begin(); it != nodes.end(); it++){
    cerr << it->first << "\t" << (it->second)->get_name() << "\n";
    }*/
  
  node*  anc = nodes[ alignment_species[0] ];
  double subtree_length = 0; 
  for( int i = 0; i < alignment_species.size(); i++  )
    for( int j = i+1; j < alignment_species.size(); j++  ){
      node* child1 = nodes[ alignment_species[i] ];
      node* child2 = nodes[ alignment_species[j] ];
      node*  p = lca( child1, child2);
      double s = sum_of_subtree_edges( p );
      //cerr << child1->get_name() << "\t" << child2->get_name() << "\t" << s << "\n";
      if( s > subtree_length ){     
	subtree_length = s;
	anc = p;
      }
	
      }
  //cerr << anc->get_name() << "\t" << subtree_length << "\n";
  return anc;
}

vector<node*> binary_tree::get_ancestors( vector<node*>& vec, node* u){
  vec.push_back(u);
  if( u == root) return vec;
  else{
    return get_ancestors( vec, u->get_father() );
  }
}


node* binary_tree::set_father( string s, node* v, node* w, double d){
  node* u;
  u = new node();
  u->set_name( s );  // father 
  u->set_distance( d );
  u->set_left_child(v);
  u->set_right_child(w);
  v->set_father( u );
  w->set_father( u );
  nodes[ s ]= u;
  return u;
}


/**
 * 
 * get node from the nodes map<strin,node*> or add a new
 * node corresponding to the name of the subalign
 *
 */

node* binary_tree::get_node( string subalign ){
  node * v;
  string name = get_species_name( subalign );

  if( nodes.find( name )== nodes.end() ){  
       v = new node();
       double dist = get_branch_length( subalign );       
       v->set_name( name );   
       v->set_distance( dist );
       nodes[ name ] = v;
       
     }
  else
    v = nodes[ name ];
  //cerr << get_branch_length( subalign );
  return v;

}


/*************************************************************************
 *************************************************************************
 ***                                                                   ***
 ***               2) Global Functions                                 ***
 ***                                                                   ***
 *************************************************************************
 *************************************************************************/



/**
 * Parser:
 * Function: get_alignment_order
 *
 * Description: Parses a complete string like 
 * 
 *              "(((maus:1 ratte:1):1 mensch:2):1 (Kuh:1 Pferd:1):3)"
 *
 *               and splits into distinct alignments 
 **/


vector<string> get_alignment_order( string tree ){  // EXTENSION WITH SCORES
list<int> positions;
  vector<string> order;
  //nodes = map<string,node>();
  
  for( int i = 0; i < tree.size(); i++ ){  // parse tree and store matching bracket pairs in order
    char c = tree[i];                      // (...):score is read by parsing to next not_paranthesis
    if( c == '(' )
      positions.push_back( i );
    
    if( c == ')' ){
      i++;
      while( i <=  tree.size() &&  !is_space(tree[i]) && !is_paranthese(tree[i])  )
	i++;
      i--;
      int last = getLast( positions);
      string s = tree.substr(last, i -last +1);
      order.push_back(s);
    }
  }
  return order;
}





/**
 * Returns and removes the last position of an opening
 * paranthesis in the tree
 *
 **/

int getLast( list<int>& positions ){

  list<int>::const_iterator it = positions.end(); it--;
  int last = *it;
  positions.pop_back();
  return last;

}

/**
 * Parser:
 * Splits a subanlignment into 2 components
 *
 **/

string* split_subtree( string subtree){
  string* subaligns = new string[2]; 
  int count = 0;
  int pos[4];
  int j = 0;
  for( int i = 0; i < subtree.size(); i++  ){
    if( is_space( subtree[i]) || is_paranthese( subtree[i] )  ){  // SEARCH FOR FIRST SPACE
      pos[j++] = i+1;
    }
  }
  
  subaligns[0] = subtree.substr( pos[0], pos[1]-pos[0]-1);   // EXTRACT POSITIONS
  subaligns[1] = subtree.substr( pos[1], pos[2]-pos[1]-1);
  return subaligns;
}


/**
 * Parser:
 * Returns the branchlenth from a string like (rat.txt:0.04  mouse.txt:0.05):1.000
 * Iterates backward over the string until a : is found and returns the second half
 */

double get_branch_length( string s){

  int i ;
  for( i =  s.size()-1; i >= 0; i-- ){
    if( s[i] == ':'){
      break;
    }
    
 }
  if( i < s.size() && i >= 0 ){
    string tmp =  s.substr(i+1,s.size()-i);
    return  atof(tmp.c_str());
  }
  else return 1.0;
  }


/**
 * Parser:
 * Returns the species name from a string like mouse.txt:1.000 or alignment2
 *
 */

string get_species_name( string s){
   int i ;
   for( i = 0; i < s.size(); i++ ){
     if( s[i] == ':'){
       break;
     } 
   }
   if( i < s.size() ){
     string tmp =  s.substr(0,i);
     return  tmp;
   }
   else return s;

  }
