#include <iostream>
#include <utility>
#include <map>
#include <cstddef>  // size_t
#include <stdlib.h>   // string to double atof
#include <list>
#include <string>
#include <fstream>
#include "genome.h"
#include "localign.h"
#include "flow.h"

using namespace std;


 /*         

            j in B
            | 
            V
         ...........
         ...........
A) i ->  ...........
         ...........
         ...........
  */


void local_aligner::smith_waterman( const alignment& X, const alignment& Y, alignment_set& alignments, scoring_scheme& scorer){
 
  A = X;
  B = Y;
  int n = A.size();
  int m = B.size();
  list<peak> local_maxima;
  if( adapt_scores )
    scorer.phylogentic_parameter_adaptation( X, Y ) ;
  
  S = scoring_matrix( n,m );
  T = trac_back_matric( n, m);
  
  for( int i = 0; i < n+1; i++ ) S[i][0] = 0;   // INITIALIZATION
   S[0] = vector<double>(m+1,0); 

 

  for( int i = 1; i < n+1; i++ ){
    for( int j = 1; j < m+1; j++ ){
      
      double a = S[i-1][j] + scorer.get_gap();                          //  UP
      double b = S[i][j-1] + scorer.get_gap();                          //  LEFT
      double c = 0;                                                     //  NEW
      double d = S[i-1][j-1] + scorer.score( A.get(i-1), B.get(j-1) );  // MATCH
      
      double tmp1 = max(a,b);
      double tmp2 = max(c,d);
      S[i][j] = max( tmp1, tmp2);
      
      if( S[i][j] == a ) T[i][j] = 'U'; // UP
      else if( S[i][j] == b ) T[i][j] = 'L'; // LEFT
      else if( S[i][j] == c ) T[i][j] = '.'; // START
      else T[i][j] = 'M'; // MATCH if( S[i][j] == d ) 
      
      if( S[i][j] > scorer.get_threshold() && T[i][j] == 'M'  ){
	local_maxima.push_back( peak( i, j, S[i][j]) );
      }
    }
  }

     

       /*********************************
       *
       * Optionally print Matrix
       *
       *********************************/
  
  /*for( int i = 0; i < n+1; i++ ){
        for( int j = 0; j < m+1; j++ )
	  cerr <<  S[i][j] <<"\t";
        cerr << "\n";
      }  
      */



     /****************************************
      *      Extract local alignments        
      * 
      * 1) list<peak> local_maxima
      *
      ***************************************/

    local_maxima.sort( comparePeaks);

    for( list<peak>::const_iterator it = local_maxima.begin(); it != local_maxima.end(); it++  ){  
         
      peak         p        = *it; 
      double       s        = p.get_score(); 
      list<column> pair_list; 
      double       last_min = last_minimum( p.get_i() , p.get_j() , s, s);
      
      if( s - last_min >= scorer.get_threshold() ){
	pair_list = trace_back( p.get_i() , p.get_j() , pair_list, last_min, s) ;
	
	if( pair_list.size() > 0 ){  
	  alignment newAlignment = alignment( pair_list, p.get_score() - last_min );	  
	  alignments.add( newAlignment );
	}
	
      }
    }
}


/**
 *
 * Return the last minimum, before the trace_back score exceeds
 * the final score;
 * 
 * example trace score:
 *
 *    p1   p2
 *    /\ 
 *   /  \  /\
 *  /    \/  \
 * /          \______
 *       LM := last minimum
 *
 **/

double local_aligner::last_minimum( int i, int j , const double& final_score, double min){


  //cerr << i << " " << j << " " << min << "\n";  
  char trace = T[i][j];
  double current_score = S[i][j];
  if( current_score < min ) min = current_score;

  if( current_score <= final_score ){
    if( trace == 'M'){ i--; j--; }
    if( trace == 'U' ){ i--; }
    if( trace == 'L' ){ j--; }
    if( trace != '.'){
      return last_minimum( i,j, final_score, min);
    }
  }
  return min;
  }



list<column> local_aligner::trace_back( int i, int j , list<column>& newAli, const double& min ,const double& max ){

  char trace = T[i][j];
  double current_score = S[i][j];
  

  if( ( current_score > min ) &&  
	//(current_score == min && (i == 1 || j == 1)) ) &&
         trace != '.' ){
    if( trace == 'M'){
      T[i][j] = '.';  // delete the trace (M := MATCH);
      i--; j--;
      double col_score = current_score - min;
      column c = column(A.get(i), B.get(j) ); // create columns from two others
      c.set_ssp( col_score );
      newAli.push_front(c);
      
    
    }
    if( trace == 'U' ){
      T[i][j] = '.';  // delete the trace (U := UP) ;
      i--;
      double col_score = current_score-min;
      column c1 = column( A.get(i) );
      column c2 = column( B.getSpecies() );
      column c  = column( c1, c2 );
      c.set_ssp( col_score );
      newAli.push_front(c);
      
   
    }
    if( trace == 'L' ){
      T[i][j] = '.';  // delete the trace (L := LEFT);
      j--;
       double col_score = current_score-min;
      column c1 = column( B.get(j) );
      column c2 = column( A.getSpecies() );
      column c  = column( c2, c1 );
      c.set_ssp( col_score);
      newAli.push_front(c);
      
    
    }
    if( trace != '.'){
      return trace_back( i,j, newAli, min, max);
    }
  }
  return newAli;
}



/*************************************************************
 *************************************************************
 **
 **               scoring_scheme
 **
 *************************************************************
 *************************************************************/


/******************************************
 * 
 * Constructor 1 : Default Constructor
 *
 ******************************************/
scoring_scheme::scoring_scheme(){
  base_gap = -2;
  base_missmatch = -3;
  base_match = 2;
  base_threshold = 4;
  mode = "id";  
}


/******************************************
 *
 * Constructor 2 : Constructor for modes: 
 *   - blast
 *   - orthology
 *
 * Read Homology Data int Scoring Table 
 *
 ******************************************/
scoring_scheme::scoring_scheme( string modus, const char* file ){
  
  base_gap        = -2;
  base_missmatch  = -3;
  base_match      = 2;
  base_threshold  = 4;
  ifstream infile( file );
  mode = modus;
 
  string s;
  cerr <<  "Reading " << file << "\n" ;

  while( getline(infile,s) ){    // * creates an 1-species-alignment for each chromosome *

	vector<string> vec = split( s );
	if( mode == "blast" ){
	  string query = vec[0];
	  string target = vec[1];
	  int bits = atoi(vec[2].c_str() );
	  scoring_table[ query ][ target ] = bits;
	  //cerr << query << " " <<  target << " " << bits<< "\n";
	}
	if( mode == "orthologs"){           
          string query = vec[0];
	  string target = vec[1];
	  int bits = 1000;
	  scoring_table[ query ][ target ] = bits;
	  scoring_table[ target ][ query ] = bits;
	}
  }
}

/******************************************
*******************************************
**
** CONSTRUCTOR 3 : Constructor for modes: 
**  - phylo
**
** Read homology data for phylogentic alignments ( Mode:= phylo )
**
********************************************
*******************************************/


scoring_scheme::scoring_scheme( string modus, const char* file, string alignment_tree, string spec_tree ){
  base_gap        = -2;
  base_missmatch  = -3;
  base_match      = 2;
  base_threshold  = 4;
  unsigned int i  = 0;
  ifstream infile( file );
  mode = modus;
  string s;
  // read_species( alignment_tree );  // (CALL OF FUNCTION READ SPECIES)
  set_species_tree( spec_tree );
  cerr <<  "Reading Homology Data: " << file << " " ;

  while( getline(infile,s) ){   
    if( i % 100000 == 0 )cerr << "." ;
    i++;
    vector<string> vec = split( s );
    	if( mode == "phylo" ){
	  string query = vec[0];
	  string target = vec[1];
	  	  
	  int bits = atoi(vec[2].c_str() );
	  scoring_table[ query ][ target ] = bits ; // not symmetric


	}
	
  }

   cerr <<  "Reading " << file << " " ;
  for(map<string, map<string,int> >::const_iterator it = scoring_table.begin(); it!= scoring_table.end();it++ ){
    string id = it->first;
    int score  =  scoring_table[ id ][id];
    if( score == 0 )  scoring_table[ id ][id] = 2000;
    
  }

  cerr <<  "\nFinished\n\n" ;
}

/*************************************************************************
 *************************************************************************
 *
 * ASSOCIATE GENES TO SPECIES FOR PHYLOGENETIC PROGRESSIVE ALIGNMENTMENT
 *
 *
 * Add On for the constructor of the phylogenetic alignment
 * For the phylogentic alignment, each gene corresponds to a species
 * This function maps associated species from the alignment tree to their genes
 *
 *  
 *************************************************************************
 *************************************************************************/


// obsolete
void scoring_scheme::read_species( string s ){    
  
  cerr <<  "Associate Genes with species " << s << "\n";
  species_of = map<string,string>();
  int j = -1;
  vector<string> species_set; 

  for( int i = 0; i < s.size(); i++){  // split "() "
    	
    if( is_paranthese(s[i]) || is_space(s[i])  || i ==  s.size() -1){
      if(  i ==  s.size() -1){
	 i++; // the last string has to be saved separatedly
	//cerr <<  s.substr(j+1,i-j-1) << "\n";
	   species_set.push_back( s.substr(j+1,i-j-2));
	  j = i;
      }
      else{
	   species_set.push_back( s.substr(j+1,i-j-1));
	  j = i;
	}
	}
    }
  
  for( int i = 0; i < species_set.size(); i++ ){ // read data into map
   string species = species_set[i];
   if( species != ""){
     ifstream infile( species.c_str() );
     string s;
     getline(infile,s);
     vector<string> vec = split( s );
     if( vec[0] == "#alignment" ){
       for(  int j = 1; j < vec.size(); j++)
	 read_single_species( vec[j] , species);
     }   
     else{
       read_single_species( species , "" );
     }
   }
 }
  cerr <<  "Finished \n" ;
}

void scoring_scheme::read_single_species( string species , string file){

  if( species !=  "" ) 
    cerr << "Recursively reading " << species << " from " << file << "\n";
  else{
    cerr << "Reading " << species << "\n";
  }
  ifstream infile( species.c_str() );
     string s;
     
     while( getline(infile,s) ){
       vector<string> vec = split( s );
       string gene = vec[0];
       species_of[ gene ] = species;
     }

}



/******************************************************************
 ******************************************************************
 **
 **  Scoring Functions:
 **
 **   S(a,b) = sum ( s(x in a, y in b) ) / |a| |b|  
 **                                                        (1)
 **   s(x,y) = 1-1/ s_bit(x,y) + 1-1 s_bit(y,x)
 **
 **
 **   -    
 **
 **  s(a,b) + s(b,a)
 ** -----------------  < 1
 **  s(a,a) + s(b,b)
 **
 ******************************************************************
 ******************************************************************/


void scoring_scheme::phylogentic_parameter_adaptation( const alignment& X, const alignment& Y ){

  if( mode == "blast" ||  mode == "orthologs"  || mode == "id"  ){
    miss      = base_missmatch;
    gap       = base_gap;
    threshold = base_threshold;
    match     = base_match;
  }
  else{
    if( mode == "phylo" ){
      
      vector<string> spec_X = X.getSpecies();
      vector<string> spec_Y = Y.getSpecies();
      double subtree_size = 0; 
      double scaling_factor = 0;
      
      for( int i = 0; i < spec_X.size() ; i++ ){   // this computes a mean scaling factor for both sets of alignments 
	for( int j = 0; j < spec_Y.size(); j++ ){
	  
	  subtree_size +=  species_tree.get_subtree_size( spec_X[i] , spec_Y[j] );
	  
	}
      }
      
   
      subtree_size = subtree_size / (spec_X.size() * spec_Y.size());
      
      if( spec_X.size() * spec_Y.size() >= 1 ){   // for two species  the subtree size == 1 
	miss      = 1 * base_missmatch * subtree_size;
	gap       = 1 * base_gap * subtree_size;
	threshold = base_threshold * subtree_size; // *subtree_size;
      }
      else{   // at upper nodes 
	miss      = base_missmatch * ( 1 - subtree_size );
	gap       = base_gap * ( 1 -  subtree_size);
	threshold = base_threshold * subtree_size; // *subtree_size;
      }
      cerr << "Distance: "<< subtree_size << "\nTheta:  "<< threshold << "\tGap: "  << gap << "\tMissmatch: " << miss << "\n";
    }
  }

}


double scoring_scheme::scoreBLAST( const column& a, const column& b ){
 
  map<string,location> A = a.get_locs();
  map<string,location> B = b.get_locs();
  double sum = 0;

  for( map<string,location>::iterator it = A.begin(); it != A.end(); it++){
    
    location  pa = it->second;
    
     for( map<string,location>::iterator jt = B.begin(); jt != B.end(); jt++){
 
       
       location  pb = jt->second;
      
       if( pa.getStrand() == pb.getStrand() ){ 
	 if( mode == "phylo"){
	   sum += (scoreGenes( pa, pb ) * 2);
	 
	 }
	 else{
	   sum += scoreGenes( pa, pb );
	   sum += scoreGenes( pb, pa );
	 }
       }
       else{
	 sum += miss/2;
       }       
     }
  }
  
  double comparisons = A.size() * B.size();
  // if(sum >0 ) {
  //  cerr << sum << "\t" <<  comparisons << endl;
  //}
  return sum/comparisons;

}


/*
 * This Function implements the score between two genes for multiple scoring modes
 * -blast
 * -orthology
 * -phylo
 *
 **/

double scoring_scheme::scoreGenes( const location& loc_a, const location& loc_b){ // faster then const reference

  string id_a = loc_a.getID();
  string id_b = loc_b.getID();

  if( mode == "phylo"){
    
    
    if( scoring_table.find(id_a) != scoring_table.end() &&    // test 1. key
	scoring_table[ id_a ].find( id_b) != scoring_table[ id_a ].end()){  // test 2. key
      
      double sab  =  scoring_table[ id_a ][id_b];
      double sba  =  scoring_table[ id_b ][id_a];
      double saa  =  scoring_table[ id_a ][id_a];
      double sbb  =  scoring_table[ id_b ][id_b];
      double ratio = (sab + sba) / (saa+sbb);
    
      if( ratio > 1 ) ratio = 1;   // if no self similarity given
      // if(ratio > 0.5 )
      //cout << id_a << "\t" << id_b << "\t" << ratio << "\td:" << distance( loc_a.getSpecies(), loc_b.getSpecies()) << endl;
      return  ratio * distance( loc_a.getSpecies(), loc_b.getSpecies() ) ; 
    }
    else 
      return miss/2;
  }
  
  else{
    if( mode == "blast" || mode== "orthology")
      if( scoring_table.find(id_a) != scoring_table.end() &&    // test 1. key
	  scoring_table[ id_a ].find( id_b) != scoring_table[ id_a ].end()){  // test 2. key
	return 1 - 1/scoring_table[ id_a ][id_b];
      }
      
    return miss/2;
  }
  
}

double scoring_scheme::scoreID( const column& a, const column& b ){

  map<string,location> A = a.get_locs();
  map<string,location> B = b.get_locs();

  double sum = 0;
  for( map<string,location>::iterator it = A.begin(); it != A.end(); it++){
     for( map<string,location>::iterator jt = B.begin(); jt != B.end(); jt++){

       string id_a = (it->second).getID();  // it->second denotes a location
       string id_b = (jt->second).getID();
       
       if( id_a == id_b ){
	 sum += match;
       }
       else{
	 sum += miss;
       }
     }
  }

  double comparisons = A.size() * B.size();
  return sum/comparisons;
}


 double scoring_scheme::score(const column& a, const column& b){

  if( mode == "phylo"  ||  mode == "orthologs" ||  mode == "blast" )
    return scoreBLAST(a, b);
  else if( mode == "id")
    return scoreID(a, b);
  
}



/**
 *
 * Get all homologous genes for a query ID
 *
 **/

list<string> scoring_scheme::get_all_homologs( string query){

  list <string> ret; 
  if( mode == "id" ){
    ret.push_back( query );
  } 
  else{
 
    //shmap< shmap<double>::map >::map::const_iterator target_map = scoring_table.find(query);  // Hash_map version
        map<string, map<string,int> >::const_iterator target_map = scoring_table.find(query);
    
    if( target_map != scoring_table.end() ){   
      
      map<string, int> targets = target_map->second;
       for( map<string, int>::const_iterator jt = targets.begin(); jt != targets.end(); jt++ ){
      //shmap<double>::map targets = target_map->second;
      //for( shmap<double>::map::const_iterator jt = targets.begin(); jt != targets.end(); jt++ ){
	string homolog = jt->first;
	ret.push_back( homolog );
      }
    }
  }
  return ret;
}


/**
 * get_all_homologs
 *
 * returns a unique list of gene Identifier which are
 * homologs of the genes in the alignments.
 *
 * a) for each gene, get all homologs
 * b) make unique  list
 *
 **/

list<string> scoring_scheme::get_all_homologs( alignment a ) {

  list <string> ret; 
  list <string> IDs = a.getIDs();  
  
  for( list<string>::const_iterator it = IDs.begin(); it != IDs.end(); it++ ){
    list <string> homologs = get_all_homologs( *it );
    copy( homologs.begin(), homologs.end(), back_inserter( ret) );
  }
 
  list<string>::iterator end = unique( ret.begin(), ret.end() );
  list<string> unique_list;
  copy( ret.begin(), end, back_inserter( unique_list ) );
  return unique_list;
}

/*****************************************************************
 *
 * Global functions
 *
 ****************************************************************/

vector< vector<char> >  trac_back_matric( int n, int m ){
  
  vector< vector<char> > T( n+1 , vector<char>(m+1, '.') );
  return T;
}

vector< vector<double> > scoring_matrix( int n, int m ){
  vector< vector<double> > S(n+1,  vector<double>(m+1,0));
  return S;
}



/*****************************************************************************
 *
 *  peak ( local alignment)
 *
 *****************************************************************************/

peak::peak( int a, int b, double c){
  i = a;
  j = b;
  score =c;

}


bool comparePeaks( const peak& a, const peak& b){
  return a.get_score() > b.get_score();
}
