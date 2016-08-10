
#include<list>
#include<vector>
#include<string>
#include<iostream>
#include<map>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include "flow.h"
#include "genome.h"
#include "localign.h"
#include "species_tree.h"

using namespace std;








string program_flow::align( string subtree , int nr ){
 
  
  string line = "* Aligning subtree " + intToString(nr + 1) + " : " +  subtree + " *";
   string frame(line.size(), '*' );
   cerr << "\n" << frame << "\n";
   cerr << line << "\n";
   cerr << frame << "\n\n";
  int pos[4];
  int j = 0;
 

  /***************
  *
  * Get Alignment Names from PHYLOGENETIC TREE STRING
  *
  ****************/

  int count = 0;
  for( int i = 0; i < subtree.size(); i++  ){
    if( is_space( subtree[i]) || is_paranthese( subtree[i] )  ){  // SEARCH FOR FIRST SPACE
      pos[j++] = i+1;
      //cerr << i+1 << " ";
    }
  }
  cerr << "\n"; 
  
  string alignment1 = subtree.substr( pos[0], pos[1]-pos[0]-1);   // EXTRACT POSITIONS
  string alignment2 = subtree.substr( pos[1], pos[2]-pos[1]-1);
  
  
  string anr = intToString( nr );
  string ret = "alignment" + anr;

  
  alignment_set A, B, C;

  /***************************************************** 
  *
  *  The alignment has either to be read from file or has been computed
  *  and stored in the Alis map
  *
  ******************************************************/

 
  if( Alis.find(alignment1) != Alis.end() )   // Get Alignment: A
    A = Alis[ alignment1 ];  
  else
    A = alignment_set( alignment1.c_str() );

 
  if( Alis.find(alignment2) != Alis.end() )   // Get Alignment: B
    B = Alis[ alignment2 ];
  else 
    B = alignment_set( alignment2.c_str() );   

  // cerr << "A" << A.to_string();
  // cerr << "B" << B.to_string();

  //B.test("ENSCAFG00000019642");

  /***********************************
   *
   * Compute New Alignment:  C
   *
   ************************************/
  
	vector<alignment> a = A.get_alignments();
	//vector<alignment> b = B.get_alignments();
	
	local_aligner  L;
	for( int i = 0; i < a.size(); i++ ){
	  
	  if( i > 0)
	    cerr << i << "/" << a.size() << "\n";
	  alignment q = a[i];
	  
	  list<string> homologs = scorer.get_all_homologs( q );   // * christoph's
	  vector<alignment> b = B.get_alignments_for( homologs ); // * heuristic
	  if( scorer.get_mode() == "id" ){
	     b = B.get_alignments();
	    }
	  
	  
	  for( int j = 0; j < b.size(); j++ ){
	    
	    
	    L.smith_waterman( a[i], b[j], C, scorer);
	    L.parameter_fixed();
	    //local_aligner R;
	    alignment rev = b[j];
	    L.smith_waterman( a[i], rev.reverse_copy() , C, scorer);
	    
	  }
	}
	
	cerr << a.size() << "/" << a.size() << "\n";

	if( C.size() > 0 ){
	  C.sort( filter , max_coverage , min_length); // filter for n best alignments
	  
	  if( print_all ){
	    if( output != "" ){
	       ofstream myfile;
	       myfile.open (output.c_str());
	       myfile << printHeader( C ); 
	       myfile << C.to_string(); 
	       myfile.close();
	       
	    }
	    else{
	      cout << printHeader( C );
	      cout << C.to_string();
	    }
	  }
	  Alis[ret] = C;
	  //cerr <<  alignment1 <<  " - "  << alignment2 << " -> " << ret << "\n"; 
	}
 return ret;
  
}



program_flow::program_flow(int argc,  char** argv){


  scorer = scoring_scheme();
  if( argc < 2 ){
    print_help();
  }
  else{
    
    int check = get_parameters( argc, argv); // homology and species tree
    
    if( check >= 2 ){
      
      vector<string> order = get_alignment_order( tree );
      
      for( int i = 0; i < order.size(); ++i ){
	
	      string s    =  order[i];

	      if( i == order.size()-1 )    // in the last round print_all is forced to be set to 1 
	        print_all = 1;

	      string repl =  align(s, i);
	
	      for(int j = i; j < order.size(); ++j ){
	        order[j] = replace( order[j], s, repl );
	      }	
      }     
    }
    else{
      print_help();
    }    
  }  
}

/**
 * Use the first alignment to print an alignment header
 * containing alphabetically sorted species named 
 *
*/


string printHeader( alignment_set A ){

  string head = "#alignment ";
  int i = 0;
  vector<alignment> vec = A.get_alignments ();
  alignment a = vec[0];
  vector<string> specs = a.getSpecies();
  for( i; i < specs.size() - 1; i++ ){
    head = head + specs[i] + " ";
  }
  
  head = head + specs[i] +  "\n"; 
  return head;
}


int program_flow::get_parameters( int argc,  char** argv ){

  double t = 4;   // default parameters
  double g = -2;
  double m = -3;
  filter = 100;
  max_coverage = 2;
  min_length = 0; 
  print_all = 1;

    int check = 0;
    for( int i = 1;  i < argc -1 ; i++ ){
      string argument = argv[i];
      string next = argv[i+1];
      if( argument == "-t" ) {
	tree = next;
	check++;
      }
      if( argument == "-h" ) {
	modus = next;
	if( modus != "id"  ){
	                            /*****************************************************/
	  if( modus == "phylo" ){   /***          MODUS: ---PHYLO----                  ***/
                        	    /*****************************************************/
	    scorer  = scoring_scheme( modus, argv[i+2] , tree,  argv[i+3] );//"(((maus:1 ratte:1.2):2 mensch:1.5):3 (Kuh:1.1 Pferd:1.5):2.2)");
	    //cerr <<  "return 0 ****\n";
	    
	  }
	  else{            
	    scorer = scoring_scheme( modus, argv[i+2] );
	  }
	}
	check++;
      }
      
      if( argument == "-last" || next == "-last") { print_all = 0; }
      if( argument == "-thr" ) { t  = atof( next.c_str()) ;}
      if( argument == "-gap" ) { g  = atof( next.c_str()) ;}
      if( argument == "-mis" ) { m  = atof( next.c_str()) ;}
      if( argument == "-filter" ) { filter  = atoi( next.c_str()) ;}
      if( argument == "-coverage" ) { max_coverage  = atoi( next.c_str()) ;}
      if( argument == "-length" ) { min_length  = atoi( next.c_str()) ;}
      if( argument == "-o" ) { output  = next;}
    }
    scorer.set_gap(g);
    scorer.set_miss(m);
    scorer.set_threshold(t);
    
    return check;
}



void program_flow::print_help(){

  string help = "program -t guide-tree -h homology_type ... \n";
  help += "guide-tree:\n\t-t \"((rat.txt mouse.txt ) human.txt)\"\n";
  help += "Homology:\n\t-h id\n\t-h blast [file]\n\t-h orthologs [file]\n";
  help += "\t-h phylo [file] [weighted_tree]\n";
  help += "\nAlignment Parameters:\n\t-thr\tthreshold (4)\n\t-gap\tgap (-2)\n\t-mis\tmismatch (-3)\n";
  help += "\nFilter options:\n\t-filter [int] best alignments or only unique assignments n=0 (100)\n";
  help += "\t-coverage [int] each gene may occur only c times in alignments (2)\n";
  help += "\t-length [int] minimum alignment length treshold (1)\n";
  help += "\t-last prints only the alignments at the last step\n";
  help += "\nOutput:\n\t-o output file\n";
    cerr << help;

}



string replace( string s, string search_string, string with ){

  string ret = "";
  int n = search_string.size();
  bool replaced =0;
  //cerr << "replace in " << s << ": " << search_string << " with " << with<< "\n";
  for( int i = 0; i+n <= s.size(); i++  ){
    
    //cerr << i << ": " << s.substr(i,n) << "\n";
    if( search_string == s.substr(i,n)){
      //cerr << "REPLACED\n";
      ret += with;
      i += n-1;
      n = 0;
      replaced = 1;
    }
    else{
      ret += s.substr(i,1);
    }

  }
  
  if( replaced ){
    //cerr << "returning:"<< ret << "!\n";
  return ret;
  }
  else{
    //cerr << "returning:"<< s << "!\n";
  return s;
  }
}
