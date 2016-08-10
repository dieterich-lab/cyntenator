#include <string>
#include <vector>
#include <list>
#include<iostream>
#include<algorithm>
#include<iterator>    //iterator
#include<fstream>
#include <stdlib.h>    // atoi
#include <sstream>     // int to string
#include <utility>    // make_pair
#include"genome.h"


using namespace std;


/**
* 
* Returns true if the passed character is a whitespace
*
*/

bool is_space( char c ){

if( c=='\t' || c=='\n' || c==' ')
	return 1;
return 0; 
	
}


bool is_paranthese( char c ){
if( c=='(' || c==')' || c=='[' || c==']' || c=='{' || c=='}' )
	return 1;
return 0; 
}

string intToString( int i ){

string s;
stringstream out;
out << i;
s = out.str();
 return s;
}

string doubleToString( double i ){

string s;
stringstream out;
out << i;
s = out.str();
 return s;
}

/**
*  
* Splits a String at all white space positions
* and returns a vector of words
*
*      j      i
*      |     |  
*      V     V
*  AAA  BBBBB   CCCC   
* j   i                => save "AAA"      
*     ji               => save ""           
*/
 // moo koenig version
/*
vector<string> split(string& s ){
int j = -1;
vector<string> vec; 

for( int i = 0; i < s.size(); i++){  
    if( is_space(s[i])  || i == s.size() - 1){
      if(  i ==  s.size() -1){
	 i++; // the last string has to be saved separatedly
	  vec.push_back( s.substr(j+1,i-j-1));
	  j = i;
      }
      else{
	  vec.push_back( s.substr(j+1,i-j-1) );
	  j = i;
	}
	}
    }

 

return vec;	
}
*/

vector<string> split(string& s ){
int j = -1;
 vector<string> vec; 
 
 for( int i = 0; i < s.size(); i++){

   // j marks position of the last white space 
   
   if( is_space(s[i]) ){  
     if( i > 0 && !is_space(s[i-1])){
       vec.push_back( s.substr(j+1,i-j-1) );
       j = i;
     }
     else{
       j++;
     }
   }

   if( i == s.size()-1 && !is_space(s[i]) ){
     vec.push_back( s.substr(j+1,i-j) );
   }
 }
 
   

 // for( int i = 0; i < vec.size(); i++){
 //   cerr << i << "\t|" << vec[i] << "|\n";
 //   }
 return vec;
}
 



/* **********************************************************************************
 *
 *                                     ALIGNMENT SET 
 *
 *
 * The constructor alignment_set( const char* file) Reads either an annotated sequence
 * file, a gene annotation file or an existing alignment file with the species names
 * denoted in the first row: #alignment human mouse ...
 * In case of the phylogenetic alignment modus, the species are read in the constructor
 * scoring_scheme::scoring_scheme() , which evokes the read species function. This 
 * associates each gene with a species in an 1:1 mapping. ATTENTION: in order to perform 
 * self-alignments gene identifiers and species names have to be modified.
 *
 ************************************************************************************/


alignment_set::alignment_set( const char* file ){

  //number_of_species = 1;     // read one species with multiple seqs from file
  //species.push_back( file );
  nr = 0;
  ifstream infile( file );

  if (infile.fail())
  {
    std::cerr << "Previous Alignment was empty (file: '" << file << " not found)! Aborting computation for this subtree!\n";
    return;
  }

  map <string, alignment> tmp;
 
  string input_type = "genome";
  string s;
  cerr <<  "Reading " << file << "\n" ;
  
  /*
   * READ THE FIRST LINE AND DETERMINE INPUT FORMAT
   * 
   * the difference between SEQUENCE and GENOME is,
   * that the gene lists are unique ID lists.
   * 
   *
   */

  getline(infile,s);  
  input_type=determine_input_type_from_first_line(s);
  vector<string>  species_set = split( s );	




if( input_type=="sequence" ){

   /*
   * for each chromosome,
   * one-species-alignment is created.
   * The file is read and split into locations, which
   * are saved in the preliminary alignment set "tmp"
   *
   */

  while( getline(infile,s) ){    // * creates an 1-species-alignment for each chromosome *
	
		vector<string> vec = split( s );
		location l = location(vec, file);      
		string seq = l.getSequence();
		tmp[ seq ].add( l , file );  // add location to the chromosome associated alignment 
  		}
	}

 if( input_type=="genome" ){  
	
	/*
	* Save only one location for each gene whereby
	* gene locations are merged as the maximum spanning	
	* segment for all transcripts. Then all unique
	* locations are saved in the alignments, if genomes
	* are not unique, genes are merged to max_spanning
	* segments
	*
	*/

	map<string, location> unique_locs;
	 while( getline(infile,s) ){  
		vector<string> vec = split( s );	
		location l = location( vec, file ); // also stores the species
		string id = l.getID();
		if( unique_locs.find(id) == unique_locs.end() ){
			unique_locs[ id ] = l; 		
			}
		else{
			unique_locs[ id ] = max_span(l, unique_locs[ id ]);
			}
		}
	
	for(map<string, location>::const_iterator it = unique_locs.begin(); it != unique_locs.end(); it++ ){
		pair<string, location > p = *it;
		string seq = (p.second).getSequence();
		tmp[ seq ].add( p.second , file );  // add location to the chromosome associated alignment
		}
	}
	

  /*
   * Sort the positions in each alignment 
   * and save the alignment in the vector 
   * "alignments"
   *
   */

 if( input_type=="genome" || input_type=="sequence"  ) {
  for( map<string, alignment>::const_iterator it = tmp.begin(); it != tmp.end(); it++ ){
    pair<string, alignment> p = *it;   // * sorts and stores the alignments *
    (p.second).sort_alignment();
    add( p.second );
    }
 }
  


 /*
  * Read Alignment Format                               
  *                                                     
  * #alignment                                          
  * Alignment 1 110.733                                 
  * 1.99855 ENSCAFG00000011388 (+)  ENSG00000146950 (+)     
  * 3.98784 ENSCAFG00000011393 (+)  ENSG00000189118 (+)    
  *                                                     
  *                                                     
  */

 if( input_type=="alignment" ){
   
   alignment       a;
   int             number_of_species = 0;
   int             row_number = 0;
   vector<string>  species_names;
   vector<string>  identifier;
   vector<bool>    strands;
   double          ssp = 0;
   list<column>    cols;
   bool            first = 1;
   double          score = 0;
   for( int i = 1; i< species_set.size();i++  )
      species_names.push_back( species_set[i]);
   
   
    while( getline(infile,s) ){
      
      //cerr << s;
     vector<string> vec = split( s );
     if( vec.size() > 0 ){         
       
       /*
       for( int k=0; k < vec.size() ; k++ ){
	 cerr << vec[k] << "\t";
       }
       cerr << endl;
       
	* ALIGNMENT HEADER
	*
	* IF NOT THE FIRST THE ALIGNMENT,
	* SAVE THE CURRENT ALIGNMENT AND
	* REINITIALIZE THE ALIGNMENT PARAMETERS
	*
	*/
 
       if( vec[0] == "Alignment"  ){   // ALIGNMENT HEADER  
	 if( !first ){                   
	   a = alignment( cols, score ); 
	   add(a);                       
	 }
	                                
	 cols  = list<column>();                       
	 score = atof(vec[2].c_str());
	 row_number = 0 ;
	 first = 0;
	 ssp = 0;
       }

       else{                         // READ ALIGNMENT COLUMN
	 identifier = vector<string>();
	 strands    = vector<bool>();
	 number_of_species = (vec.size() - 1) / 2;

	 double ssp = atof(vec[0].c_str());

      
	 for( int i = 1; i < vec.size(); i += 2 )                // read gene identifier
	   identifier.push_back( vec[i] );
	 
	 for( int i = 2; i < vec.size(); i += 2 )                // read strand information
	   if( vec[i] == "(+)" ||  vec[i] == "(+)\t"  )
	     strands.push_back( true );
	   else
	     strands.push_back( false );
	       
         column c( number_of_species, species_names , identifier, strands, row_number, ssp);
	 cols.push_back(c);
	 row_number++;

       }
     }     
    }

   if( !first ){
	   a = alignment( cols, score );                        // SAVE THE LAST ALIGNMENT 
	   add(a); 
	 }
 }

}


location  max_span(const location& l1, const location& l2){
	
	location mspan = location(l1);
	if( mspan.getStart() > l2.getStart() ){
		mspan.setStart( l2.getStart() );
	}
	if( mspan.getEnd() < l2.getEnd() ){ 
                mspan.setEnd( l2.getEnd() );
        }

	return mspan;

}



/*******************************************************
*
* Input Type
*	- genome := non_redundant gene sets, leftmost position
*	- sequence := redundant set of alignment positions ( redundant )
*	- alignment : output of previous alignment
*
*******************************************************/

string determine_input_type_from_first_line( string s ){

	string input_type="";
	
	vector<string> vec = split( s );
	
	if( vec[0] == "#genome" ){
        	input_type="genome";
        	}
	if( vec[0] == "#sequence" ){
                input_type="sequence";
                }
	if( vec[0] == "#alignment" ){
                input_type="alignment";
                }
	return input_type;

}


/********************************************************
 *                                                      *
 * Function:   add                                      *
 * Parameter:  alignment ali                            *
 * Return :    void                                     *
 *                                                      *
 * Add an alignment to the current alignment set        *
 * and update the gene2alignment data structure         *
 *                                                      *
 ********************************************************/

void alignment_set::add( alignment& ali ){
  alignments.push_back( ali );  
  list<string> IDs = ali.getIDs();
  for( list<string>::const_iterator it = IDs.begin(); it != IDs.end(); it++ ){
    list<int> L = genes2alignments[ *it ];
    L.push_back( nr );  
    genes2alignments[ *it ] = L;
    
  }
nr++;
}

/********************************************************
 *                                                      *
 * Function:   get_alignments_for                       *
 * Parameter:  list<string> L)                          *
 * Return :    vector<alignments>                       *
 *                                                      *
 * Retrieves all alignments for the query set L         *
 * and returns a unique alignment set                   *
 *                                                      *
 *                                                      *
 ********************************************************/

vector<alignment> alignment_set::get_alignments_for( list<string> L){
  vector<alignment> ret;
  list<int> tmp;
 
  /*
   * for each id and 
   *      for each corresponding alignment
   *          save the address in tmp 
   */

  for( list<string>::const_iterator it = L.begin(); it != L.end(); it++ ){
    string id = *it;
    map< string, list<int> >::const_iterator jt = genes2alignments.find( id );
    if( jt != genes2alignments.end() ){
      list<int> l =  jt->second;
      
      copy( l.begin(), l.end(), back_inserter(tmp) );
    } 
  }
      
  
      /*
       * make unique
       */


    //list<int>::iterator end = unique( tmp.begin(), tmp.end(), compareInts  );
  map<int,bool> found;
  for( list<int>::const_iterator it = tmp.begin(); it != tmp.end(); it++ ){
    int index = *it;
    if( found.find(index) == found.end() ){
      found[index] = true;
      ret.push_back( alignments[ index ] );
    } 
  }
   
    return ret;
  
}



string alignment_set::to_string(){

  if( alignments.size() == 0) return "";

  vector <int> maxlen;

  list<alignment> l;    // l contains all sorted alignments
  for( int i = 0 ; i != alignments.size(); i++ ){
    l.push_back( alignments[i] );
    alignments[i].updateLengths(  maxlen );
  }
  l.sort( compareAlignments );
 
  string ret;
  int i = 0;
  //for( int i = 0; i < alignments.size(); i++){
  for( list<alignment>::const_iterator it = l.begin(); it != l.end(); it++ ){
    
      alignment a = *it;  
    ret += "Alignment "; 
    ret += intToString( i+1 );
    ret += " ";
    ret += doubleToString( a.get_score());
    ret += "\n";
    //ret += a.to_string();
    ret += a.to_padded_string( maxlen );
    ret += "\n";
    i++;
  }
  
  return ret;
}



void alignment::updateLengths( vector<int>& maxlen){

 for( int i = 0; i < cols.size(); i++ ){
    column c = cols[i];
    string line = c.to_string( species );
   
    vector<string> words = split( line );
    for( int j = 0; j < words.size(); j++ ){
      int n = words[j].size();
      if( maxlen.size() <= j){
	maxlen.push_back( n+2 );
      }
      else{
	if( maxlen[j] < n+2 )
	  maxlen[j] = n+2;
	
      }
    }
  }

}



/*
 *  Create a list copy of the alignment vector and sort it
 *
 **/

list<alignment>  alignment_set::sort_list(){
    list<alignment> l ;
    for( int i = 0 ; i != alignments.size(); i++ ){
      l.push_back( alignments[i] );
    }
    l.sort( compareAlignments );
    
    return l;
  }


 /* iterate over all genes and compute the maximum coverag of any gene in the alignment */
int alignment_set::max_cov( list<string>& ids, map<string,int>& gene_coverage  ){

  int maxcount = 0;
	for( list<string>::const_iterator it =ids.begin();it != ids.end();it++){
	  string gene_id = *it;
	  //cerr << ":" << gene_id << ":" << endl;
	  if( gene_id != "-")
	      if(gene_coverage.find( gene_id ) != gene_coverage.end() ){ 
		if( maxcount < gene_coverage[ gene_id ])
		  maxcount =  gene_coverage[ gene_id ];
	      }	 
	}
   return maxcount;
}

/*
 * update coverages 
 * increment the coverage for each gene in the list ofe ids
*/

void alignment_set::update_cov( list<string>& ids,  map<string,int>& gene_coverage ){

    for( list<string>::const_iterator it =ids.begin(); it != ids.end(); it++){
	    string gene_id = *it;
	    //cerr << ":" << gene_id << ":" << endl;
	    if( gene_id != "-")
	      if( gene_coverage.find( gene_id ) != gene_coverage.end() ){ // increment the coverage
		      gene_coverage[ gene_id ]++;
		      //cout << gene_id << "\t" << gene_coverage[gene_id ] << endl;
	      }
	      else{                                  // insert this gene in the table
		      gene_coverage.insert( make_pair( gene_id, 1) );
		      gene_coverage[ gene_id ] = 1;
		       //cout << gene_id << "\t" << gene_coverage[gene_id ] << endl;
	      }
    }
}

/**
 *  Store up to n alignments for which each gene does not occur
 *  more than maximal_coverage times.
 *
 **/
void alignment_set::sort(int n, int maximal_coverage, int length_threshold ){

  list<alignment>    l                    = sort_list();    // l contains all sorted alignments
  map<string,int>    gene_coverage        = map<string,int>();  // count how often a gene occurs
                     genes2alignments     = map< string, list<int> >(); 
                     nr                   = 0;
                     alignments           = vector<alignment>();   // reinitialize alignments
  int                lsize                = l.size();
  list<alignment>::iterator iter          = l.begin();
  cerr <<"Sort and filter (" << n << ") " << l.size() << " Alignments\n";
 
 /* if the maximum coverage is zero, then we can save the alignment and increment gene coverages */
 
    for(int j = 0; j < lsize ; j++ ){         // iterate over all alignments   
      list<string> ids  = iter->getIDs();
      int maxcount = max_cov( ids , gene_coverage );                          
      if( maxcount == 0 && iter->size() >= length_threshold ){    
	  add( *iter );          // copy alignments consisting only of new genes
	  iter = l.erase(iter);  // delete the first element and place the iterator on the next element
	  update_cov( ids,  gene_coverage);
	}
	else{
	  iter++;
	}
    }
      
     cerr <<"Assign remaining " << l.size() << " non-unique Alignments\n";
     iter          = l.begin();
     if( n > 0 ) {
       while( iter != l.end()  ){	 
	 list<string> ids  = iter->getIDs();
	 int maxcount = max_cov( ids , gene_coverage );                                	 
	 if( maxcount < maximal_coverage && alignments.size() < n  && iter->size() >= length_threshold ){      // add the alignment
	   add( *iter );
	   iter = l.erase(iter);  // delete the first element and increment the iterator
	   update_cov( ids,  gene_coverage);
	 }
	 else
	   iter++;
       }
     } 
}


/* **********************************************************************************
 *
 *                                     ALIGNMENT 
 *
 ***********************************************************************************/



alignment::alignment(){
  score = 0;
  species = vector<string>() ;
  cols = vector<column>(); 

}

/*
 * The constructor for a list of aligned pairs
 * as computed in trace_back of the localigner
 *
 */

alignment::alignment( list<column> L , double s){

  score = s;
  for( list<column>::const_iterator it = L.begin(); it != L.end(); it++ ){

    if( it == L.begin() ){  // * initialize species *
      map<string,location> AlignedPositions = it->get_locs();  // get all locs at first position
      for(map<string,location>::const_iterator jt = AlignedPositions.begin(); jt != AlignedPositions.end(); jt++   ){
	species.push_back( jt->first );
      }
    }
    cols.push_back( *it );
    
    
  }


}

/*
 * WARNING: this function works only for the initialization
 * of an alignment. It sorts the aligned columns according
 * to the locations of the location of the first position
 * in the AlignedPositions map.
 * It is not guaranteed that this function works for alignments
 * with more than one species.
*/

void alignment::sort_alignment(){   // sort can be defined only once
 
  sort( cols.begin(), cols.end(), compareCol ); 

  }


string alignment::to_string(){

  string ret;
 
  for( int i = 0; i < cols.size(); i++ ){
    column c = cols[i];
    ret += c.to_string( species );
  }

  return ret;
}

string alignment::to_padded_string( vector<int> maxlen){

  string ret;
 
  for( int i = 0; i < cols.size(); i++ ){
    column c = cols[i];
    string line = c.to_string( species );
    vector<string> words = split( line );
    for( int j = 0; j < words.size(); j++ ){
      int m = words[j].size();
      ret += words[j] + string( maxlen[j]-m , ' ');
      }
    ret += "\n";
  }
  
  return ret;
}


/*
 *  Adds a location to an alignment
 *  if the corresponding species/file is not yet in the alignment
 *  it is added as well
*/

void alignment::add( location l , string spec ){
  
  add( spec);
  column c = column( l, spec );
  cols.push_back( c );

}

/*
 * Uses a linear search to look for the spec 
 * in the species vector and adds it to the 
 * alignment species vector if it is new.
 *
 */

void alignment::add( string spec ){
  bool found = 0;
  for( int i = 0; i < species.size(); i++ ){
    if( species[i] == spec)
      found = 1;
  } 
  if( found == 0)
    species.push_back( spec);

}


bool compareCol( const column& a,  const column& b){
  
  location l1 = a.first();
  location l2 = b.first();
  return compareLoc( l1, l2 );

}

alignment alignment::reverse_copy(){

  alignment ret;
  ret.score = score;
  for( vector<string>::const_iterator it = species.begin();
       it != species.end();
       it++
       )
    ret.species.push_back( *it );
   
  vector<column>::const_iterator jt = cols.end();
  jt--;
  bool continue_copy = 1;
  while( continue_copy ){
    column c = *jt;
    ret.cols.push_back( c.reverseCopy() );
    if( jt == cols.begin() )
      continue_copy = 0;
    else
      jt--;

  }

  return ret;


}


/*
 * GET ALL IDS FOR EACH COLUMN
 *
 *
*/

list<string> alignment::getIDs(){

list<string> ret; 
for( int i = 0; i < cols.size(); i++){
        list<string> ids = cols[i].getIDs();
        copy( ids.begin(), ids.end(), back_inserter(ret) );     
        }
return ret;
}

/************************************************************************************
 *
 *                                     COLUMN
 *
 ***********************************************************************************/

column::column(){
  ssp = 0; 
 numberOfSequences = 0;
 AlignedPositions = map<string,location>();
}


/**
 *  
 * This is the constructor to read an alignment column
 * from a previously written alignment 
 *
 **/

 column::column( int species_number , vector<string> spec_nam, vector<string> genes, vector<bool> oris, int rn, double ssp2){

   string sequence = "alignment";
   int start = rn;
   int end = rn + 1;
   vector<location> locs;
   ssp = ssp2;
   numberOfSequences = species_number;

   for( int i = 0; i < species_number; i++ ){
     location l( genes[i], sequence, start, end, oris[i], spec_nam[i]  );
    
     AlignedPositions[ spec_nam[i] ] = l;  // species -> location
   }
   
}

column::column( const column& a ){

  int sp = 0;
  
  for(map<string,location>::const_iterator it = a.AlignedPositions.begin(); it != a.AlignedPositions.end(); it++   )
    add( it->second, it->first );
  sp++;
  ssp = a.get_ssp();
  numberOfSequences = sp;
}



/*
 *
 * Gap constructor
 *
 */

column::column( vector<string> spec ){
 
  ssp = 0;
  for( int i = 0; i < spec.size(); i++ ){
    location l = location( "-" );
    add( l, spec[i] );
  }

}


column::column( column a, column b ){

  for(map<string,location>::const_iterator it = a.AlignedPositions.begin(); it != a.AlignedPositions.end(); it++   ){
    add( it->second, it->first );
  }
  for(map<string,location>::const_iterator it = b.AlignedPositions.begin(); it != b.AlignedPositions.end(); it++   ){
    add( it->second, it->first );
  }

  

}

void column::add(location l, string species){

  numberOfSequences++;
  pair<string,location>  p(species, l);
  AlignedPositions.insert( p );

}




column::column(location l, string species){

  numberOfSequences++;
  pair<string,location>  p(species, l);
  AlignedPositions.insert( p );

}


column column::reverseCopy(){

  column ret;
  ret.numberOfSequences = numberOfSequences;
  for( map<string,location>::const_iterator i = AlignedPositions.begin() ; i != AlignedPositions.end(); i++){
    
    string spec = i->first;
    location l = i->second;
    ret.add( l.reverseCopy() , spec);
    
  }

  return ret;

}


string column::to_string( vector<string> species){
  
  string line;
  line += doubleToString(ssp) + "\t";
  
  
  
  for( int i = 0; i < species.size(); i++ ){
    
    string spec = species[i];
    map<string,location>::iterator it;
    it = (AlignedPositions.find( spec ));
    if( it != AlignedPositions.end() ){
      //line += spec + " ";                // print species
      //line +=  intToString( (it->second).getStart() ) + "\t"; // print location 
       line +=  (it->second).getID();
       bool strand = (it->second).getStrand();
       
       if( (it->second).getID() != "-")
	 if( strand )
	   line +=  " (+)";
	 else
	   line +=  " (-)";
       else{
	 line +=  " (.)";
       }
       
       line += "\t";
    }
    

  }
  line += "\n";

  return line;
}

string column::to_string(){
  
  string line;
  for( map<string,location>::const_iterator it = AlignedPositions.begin(); it != AlignedPositions.end(); it++){
    pair<string,location>  pa = *it;
    line +=  (pa.first)+ ":" + (pa.second).getID();
    line += "\t";
}
  return line;
}



/*
 * Iterate icer The Aligned Positions and return a list with all IDs
 *
*/

list<string> column::getIDs(){

  list<string> ret;

  for( map<string,location>::const_iterator it = AlignedPositions.begin(); it != AlignedPositions.end(); it++){
    pair<string,location>  pa = *it;
    ret.push_back ( (pa.second).getID());
    //ret.push_back ( (pa.first)+ ":" + (pa.second).getID());
  }
  return ret;   
}


/************************************************************************************
 *
 *                                     LOCATION 
 *
 ***********************************************************************************/


location::location( vector<string> vec , string spec ){

  location();
  species = spec;
  //  cerr << vec[0] << "\t" <<species <<"\n";
  if( vec.size() == 5){
    //cerr << vec[0] <<"\n";
    id        = vec[0];
    sequence  = vec[1];
    start     = atoi(vec[2].c_str());   // string to int conversion
    end       = atoi(vec[3].c_str());
    strand;
    if( vec[4] == "+")
      strand = 1;
    else
      strand = 0;
    //    cerr <<  vec[0] << " " << vec[4] << " " <<  strand <<"\n";
  }
  else{
    ;//cerr <<  vec[0] <<" wrong vec \n";
  }
}




location::location(){
sequence = "not_assinged";
 id = "not_assinged";
 start = 0;
 end = 0;
 strand = 1;
  
}

location location::reverseCopy(){
  location ret;
  ret.id = id;
  ret.start = - start;
  ret.end = - end;
  ret.sequence = sequence;
  ret.species = species;
  if( strand == 0)
    ret.strand = 1; 
  else
    ret.strand = 0;

  return ret;
}


string location::to_string(){

  string ori = "-";
  if( strand ) ori = "+";
  string s =  id + "\t" + sequence  + ":" + intToString(start) + ":" + intToString(end) + ":" + ori + "\t"; 
  return s;
} 


bool compareLoc( const location& a, const location& b ){

  string aseq =  a.getSequence();
  string bseq =  a.getSequence();
  
  if( aseq != bseq ){
    return aseq < bseq;
  }
  else{
    if( a.getStart() != b.getStart() ){
      return a.getStart() < b.getStart();
    }
    else{
      return a.getID() < b.getID();
    }
  }
}

bool compareAlignments( const alignment& a,  const alignment& b){
  
  return a.get_score() > b.get_score();

}


 bool equalAlignments( const alignment&  x,  const alignment& y ){

   alignment a  = x;
   alignment b  = y;

   if( a.get_score() != b.get_score() ) return false;
   if( a.to_string() != b.to_string() ) return false;
   return true;
   
}


bool compareInts( const int& a, const int& b ){

  if( a == b ) return true;
  else return false;
}
