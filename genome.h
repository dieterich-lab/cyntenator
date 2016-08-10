#ifndef GUARD_genome_h
#define GUARD_genome_h

using namespace std;

#include <string>
#include <vector>
#include <map>
#include <list>
#include<algorithm>

/*************************************************************************
 *
 *  1) alignment_set stores a vector of alignments
 *  2) alignment stores a column vector
 *  3) column stores a vector of loactions
 *  4) location is thought to store a position within a sequence
 *
 *************************************************************************/ 


/*************************************************************************
 *
 * The LOCATION class saves the sequential data
 * such as genes, TFBS, exons ...
 *
 ************************************************************************
*/
class location{

 private:
 int start;
 int end;
 string sequence;
 string id;
 string species;
 bool strand;
 

 public:
 location();
 location( string s ){ id = s; start = 0; end = 0; sequence= "-"; strand = 0;}
 location( string a, string b, int c, int d, bool e, string spec){
   id = a; sequence = b; start = c, end = d; strand = e; species = spec;
 }
 location(vector<string> vec, string spec );
 location( const location& l){
   start = l.start; end = l.end; sequence = l.sequence; id = l.id; strand = l.strand; 
   species = l.species;
 } 

 bool getStrand(){ return strand; }
 string getSpecies() const { return species; }
 string getID()       const { return id;}
 string getSequence() const { return sequence;}
 int getStart() const { return start;}
 int getEnd() const { return end;}
 void setStart(int x){ start =x; }
 void setEnd(int x){ end =x; }
 string to_string();
 location reverseCopy();

};

/*************************************************************************
 *
 * The COLUMN class saves a column within an alignment
 * 
 *************************************************************************
*/

class column{
 private: 
  int numberOfSequences;
  double ssp;
  map<string,location> AlignedPositions;
 
 public:
  column();
  column( location l, string species);
  column( int a, vector<string> spec_nam, vector<string> genes, vector<bool> oris, int rn, double ssp);
  column( vector<string> spec );    // * gap * 
  column( const column&  );         // * copy constructor *
  column( column a, column b);
  column reverseCopy();
  string to_string( vector<string> species);
  string to_string();
  double get_ssp() const{ return ssp;}
  void set_ssp( double x){ ssp = x;}
  location first() const{ return (AlignedPositions.begin())->second;}
  void add(location l, string species  );
  map<string,location> get_locs() const { return AlignedPositions; }
  list<string> getIDs();
};


/*************************************************************************
 *
 * The ALIGNMENT class saves a sequences of
 * Alignment columns
 *
 *************************************************************************
*/
class alignment{

 private:
  double score;
  vector<string> species;
  vector<column> cols;
  
 
 public:
  alignment();
  alignment(char* file);
  alignment( list<column> L , double s);

  alignment reverse_copy();
  int    size(){ return cols.size();}
  void   sort_alignment();
  void   add( location l ,  string spec );
  void   add( string spec);
  void   add( column c ){ cols.push_back(c);}
  void   updateLengths( vector<int>& maxlen);
  list<string> getIDs();
 
  column get( int i ){ return cols[i]; }
  double get_score() const{ return score; }
  vector<string> getSpecies() const{ return species; }
  string to_string();
  string to_padded_string( vector<int> maxlen);
};

/*************************************************************************
 *
 * The ALIGNMENTSET class saves a sequences of
 * Alignments
 *
 *************************************************************************
*/

class alignment_set{

  //int number_of_species;
  //vector<string>    species;
  vector<alignment> alignments;
  map< string, list<int> > genes2alignments;
  int nr;    // current number of highest alignment
  
  
  list<alignment> sort_list();            // create an ordered alignment set
  int max_cov( list<string>& ids,  map<string,int>& gene_coverage );  // maximum coverage alignment set
  void update_cov( list<string>& ids,  map<string,int>& gene_coverage );  
 public:
  //map< string, list<alignment> > genes2alignments;
  alignment_set(){} // number_of_species=0; }
  alignment_set( const char* file );
  vector<alignment> get_alignments () const {return alignments;} 
  void add( alignment& ali);
  void test(string id);
  int size() const{ return alignments.size(); }
  string to_string();
  void sort(int n, int max_cov , int length_threshold );
  vector<alignment> get_alignments_for( list<string> L);

};




string intToString( int i );
bool compareInts( const int& a, const int& b );
bool compareLoc( const location& a, const location& b );
bool compareCol( const column& a, const column& b );
bool compareAlignments( const alignment& a, const alignment& b );
bool equalAlignments( const alignment& a,  const alignment& b );
bool is_space( char c );
bool is_paranthese( char c );
location  max_span(const location& l1, const location& l2);
string determine_input_type_from_first_line( string s );
vector<string> split(string& s );
#endif

