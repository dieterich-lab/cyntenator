#ifndef GUARD_mymatrix_h
#define GUARD_mymatrix_h
#include <vector>
#include <map>
#include <string>



using namespace std;




template <class T>
class matrix{
  allocator<T> A;
  T* mp;
  bool init;
  int rows;
  int cols;
 public:
  
  matrix(){
    init = 0;
    }
  
  matrix(int n, int m){ 
    int size = n*m;
    //mp = new T[ size ]; 
    mp = A.allocate(size );
    rows = n; 
    cols = m ;
    init = 1;
  }

 void  alloc(int n, int m){ 
    //mp = new T[ n*m ]; 
    mp = A.allocate( n*m  );
    rows = n; 
    cols = m ;
    init = 1;
  }

  ~matrix(){ 
    cerr << "Destructor\n"; // << mp[0] << " first element\n"; 
    if(init ){ 
      A.deallocate( mp, (rows*cols) );
      //delete[] mp;
      }
    init = 0;
    }



  T get(int i, int j){
    return mp[ i*rows + j ];}

  void set(int i, int j, T value){
    mp[ i*rows + j ] = value;}

  int size() { return rows; }
  int length() { return cols; }

};
