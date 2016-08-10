#include <ext/hash_map>

using std::string;

struct string_hasher_op {
   size_t operator()(const string &s) const {
      int n= s.size();
      size_t x1=0;
      for (int i=0; i<n; i++){
         x1= 31*x1 + s[i];
      }
      return x1;
   }
};

struct string_hasher_eql_op {
   bool operator()(const string &a, const string &b){
      return a==b;
   }
};

// http://www.artima.com/forums/flat.jsp?forum=226&thread=142909&start=45&msRange=15

template<class T> class shmap {
public:
   typedef __gnu_cxx::hash_map<string,T,string_hasher_op,string_hasher_eql_op> map;
   typedef typename map::const_iterator cit;
   typedef typename map::iterator it;
};
