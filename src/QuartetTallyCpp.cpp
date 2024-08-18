#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

//' Build quartet table from distances
//'
//' This is a C++ function, called from quartetTable, to fill in the quartet counts.
//' From a list of topological distance matrices (1 for each gene tree) it determines all
//' gene quartets. It is not intended to be used as a stand-alone function, and hence not fully
//' documented. The faster looping in C++ over R gives substantial time improvements
//'
//' @param dList a list of distance matrices
//' @param M number of sets of 4 taxa
//' @param nt number of gene trees/distance matrices
//' @param Q matrix to fill out as table of quartet counts
//' @param random if 0 compute for all sets of 4 taxa, otherwise for M random ones
//' @param progressbar if TRUE, display progress bar
//' @seealso \code{\link{quartetTable}}, \code{\link{quartetTableParallel}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List quartetTallyCpp(Rcpp::List dList,
                                  int M,
                                  int nt,
                                  Rcpp::NumericMatrix Q,
                                  int random,
                                  bool progressbar = false) {

int warnMissing=0; //flag for some taxa missing
int  numColsQ=0; //number of columns of Q
Rcpp::StringVector qColNames; // column names of Q
Rcpp::StringVector qtt = {"12|34","13|24","14|23","1234"};
std::vector<int> tqc= {0,0,0,0}; //column indices for qtt
Rcpp::StringVector currentTaxa={"a","a","a","a"};
Rcpp::CharacterMatrix Qtax(M,4); //for storing names of taxa in sets

Rcpp::NumericMatrix d; //distance matrix
Rcpp::StringVector dNames; //names on distance matrix
Rcpp::IntegerVector dpos={0,0,0,0}; // indices for current taxa

double a, b, c; //4-point condition quantities
std::vector<double> z= {0,0,0}; //for sorting a,b,c
int  tc1, tc2, tc3, tc4; // taxon numbers in quartet


qColNames = colnames(Q);  //get column names
numColsQ = qColNames.length(); //   and number of them
tqc = {numColsQ-4,numColsQ-3,numColsQ-2,numColsQ-1}; // indices for quartet count columns
int N=numColsQ-4; //number of taxa

if (random==0) //if not considering random quartets, put them all into table
 {
   int w=0;
   for (int i = 0; i<(N - 3);i++)
     for (int j = (i + 1); j<(N - 2);j++)
       for (int k = (j + 1); k<(N - 1);k++)
         for (int l = (k + 1); l<N;l++)
         {
           Q(w,i)   = 1;
           Q(w,j)   = 1;
           Q(w,k)   = 1;
           Q(w,l)   = 1;
           Qtax(w,0) = qColNames[i];
           Qtax(w,1) = qColNames[j];
           Qtax(w,2) = qColNames[k];
           Qtax(w,3) = qColNames[l];
           w++;
         }
 }
 else {  // If random set of 4 taxa....
   int w = 0, j = 0;
   Rcpp::IntegerVector q;

   while (w < random)  {
     q = Rcpp::sample(N, 4); // 4 random taxa, numbered from 0 up
     std::sort(q.begin(), q.end() );

     j = 0; // check earlier rows of Q for repeated set
     while (j < w) {
       if ( (Q(j,q[0])==1) && (Q(j,q[1])==1) && (Q(j,q[2])==1) && (Q(j,q[3])==1) )
         j = w+1; // flag a match, and get out of loop
       else
           j++; // go to next row
     }
     if (j == w) { // accept this choice of 4 taxa
       Q(w,q[0]-1)=1;
       Q(w,q[1]-1)=1;
       Q(w,q[2]-1)=1;
       Q(w,q[3]-1)=1;
       Qtax(w,0) = qColNames[q[0]-1];
       Qtax(w,1) = qColNames[q[1]-1];
       Qtax(w,2) = qColNames[q[2]-1];
       Qtax(w,3) = qColNames[q[3]-1];
       w++;
     }
   }
 }

 // Initializes the progress bar
 Progress p(nt, progressbar);

for (int i = 0; i < nt; ++i) { // run through list of topological distance tables
   d = Rcpp::as<Rcpp::NumericMatrix>(dList[i]);
   dNames = colnames(d);

   for (int j = 0; j < M; ++j) { // for each set of 4 taxa

     dpos=Rcpp::match(Qtax(j,Rcpp::_),dNames); // find taxa position for distance table
     if ( Rcpp::any( Rcpp::is_na(dpos) ) )
       warnMissing=1; // flag one or more missing taxon
     else {
       tc1=dpos[0]-1;
       tc2=dpos[1]-1;
       tc3=dpos[2]-1;
       tc4=dpos[3]-1;

       a = d(tc1, tc2) + d(tc3, tc4); // compute 4-point condition quantities
       b = d(tc1, tc3) + d(tc2, tc4);
       c = d(tc1, tc4) + d(tc2, tc3);

       z[0]=a; // copy for sorting
       z[1]=b;
       z[2]=c;
       std::sort (z.begin(), z.end()); // using default comparison (operator <):
       if (z[0] == z[1])
         Q(j,tqc[3])++; // an unresolved quartet 1234
       else if (z[0] == a)
         Q(j,tqc[0])++; // 12|34
       else if (z[0] == b)
         Q(j, tqc[1])++; // 13|24
       else
         Q(j, tqc[2])++; // 14|23
     }
   }
   // increments the progress bar
   if (Progress::check_abort() )
     return -1.0;

   p.increment(); // update progress
 }

return Rcpp::List::create(Rcpp::Named("table") = Q,
                          Rcpp::Named("missingFlag") = warnMissing);
}


