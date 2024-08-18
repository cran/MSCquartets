#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

//################# Helper function ################################


// Determine which entries in a vector have a given value
//
//
//
Rcpp::IntegerVector whichVal( Rcpp::IntegerVector x, int W) {
  int nx = x.size();
  Rcpp::IntegerVector y;
  for(int i = 0; i < nx; i++) {
    if (x[i]==W) y.push_back(i);
  }
  return y;
}

//#################################################################################
//' Initialize vector of B quartets
//'
//' This is a C++ function, called from \code{TINNIKdist}, to initialize for
//' inference of B and T quartets.
//'
//' @param pTable a quartet table with p-values
//' @param m number of rows in pTable
//' @param alpha critical value for tree test
//' @param beta critical value for star tree test
//' @param colptest column with p value for tree test
//' @param colpstar column with p value for star tree test
//' @param Bquartets 0/1 vector encoding initial B quartets
//' @seealso \code{\link{quartetTable}}, \code{\link{quartetTableParallel}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List initBquartets(Rcpp::NumericMatrix pTable,int m, double alpha, double beta,
                         int colptest, int colpstar,Rcpp::StringVector  Bquartets) {
  Rcpp::NumericVector L1 (m);
  Rcpp::List ret;

  int j=0; // index for list of new Bquartets
  for(int i =0; i < m ; i++)
  {
    if ((pTable(i, colptest-1 ) < alpha) || (pTable(i,colpstar-1) > beta)) {
      Bquartets[i] = true;
      L1[j++]=i+1;
    }
  }

  ret["Lout"] = L1;
  ret["Bout"] = Bquartets;
  return ret;
}


//#########################################

//' Main loop of B-quartet inference
//'
//' This is a C++ function, called from \code{TINNIKdist}, to
//' infer B and T quartets.
//'
//' @param pTable a quartet table with p-values
//' @param C precomputed binomial coefficients
//' @param Cn4 precomputed binomial coefficient
//' @param n number of taxa
//' @param Bquartets 0/1 vector of initial Bquartets
//' @param L1 vector of recently inferred B quartets
//' @param lenL1 lnegth of L1
//' @param Nrule1 count of inference from rule 1
//' @param Nrule2 count of inference from rule 2
//' @param cuttops inferred cut topologies
//'
//' @seealso \code{\link{quartetTable}}, \code{\link{quartetTableParallel}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::LogicalVector BQinference(Rcpp::IntegerMatrix pTable,Rcpp::IntegerMatrix C,
                                int Cn4,int n, Rcpp::LogicalVector  Bquartets,
                                Rcpp::IntegerVector L1, int lenL1, int Nrule1,
                                int Nrule2,Rcpp::IntegerVector cuttops) {

  Rcpp::LogicalVector ret;
  Rcpp::IntegerVector L2={};
  Rcpp::IntegerVector idx = seq(0,n-1);
  int np1= n-1;

  while ( !(lenL1 == 0) ) {

    Rcpp::IntegerVector L2={};// start list of new B-quartets

    for (int i =0; i < lenL1 ; i++) { //for last found B-quartets
      Rcpp::IntegerVector taxind = head(pTable( L1[i]-1, _),n);// get taxon entries from row for quartet
      Rcpp::IntegerVector current4 = whichVal(taxind,1); //#determine taxa in this B-quartet
      Rcpp::IntegerVector others = whichVal(taxind,0); //# and all those not in it

      for (int j =0; j < 4 ; j++) {
        int current4th = current4[j]; //drop one taxa from current quartet
        Rcpp::IntegerVector current3 = current4;
        current3.erase(j);

        for (int k =0; k < (n-4) ; k++) {
          int new4th =  others[k];
          Rcpp::IntegerVector new4=current3;
          new4.push_back(new4th);
          std::sort(new4.begin(), new4.end()); //Create new set of 4 with 3 from current quartet
          int new4index = Cn4 - C(np1 - new4[0], 3) - C(np1 - new4[1], 2) - C(np1 - new4[2], 1) - C(np1 - new4[3], 0) -1;// #index for new4

          if (Bquartets[new4index] == true) { //the new set was a B-quartet, so apply inference rules

            //Get three quartets that can be used in inference rule
            Rcpp::IntegerVector set1={current4th, new4th, current3[0], current3[1] };
            std::sort(set1.begin(),set1.end());
            Rcpp::IntegerVector set2={current4th, new4th, current3[0], current3[2] };
            std::sort(set2.begin(),set2.end());
            Rcpp::IntegerVector set3={current4th, new4th, current3[1], current3[2] };
            std::sort(set3.begin(),set3.end());

            //store these in a matrix for later easy access
            Rcpp::IntegerMatrix sets (3,4);
            sets(0,_) = set1; sets(1,_)=set2; sets(2,_)=set3;

            //Compute indices to look up these sets
            Rcpp::IntegerVector indexsets (3);
            indexsets[0] = Cn4-1 - C(np1 - set1[0], 3) - C(np1 - set1[1], 2) - C(np1 - set1[2], 1) - C(np1 - set1[3], 0); //#index for set1
            indexsets[1] = Cn4-1 - C(np1 - set2[0], 3) - C(np1 - set2[1], 2) - C(np1 - set2[2], 1) - C(np1 - set2[3], 0); //#index for set2
            indexsets[2] = Cn4-1 - C(np1 - set3[0], 3) - C(np1 - set3[1], 2) - C(np1 - set3[2], 1) - C(np1 - set3[3], 0); //#index for set3

            Rcpp::LogicalVector Bs = Bquartets[indexsets];
            if (!(Bs[0] && Bs[1] && Bs[2])) { //if at least one of the 3 sets is not a B-quartet

              bool makeB = 0; // flag for creating new B quartet

              if (Bs[0] || Bs[1] || Bs[2]) { //Use inference rule 2
                Nrule2++;
                makeB = 1;
              } else { // Use inference rule 1, maybe
                for (int s=0; s<3; s++) {
                  Rcpp::IntegerVector sestS = sets(s,_);

                  //locate taxa that need to be checked for non-cherry-ness
                  Rcpp::IntegerVector staxa (2);
                  staxa[0]=whichVal(sestS,current4th)[0];
                  staxa[1]=whichVal(sestS,new4th)[0];
                  std::sort(staxa.begin(),staxa.end());

                  int cutindex = cuttops[indexsets[s]]; //#get inferred topology

                  // Check if rule 1 applies
                  if ( ( (cutindex == 1)&&( !((staxa[0]==0) && (staxa[1]==1)) )&&(!((staxa[0]==2) && (staxa[1]==3)) ) )  ||
                       ( (cutindex == 2)&&( !((staxa[0]==0) && (staxa[1]==2)) )&&( !((staxa[0]==1) && (staxa[1]==3)) ) )  ||
                       ( (cutindex == 3)&&( !((staxa[0]==0) && (staxa[1]==3)) )&&( !((staxa[0]==1) && (staxa[1]==2)) ) ) ) {

                    makeB = 1;
                    Nrule1++;
                  }
                }
              }

              if (makeB) {
                for (int kk=0; kk<3; kk++) {
                  //# if new B-quartet
                  if ( !(Bquartets[indexsets[kk]] == true) ) {
                    Bquartets[indexsets[kk]] = true; //# flag it
                    L2.push_back(indexsets[kk]+1);
                  }//end if
                }//end for
              }//end if


          } //end if

        }//end if

      }//end for

    }//end for

    }//end for

    L1 = L2; //#move list of new B-quartets
    lenL1 = L1.size();

}// end while

return Bquartets;
}

