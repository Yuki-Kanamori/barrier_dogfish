#include <TMB.hpp>
#include <Eigen/Eigenvalues>

// Function to import R list for user-defined Options_vec and Options, packaged as list Options_list in TmbData
//template<class Type>
//struct options_list {
//  vector<int> Options_vec;
//  vector<int> Options;
//  matrix<int> yearbounds_zz;
//  matrix<int> Expansion_cz;
//  matrix<int> overlap_zz;
//  options_list(SEXP x){ // Constructor
//    Options_vec = asVector<int>(getListElement(x,"Options_vec"));
//    Options = asVector<int>(getListElement(x,"Options"));
//    yearbounds_zz = asMatrix<int>(getListElement(x,"yearbounds_zz"));
//    Expansion_cz = asMatrix<int>(getListElement(x,"Expansion_cz"));
//    overlap_zz = asMatrix<int>(getListElement(x,"overlap_zz"));
//  }
//};

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace Eigen;
  using namespace density;

  // Dimensions
  //DATA_INTEGER(n_i);         // Number of observations (stacked across all years)
  DATA_INTEGER(n_s);         // Number of "strata" (i.e., vectices in SPDE mesh)
  //DATA_INTEGER(n_g);         // Number of extrapolation-grid cells
  //DATA_INTEGER(n_t);         // Number of time-indices
  //DATA_INTEGER(n_c);         // Number of categories (e.g., length bins)
  //DATA_INTEGER(n_e);         // Number of error distributions
  //DATA_INTEGER(n_p);         // Number of dynamic covariates
  //DATA_INTEGER(n_v);          // Number of tows/vessels (i.e., levels for the factor explaining overdispersion)
  //DATA_INTEGER(n_l);         // Number of indices to post-process
  //DATA_INTEGER(n_m);         // Number of range metrics to use (probably 2 for Eastings-Northings)

  // Config
  //DATA_STRUCT( Options_list, options_list );
  
  // SPDE objects
  // DATA_STRUCT(); Get data list object from R and make it available in C++.
  DATA_STRUCT(spde,spde_t);

  PARAMETER(logkappa1);
  //PARAMETER(logkappa2);


  // Random field probability
  Eigen::SparseMatrix<Type> Q1( n_s, n_s );
  //Eigen::SparseMatrix<Type> Q2( n_s, n_s );
  GMRF_t<Type> gmrf_Q;

  Q1 = Q_spde(spde, exp(logkappa1));
  //Q2 = Q_spde(spde, exp(logkappa2));

  return Q1;
  //return Q2;
}