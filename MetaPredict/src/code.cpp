//put_na
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame put_na(DataFrame recon, DataFrame pred) {
  for (int column = 1; column < recon.length(); ++column) {
    IntegerVector reconCol = recon[column];
    NumericVector predCol = pred[column];

    for (int x = 0; x < reconCol.length(); ++x) {
      if (reconCol[x] == 1) {
        predCol[x] = NA_REAL;
      } else {
        predCol[x] = predCol[x];
      }
    }
  }
  return pred;
}



#include <Rcpp.h>
using namespace Rcpp;

// Anywhere the reconstructed result is NA (e.g., the module is not present), replace the NA with the prediction score
// [[Rcpp::export]]
DataFrame put_pred(DataFrame result, DataFrame pred) {
  for (int column = 1; column < result.length(); ++column) {
    CharacterVector resultCol = result[column];
    CharacterVector predCol = pred[column];

    for (int x = 0; x < resultCol.length(); ++x) {
      if (resultCol[x] == NA_STRING) {
        resultCol[x] = predCol[x];
      } else {
        resultCol[x] = resultCol[x];
      }
    }
  }
  return result;
}



#include <Rcpp.h>
using namespace Rcpp;

// If the result column value is "1", change it to the string "Present", otherwise change it to NA
// [[Rcpp::export]]
DataFrame reclassify(DataFrame result) {
  for (int column = 1; column < result.length(); ++column) {
    CharacterVector resultCol = result[column];

    for (int x = 0; x < resultCol.length(); ++x) {
      if (resultCol[x] == "1") {
        resultCol[x] = "Present";
      } else {
        resultCol[x] = NA_STRING;
      }
    }
  }
  return result;
}



#include <Rcpp.h>
using namespace Rcpp;

// If the kegg_df column value is NA, change it to 0, otherwise leave it
// [[Rcpp::export]]
DataFrame na_to_zero(DataFrame kegg_df) {
  for (int column = 0; column < kegg_df.length(); ++column) {
    IntegerVector keggCol = kegg_df[column];

    for (int x = 0; x < keggCol.length(); ++x) {
      if (keggCol[x] == NA_INTEGER) {
        keggCol[x] = 0;
      } else {
        keggCol[x] = keggCol[x];
      }
    }
  }
  return kegg_df;
}



#include <Rcpp.h>
using namespace Rcpp;

// Round values of prediction columns in prediction tibble
// [[Rcpp::export]]
DataFrame round_predictions(DataFrame pred_df) {
  for (int column = 1; column < pred_df.length(); ++column) {
    NumericVector predCol = pred_df[column];
    predCol = round(predCol, 4);
  }
  return pred_df;
}
















