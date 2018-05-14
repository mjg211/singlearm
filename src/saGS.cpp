#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix saGS(int J, double pi0, double pi1, double alpha, double beta,
                   int Nmin, int Nmax, int futility, int efficacy, int equal_n,
                   int ensign, int summary) {
  NumericMatrix dbinomial_pi0(Nmax, Nmax - 1), dbinomial_pi1(Nmax, Nmax - 1);
  for (int j = 0; j <= Nmax - 2; j++) {
    for (int i = 0; i <= j + 1; i++) {
      dbinomial_pi0(i, j) = R::dbinom(i, j + 1, pi0, 0);
      dbinomial_pi1(i, j) = R::dbinom(i, j + 1, pi1, 0);
    }
  }
  NumericVector A_pi0(J), A_pi1(J), R_pi0(J), R_pi1(J);
  NumericMatrix C_pi0(Nmax + 1, J), C_pi1(Nmax + 1, J);
  NumericMatrix feasible_designs(10000000, (J == 2 ? 15 : 20));
  int counter = 0;
  int interrupt = 0;
  if (J == 2) {
    for (int n1 = (futility + efficacy == 2 ? 2 : 1);
         n1 <= (equal_n == 1 ? floor(Nmax/2) : (Nmax - 1)); n1++) {
      C_pi0(_, 0) = NumericVector(Nmax + 1);
      C_pi1(_, 0) = NumericVector(Nmax + 1);
      for (int s = 0; s <= n1; s++) {
        C_pi0(s, 0) = dbinomial_pi0(s, n1 - 1);
        C_pi1(s, 0) = dbinomial_pi1(s, n1 - 1);
      }
      for (int a1 = (futility == 1 ? 0 : -1);
           a1 <= (futility == 1 ? (ensign == 1 ? 0 :
                                     (efficacy == 1 ? n1 - 2 : n1 - 1)) : -1);
           a1++) {
        if (a1 >= 0) {
          NumericMatrix stop_densA_pi0 = C_pi0(Range(0, a1), Range(0, 0));
          A_pi0[0]                     = sum(stop_densA_pi0);
          NumericMatrix stop_densA_pi1 = C_pi1(Range(0, a1), Range(0, 0));
          A_pi1[0]                     = sum(stop_densA_pi1);
        }
        else {
          A_pi0[0] = 0;
          A_pi1[0] = 0;
        }
        interrupt++;
        if (interrupt % 10000 == 0) {
          Rcpp::checkUserInterrupt();
          if (summary == 1) {
            Rcpp::Rcout << "...over " << interrupt  << " designs assessed..." << std::endl;
          }
        }
        if (A_pi1[0] > beta) {
          break;
        }
        for (int r1 = (efficacy == 1 ? n1 : n1 + 1);
             r1 >= (efficacy == 1 ? a1 + 2 : n1 + 1); r1--) {
          if (r1 <= n1) {
            NumericMatrix stop_densR_pi0 = C_pi0(Range(r1, n1), Range(0, 0));
            R_pi0[0]                     = sum(stop_densR_pi0);
            NumericMatrix stop_densR_pi1 = C_pi1(Range(r1, n1), Range(0, 0));
            R_pi1[0]                     = sum(stop_densR_pi1);
          }
          else {
            R_pi0[0] = 0;
            R_pi1[0] = 0;
          }
          if (R_pi0[0] > alpha){
            break;
          }
          for (int n2 = (equal_n == 1 ? n1 : max(1, Nmin - n1));
               n2 <= (equal_n == 1 ? n1 : Nmax - n1); n2++) {
            C_pi0(_, 1) = NumericVector(Nmax + 1);
            C_pi1(_, 1) = NumericVector(Nmax + 1);
            for (int s1 = max(a1 + 1, 0); s1 <= min(r1 - 1, n1); s1++) {
              for (int s2 = 0; s2 <= n2; s2++) {
                C_pi0(s1 + s2, 1) = C_pi0(s1 + s2, 1) +
                  C_pi0(s1, 0)*dbinomial_pi0(s2, n2 - 1);
                C_pi1(s1 + s2, 1) = C_pi1(s1 + s2, 1) +
                  C_pi1(s1, 0)*dbinomial_pi1(s2, n2 - 1);
              }
            }
            for (int a2 = max(0, a1 + 1);
                 a2 <= min(r1 + n2 - 2, n1 + n2 - 1); a2++){
              NumericMatrix stop_densR_pi1 = C_pi1(Range(a2 + 1, n1 + n2),
                                                   Range(1, 1));
              R_pi1[1]                     = sum(stop_densR_pi1);
              if (sum(R_pi1) < 1 - beta) {
                break;
              }
              NumericMatrix stop_densR_pi0 = C_pi0(Range(a2 + 1, n1 + n2),
                                                   Range(1, 1));
              R_pi0[1]                     = sum(stop_densR_pi0);
              if (sum(R_pi0) <= alpha) {
                NumericMatrix stop_densA_pi0 = C_pi0(Range(0, a2), Range(1, 1));
                A_pi0[1]                     = sum(stop_densA_pi0);
                NumericMatrix stop_densA_pi1 = C_pi1(Range(0, a2), Range(1, 1));
                A_pi1[1]                     = sum(stop_densA_pi1);
                double PET1_pi0 = A_pi0[0] + R_pi0[0];
                double PET1_pi1 = A_pi1[0] + R_pi1[0];
                double ESS_pi0 = n1 + (1 - PET1_pi0)*n2;
                double ESS_pi1 = n1 + (1 - PET1_pi1)*n2;
                double Med_pi0 = (PET1_pi0 < 0.5 ?
                                    n1 + n2 : (PET1_pi0 > 0.5 ?
                                    n1 : n1 + 0.5*n2));
                double Med_pi1 = (PET1_pi1 < 0.5 ?
                                    n1 + n2 : (PET1_pi1 > 0.5 ?
                                    n1 : n1 + 0.5*n2));
                double Var_pi0 = PET1_pi0*pow(n1, 2) +
                  (1 - PET1_pi0)*pow(n1 + n2, 2) -
                  pow(ESS_pi0, 2);
                double Var_pi1 = PET1_pi1*pow(n1, 2) +
                  (1 - PET1_pi1)*pow(n1 + n2, 2) -
                  pow(ESS_pi1, 2);
                feasible_designs(counter, _) =
                  NumericVector::create(n1, n2, a1, a2, r1, sum(R_pi0),
                                        sum(R_pi1), PET1_pi0, PET1_pi1, ESS_pi0,
                                        ESS_pi1, Med_pi0, Med_pi1, Var_pi0,
                                        Var_pi1);
                counter++;
              }
            }
          }
        }
      }
    }
  }
  else {
    for (int n1 = (futility + efficacy == 2 ? 2 : 1);
         n1 <= (equal_n == 1 ? floor(Nmax/3) : (Nmax - 3)); n1++) {
      C_pi0(_, 0) = NumericVector(Nmax + 1);
      C_pi1(_, 0) = NumericVector(Nmax + 1);
      for (int s = 0; s <= n1; s++) {
        C_pi0(s, 0) = dbinomial_pi0(s, n1 - 1);
        C_pi1(s, 0) = dbinomial_pi1(s, n1 - 1);
      }
      for (int a1 = (futility == 1 ? 0 : -1);
           a1 <= (futility == 1 ? (ensign == 1 ?
                                     0 : (efficacy == 1 ?
                                     n1 - 2 : n1 - 1)) : -1);
           a1++) {
        if (a1 >= 0) {
          NumericMatrix stop_densA_pi0 = C_pi0(Range(0, a1), Range(0, 0));
          A_pi0[0]                     = sum(stop_densA_pi0);
          NumericMatrix stop_densA_pi1 = C_pi1(Range(0, a1), Range(0, 0));
          A_pi1[0]                     = sum(stop_densA_pi1);
        }
        else {
          A_pi0[0] = 0;
          A_pi1[0] = 0;
        }
        interrupt++;
        if (interrupt % 10000 == 0) {
          Rcpp::checkUserInterrupt();
          if (summary == 1) {
            Rcpp::Rcout << "...over " << interrupt  << " designs assessed..." << std::endl;
          }
        }
        if (A_pi1[0] > beta) {
          break;
        }
        for (int r1 = (efficacy == 1 ? n1 : n1 + 1);
             r1 >= (efficacy == 1 ? a1 + 2 : n1 + 1); r1--) {
          if (r1 <= n1) {
            NumericMatrix stop_densR_pi0 = C_pi0(Range(r1, n1), Range(0, 0));
            R_pi0[0]                     = sum(stop_densR_pi0);
            NumericMatrix stop_densR_pi1 = C_pi1(Range(r1, n1), Range(0, 0));
            R_pi1[0]                     = sum(stop_densR_pi1);
          }
          else {
            R_pi0[0] = 0;
            R_pi1[0] = 0;
          }
          if (R_pi0[0] > alpha){
            break;
          }
          for (int n2 = (equal_n == 1 ? n1 : 2);
               n2 <= (equal_n == 1 ? n1 : Nmax - n1 - 1); n2++) {
            C_pi0(_, 1) = NumericVector(Nmax + 1);
            C_pi1(_, 1) = NumericVector(Nmax + 1);
            for (int s1 = max(a1 + 1, 0); s1 <= min(r1 - 1, n1); s1++) {
              for (int s2 = 0; s2 <= n2; s2++) {
                C_pi0(s1 + s2, 1) = C_pi0(s1 + s2, 1) +
                  C_pi0(s1, 0)*dbinomial_pi0(s2, n2 - 1);
                C_pi1(s1 + s2, 1) = C_pi1(s1 + s2, 1) +
                  C_pi1(s1, 0)*dbinomial_pi1(s2, n2 - 1);
              }
            }
            for (int a2 = (futility == 1 ? max(0, a1 + 1) : -1);
                 a2 <= (futility == 1 ? min(r1 + n2 - 2, n1 + n2 - 1) : -1);
                 a2++){
              if (a2 >= 0){
                NumericMatrix stop_densA_pi0 = C_pi0(Range(0, a2), Range(1, 1));
                A_pi0[1]                     = sum(stop_densA_pi0);
                NumericMatrix stop_densA_pi1 = C_pi1(Range(0, a2), Range(1, 1));
                A_pi1[1]                     = sum(stop_densA_pi1);
              }
              else {
                A_pi0[1] = 0;
                A_pi1[1] = 0;
              }
              if (A_pi1[0] + A_pi1[1] > beta){
                break;
              }
              for (int r2 = (efficacy == 1 ? r1 + n2 - 1 : n1 + n2 + 1);
                   r2 >= (efficacy == 1 ? max(a2 + 2, r1) : n1 + n2 + 1);
                   r2--) {
                if (r2 <= n1 + n2){
                  NumericMatrix stop_densR_pi0 = C_pi0(Range(r2, n1 + n2 + 1),
                                                       Range(1, 1));
                  R_pi0[1]                     = sum(stop_densR_pi0);
                  NumericMatrix stop_densR_pi1 = C_pi1(Range(r2, n1 + n2 + 1),
                                                       Range(1, 1));
                  R_pi1[1]                     = sum(stop_densR_pi1);
                }
                else {
                  R_pi0[1] = 0;
                  R_pi1[1] = 0;
                }
                if (R_pi0[0] + R_pi0[1] > alpha){
                  break;
                }
                for (int n3 = (equal_n == 1 ? n1 : max(1, Nmin - n1 - n2));
                     n3 <= (equal_n == 1 ? n1 : Nmax - n1 - n2); n3++) {
                  C_pi0(_, 2) = NumericVector(Nmax + 1);
                  C_pi1(_, 2) = NumericVector(Nmax + 1);
                  for (int s2 = max(a2 + 1, 0); s2 <= min(r2 - 1, n1 + n2);
                  s2++) {
                    for (int s3 = 0; s3 <= n3; s3++) {
                      C_pi0(s2 + s3, 2) = C_pi0(s2 + s3, 2) +
                        C_pi0(s2, 1)*dbinomial_pi0(s3, n3 - 1);
                      C_pi1(s2 + s3, 2) = C_pi1(s2 + s3, 2) +
                        C_pi1(s2, 1)*dbinomial_pi1(s3, n3 - 1);
                    }
                  }
                  for (int a3 = max(0, a2 + 1);
                       a3 <= min(r2 + n3 - 2, n1 + n2 + n3 - 1); a3++) {
                    NumericMatrix stop_densR_pi1 = C_pi1(Range(a3 + 1,
                                                               n1 + n2 + n3),
                                                         Range(2, 2));
                    R_pi1[2]                     = sum(stop_densR_pi1);
                    if (sum(R_pi1) < 1 - beta){
                      break;
                    }
                    NumericMatrix stop_densR_pi0 = C_pi0(Range(a3 + 1,
                                                               n1 + n2 + n3),
                                                         Range(2, 2));
                    R_pi0[2]                     = sum(stop_densR_pi0);
                    if (sum(R_pi0) <= alpha){
                      NumericMatrix stop_densA_pi0 = C_pi0(Range(0, a3),
                                                           Range(2, 2));
                      A_pi0[2]                     = sum(stop_densA_pi0);
                      NumericMatrix stop_densA_pi1 = C_pi1(Range(0, a3),
                                                           Range(2, 2));
                      A_pi1[2]                     = sum(stop_densA_pi1);
                      double PET1_pi0 = A_pi0[0] + R_pi0[0];
                      double PET1_pi1 = A_pi1[0] + R_pi1[0];
                      double PET2_pi0 = A_pi0[1] + R_pi0[1];
                      double PET2_pi1 = A_pi1[1] + R_pi1[1];
                      double ESS_pi0 = PET1_pi0*n1 + PET2_pi0*(n1 + n2) +
                                         (1 - PET1_pi0 - PET2_pi0)*
                                         (n1 + n2 + n3);
                      double ESS_pi1 = PET1_pi1*n1 + PET2_pi1*(n1 + n2) +
                                         (1 - PET1_pi1 - PET2_pi1)*
                                         (n1 + n2 + n3);
                      double Med_pi0, Med_pi1;
                      if (PET1_pi0 == 0.5) {
                        Med_pi0 = n1 + 0.5*n2;
                      }
                      else if (PET1_pi0 + PET2_pi0 == 0.5){
                        Med_pi0 = n1 + n2 + 0.5*n3;
                      }
                      else {
                        if (PET1_pi0 > 0.5){
                          Med_pi0 = n1;
                        }
                        else if ((PET1_pi0 < 0.5) &&
                                   (PET1_pi0 + PET2_pi0 > 0.5)){
                          Med_pi0 = n1 + n2;
                        }
                        else {
                          Med_pi0 = n1 + n2 + n3;
                        }
                      }
                      if (PET1_pi1 == 0.5){
                        Med_pi1 = n1 + 0.5*n2;
                      }
                      else if (PET1_pi1 + PET2_pi1 == 0.5){
                        Med_pi1 = n1 + n2 + 0.5*n3;
                      }
                      else {
                        if (PET1_pi1 > 0.5){
                          Med_pi1 = n1;
                        }
                        else if ((PET1_pi1 < 0.5) &&
                                   (PET1_pi1 + PET2_pi1 > 0.5)){
                          Med_pi1 = n1 + n2;
                        }
                        else {
                          Med_pi1 = n1 + n2 + n3;
                        }
                      }
                      double Var_pi0 = PET1_pi0*pow(n1, 2) +
                        PET2_pi0*pow(n1 + n2, 2) + (1 - PET1_pi0 - PET2_pi0)*
                                                     pow(n1 + n2 + n3, 2) -
                        pow(ESS_pi0, 2);
                      double Var_pi1 = PET1_pi1*pow(n1, 2) +
                        PET2_pi1*pow(n1 + n2, 2) + (1 - PET1_pi1 - PET2_pi1)*
                                                     pow(n1 + n2 + n3, 2) -
                        pow(ESS_pi1, 2);
                      feasible_designs(counter, _) =
                        NumericVector::create(n1, n2, n3, a1, a2, a3, r1, r2,
                                              sum(R_pi0), sum(R_pi1), PET1_pi0,
                                              PET1_pi1, PET2_pi0, PET2_pi1,
                                              ESS_pi0, ESS_pi1, Med_pi0,
                                              Med_pi1, Var_pi0, Var_pi1);
                      counter++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  NumericMatrix output = feasible_designs(Range(0, 0 + (counter > 0 ?
                                                          counter - 1 : 0)),
                                          Range(0, (J == 2 ? 14 : 19)));
  return output;
}
