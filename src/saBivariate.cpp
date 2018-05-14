#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix saBivariateBD(int J, double alphaR, double alphaT, double beta,
                            int Nmin, int Nmax, int equal_n,
                            Rcpp::List dbivar_00_list,
                            Rcpp::List dbivar_01_list,
                            Rcpp::List dbivar_10_list,
                            Rcpp::List dbivar_11_min_list,
                            Rcpp::List dbivar_11_max_list, int summary) {
  int counter = 0;
  int interrupt = 0;
  if (J == 1) {
    NumericMatrix feasible_designs(10000000, 8);
    double R_00, R_01, R_10, R_11_min, R_11_max;
    for (int n = Nmin; n <= Nmax; n++) {
      Rcpp::NumericMatrix dbivar_00 = dbivar_00_list[n - 1];
      Rcpp::NumericMatrix dbivar_01 = dbivar_01_list[n - 1];
      Rcpp::NumericMatrix dbivar_10 = dbivar_10_list[n - 1];
      Rcpp::NumericMatrix dbivar_11_min = dbivar_11_min_list[n - 1];
      Rcpp::NumericMatrix dbivar_11_max = dbivar_11_max_list[n - 1];
      for (int aR = 0; aR < n; aR++) {
        for (int aT = n - 1; aT > 0; aT--) {
          interrupt++;
          if (interrupt % 10000 == 0) {
            Rcpp::checkUserInterrupt();
            if (summary == 1) {
              Rcpp::Rcout << "Over " << interrupt  << " designs assessed..." << std::endl;
            }
          }
          NumericMatrix reject_01       = dbivar_01(Range(aR + 1, n),
                                                    Range(0, aT - 1));
          R_01                          = sum(reject_01);
          if (R_01 <= alphaR) {
            NumericMatrix reject_10     = dbivar_10(Range(aR + 1, n),
                                                    Range(0, aT - 1));
            R_10                        = sum(reject_10);
            if (R_10 <= alphaT) {
              NumericMatrix reject_11   = dbivar_11_min(Range(aR + 1, n),
                                                        Range(0, aT - 1));
              R_11_min                  = sum(reject_11);
              if (R_11_min >= 1 - beta) {
                NumericMatrix reject_00 = dbivar_00(Range(aR + 1, n),
                                                    Range(0, aT - 1));
                R_00                    = sum(reject_00);
                NumericMatrix reject_11   = dbivar_11_max(Range(aR + 1, n),
                                                          Range(0, aT - 1));
                R_11_max                  = sum(reject_11);
                feasible_designs(counter, _) = NumericVector::create(n, aR,
                                 aT, R_00, R_01, R_10, R_11_min, R_11_max);
                counter++;
              }
            }
          }
        }
      }
    }
    NumericMatrix output = feasible_designs(Range(0, 0 + (counter > 0 ?
                                                            counter - 1 : 0)),
                                                            Range(0, 7));
    return output;
  }
  else {
    NumericMatrix feasible_designs_1(10000000, 11), feasible_designs_2(10000000, 20);
    NumericVector A_00(2), A_01(2), A_10(2), A_11_min(2), A_11_max(2);
    for (int n1 = (equal_n == 1 ? ceil(Nmin/2) : 1);
         n1 <= (equal_n == 1 ? floor(Nmax/2) : (Nmax - 1)); n1++) {
      Rcpp::NumericMatrix dbivar_00_1     = dbivar_00_list[n1 - 1];
      Rcpp::NumericMatrix dbivar_01_1     = dbivar_01_list[n1 - 1];
      Rcpp::NumericMatrix dbivar_10_1     = dbivar_10_list[n1 - 1];
      Rcpp::NumericMatrix dbivar_11_min_1 = dbivar_11_min_list[n1 - 1];
      Rcpp::NumericMatrix dbivar_11_max_1 = dbivar_11_max_list[n1 - 1];
      for (int aR1 = 0; aR1 < n1; aR1++) {
        for (int aT1 = n1 - 1; aT1 > 0; aT1--) {
          NumericMatrix reject_11_min = dbivar_11_min_1(Range(aR1 + 1, n1),
                                                        Range(0, aT1 - 1));
          A_11_min[0] = 1 - sum(reject_11_min);
          if (A_11_min[0] > beta) {
            break;
          }
          NumericMatrix reject_00 = dbivar_00_1(Range(aR1 + 1, n1),
                                                Range(0, aT1 - 1));
          A_00[0] = 1 - sum(reject_00);
          NumericMatrix reject_01 = dbivar_01_1(Range(aR1 + 1, n1),
                                                Range(0, aT1 - 1));
          A_01[0] = 1 - sum(reject_01);
          NumericMatrix reject_10 = dbivar_10_1(Range(aR1 + 1, n1),
                                                Range(0, aT1 - 1));
          A_10[0] = 1 - sum(reject_10);
          NumericMatrix reject_11 = dbivar_11_max_1(Range(aR1 + 1, n1),
                                                    Range(0, aT1 - 1));
          A_11_max[0] = 1 - sum(reject_11);
          for (int n2 = (equal_n == 1 ? n1 : max(1, Nmin - n1));
               n2 <= (equal_n == 1 ? n1 : Nmax - n1); n2++) {
            Rcpp::NumericMatrix dbivar_00_2    = dbivar_00_list[n2 - 1];
            Rcpp::NumericMatrix dbivar_01_2    = dbivar_01_list[n2 - 1];
            Rcpp::NumericMatrix dbivar_10_2    = dbivar_10_list[n2 - 1];
            Rcpp::NumericMatrix dbivar_11_min_2 = dbivar_11_min_list[n2 - 1];
            Rcpp::NumericMatrix dbivar_11_max_2 = dbivar_11_max_list[n2 - 1];
            Rcpp::NumericMatrix pdf_2_00(n1 + n2 + 1, n1 + n2 + 1);
            Rcpp::NumericMatrix pdf_2_01(n1 + n2 + 1, n1 + n2 + 1);
            Rcpp::NumericMatrix pdf_2_10(n1 + n2 + 1, n1 + n2 + 1);
            Rcpp::NumericMatrix pdf_2_11_min(n1 + n2 + 1, n1 + n2 + 1);
            Rcpp::NumericMatrix pdf_2_11_max(n1 + n2 + 1, n1 + n2 + 1);
            for (int sR1 = aR1 + 1; sR1 <= n1; sR1++) {
              for (int sT1 = 0; sT1 <= aT1 - 1; sT1++) {
                for (int sR2 = 0; sR2 <= n2; sR2++) {
                  for (int sT2 = 0; sT2 <= n2; sT2++) {
                    pdf_2_00(sR1 + sR2, sT1 + sT2) = pdf_2_00(sR1 + sR2, sT1 + sT2) + dbivar_00_1(sR1, sT1)*dbivar_00_2(sR2, sT2);
                    pdf_2_01(sR1 + sR2, sT1 + sT2) = pdf_2_01(sR1 + sR2, sT1 + sT2) + dbivar_01_1(sR1, sT1)*dbivar_01_2(sR2, sT2);
                    pdf_2_10(sR1 + sR2, sT1 + sT2) = pdf_2_10(sR1 + sR2, sT1 + sT2) + dbivar_10_1(sR1, sT1)*dbivar_10_2(sR2, sT2);
                    pdf_2_11_min(sR1 + sR2, sT1 + sT2) = pdf_2_11_min(sR1 + sR2, sT1 + sT2) + dbivar_11_min_1(sR1, sT1)*dbivar_11_min_2(sR2, sT2);
                    pdf_2_11_max(sR1 + sR2, sT1 + sT2) = pdf_2_11_max(sR1 + sR2, sT1 + sT2) + dbivar_11_max_1(sR1, sT1)*dbivar_11_max_2(sR2, sT2);
                  }
                }
              }
            }
            for (int aR2 = aR1 + 1; aR2 < n1 + n2; aR2++) {
              for (int aT2 = aT1 + n2 - 1; aT2 > 0; aT2--) {
                interrupt++;
                if (interrupt % 10000 == 0) {
                  Rcpp::checkUserInterrupt();
                  if (summary == 1) {
                    Rcpp::Rcout << "Over " << interrupt  << " designs assessed..." << std::endl;
                  }
                }

                NumericMatrix reject_11_min = pdf_2_11_min(Range(aR2 + 1, n1 + n2),
                                                           Range(0, aT2 - 1));
                A_11_min[1] = sum(pdf_2_11_min) - sum(reject_11_min);
                if (A_11_min[0] + A_11_min[1] > beta) {
                  break;
                }
                NumericMatrix reject_01 = pdf_2_01(Range(aR2 + 1, n1 + n2),
                                                   Range(0, aT2 - 1));
                A_01[1] = sum(pdf_2_01) - sum(reject_01);
                if (1 - A_01[0] - A_01[1] <= alphaR) {
                  NumericMatrix reject_10 = pdf_2_10(Range(aR2 + 1, n1 + n2),
                                                     Range(0, aT2 - 1));
                  A_10[1] = sum(pdf_2_10) - sum(reject_10);
                  if (1 - A_10[0] - A_10[1] <= alphaT) {
                    NumericMatrix reject_00 = pdf_2_00(Range(aR2 + 1, n1 + n2),
                                                       Range(0, aT2 - 1));
                    A_00[1] = sum(pdf_2_00) - sum(reject_00);
                    NumericMatrix reject_11 = pdf_2_11_min(Range(aR2 + 1, n1 + n2),
                                                              Range(0, aT2 - 1));
                    A_11_max[1] = sum(pdf_2_11_min) - sum(reject_11);
                    double PET1_00 = A_00[0];
                    double PET1_01 = A_01[0];
                    double PET1_10 = A_10[0];
                    double PET1_11_min = A_11_min[0];
                    double PET1_11_max = A_11_max[0];
                    double ESS_00 = n1 + (1 - PET1_00)*n2;
                    double ESS_01 = n1 + (1 - PET1_01)*n2;
                    double ESS_10 = n1 + (1 - PET1_10)*n2;
                    double ESS_11_min = n1 + (1 - PET1_11_min)*n2;
                    double ESS_11_max = n1 + (1 - PET1_11_max)*n2;
                    double Med_00 = (PET1_00 < 0.5 ?
                                        n1 + n2 : (PET1_00 > 0.5 ?
                                        n1 : n1 + 0.5*n2));
                    double Med_01 = (PET1_01 < 0.5 ?
                                       n1 + n2 : (PET1_01 > 0.5 ?
                                       n1 : n1 + 0.5*n2));
                    double Med_10 = (PET1_10 < 0.5 ?
                                       n1 + n2 : (PET1_10 > 0.5 ?
                                       n1 : n1 + 0.5*n2));
                    double Med_11_min = (PET1_11_min < 0.5 ?
                                       n1 + n2 : (PET1_11_min > 0.5 ?
                                       n1 : n1 + 0.5*n2));
                    double Med_11_max = (PET1_11_max < 0.5 ?
                                       n1 + n2 : (PET1_11_max > 0.5 ?
                                       n1 : n1 + 0.5*n2));
                    double Var_00 = PET1_00*pow(n1, 2) +
                                      (1 - PET1_00)*pow(n1 + n2, 2) -
                                      pow(ESS_00, 2);
                    double Var_01 = PET1_01*pow(n1, 2) +
                                      (1 - PET1_01)*pow(n1 + n2, 2) -
                                      pow(ESS_01, 2);
                    double Var_10 = PET1_10*pow(n1, 2) +
                                      (1 - PET1_10)*pow(n1 + n2, 2) -
                                      pow(ESS_10, 2);
                    double Var_11_min = PET1_11_min*pow(n1, 2) +
                                      (1 - PET1_11_min)*pow(n1 + n2, 2) -
                                      pow(ESS_11_min, 2);
                    double Var_11_max = PET1_11_max*pow(n1, 2) +
                                      (1 - PET1_11_max)*pow(n1 + n2, 2) -
                                      pow(ESS_11_max, 2);
                    feasible_designs_1(counter, _) = NumericVector::create(n1, n2, aR1,
                                     aR2, aT1, aT2, 1 - sum(A_00), 1 - sum(A_01),
                                     1 - sum(A_10), 1 - sum(A_11_min), 1 - sum(A_11_max));
                    feasible_designs_2(counter, _) = NumericVector::create(PET1_00,
                                       PET1_01, PET1_10, PET1_11_min, PET1_11_max,
                                     ESS_00, ESS_01, ESS_10, ESS_11_min, ESS_11_max,
                                     Med_00, Med_01, Med_10, Med_11_min, Med_11_max,
                                     Var_00, Var_01, Var_10, Var_11_min, Var_11_max);
                    counter++;
                  }
                }
              }
            }
          }
        }
      }
    }
    NumericMatrix feasible_designs = Rcpp::cbind(feasible_designs_1, feasible_designs_2);
    NumericMatrix output = feasible_designs(Range(0, 0 + (counter > 0 ?
                                                            counter - 1 : 0)),
                                                            Range(0, 30));
    return output;
  }
}

// [[Rcpp::export]]
NumericMatrix saBivariateCP(int J, double alphaL, double alphaG, double beta,
                            int Nmin, int Nmax, int equal_n,
                            Rcpp::List dbivar_00_list,
                            Rcpp::List dbivar_01_list,
                            Rcpp::List dbivar_10_list,
                            Rcpp::List dbivar_11_list,
                            Rcpp::List dbivar_max1_list,
                            Rcpp::List dbivar_max2_list, int summary) {
  int counter = 0;
  int interrupt = 0;
  if (J == 1) {
    NumericMatrix feasible_designs(10000000, 8);
    double R_00, R_01, R_10, R_11, R_max;
    for (int n = Nmin; n <= Nmax; n++) {
      Rcpp::NumericMatrix dbivar_00 = dbivar_00_list[n - 1];
      Rcpp::NumericMatrix dbivar_01 = dbivar_01_list[n - 1];
      Rcpp::NumericMatrix dbivar_10 = dbivar_10_list[n - 1];
      Rcpp::NumericMatrix dbivar_11 = dbivar_11_list[n - 1];
      Rcpp::NumericMatrix dbivar_max1 = dbivar_max1_list[n - 1];
      Rcpp::NumericMatrix dbivar_max2 = dbivar_max2_list[n - 1];
      for (int aR = 0; aR < n; aR++) {
        for (int aT = n - 1; aT > 0; aT--) {
          interrupt++;
          if (interrupt % 10000 == 0) {
            Rcpp::checkUserInterrupt();
            if (summary == 1) {
              Rcpp::Rcout << "Over " << interrupt  << " designs assessed..." << std::endl;
            }
          }
          NumericMatrix reject_00       = dbivar_00(Range(aR + 1, n),
                                                    Range(0, aT - 1));
          R_00                          = sum(reject_00);
          if (R_00 <= alphaL) {
            NumericMatrix reject_11     = dbivar_11(Range(aR + 1, n),
                                                    Range(0, aT - 1));
            R_11                        = sum(reject_11);
            if (R_11 >= 1 - beta) {
              NumericMatrix reject_max1 = dbivar_max1(Range(aR + 1, n),
                                                        Range(0, aT - 1));
              NumericMatrix reject_max2 = dbivar_max2(Range(aR + 1, n),
                                                      Range(0, aT - 1));
              R_max                     = max(sum(reject_max1), sum(reject_max2));
              if (R_max <= alphaG) {
                NumericMatrix reject_01 = dbivar_01(Range(aR + 1, n),
                                                          Range(0, aT - 1));
                R_01                    = sum(reject_01);
                NumericMatrix reject_10 = dbivar_10(Range(aR + 1, n),
                                                    Range(0, aT - 1));
                R_10                    = sum(reject_10);
                feasible_designs(counter, _) = NumericVector::create(n, aR,
                                 aT, R_00, R_max, R_01, R_10, R_11);
                counter++;
              }
            }
          }
        }
      }
    }
    NumericMatrix output = feasible_designs(Range(0, 0 + (counter > 0 ?
                                                            counter - 1 : 0)),
                                                            Range(0, 7));
    return output;
  }
  else {
    NumericMatrix feasible_designs_1(20000000, 18), feasible_designs_2(20000000, 18);
    NumericVector A_00(2), A_01(2), A_10(2), A_11(2), A_max1(2), A_max2(2);
    for (int n1 = (equal_n == 1 ? ceil(Nmin/2) : 1);
         n1 <= (equal_n == 1 ? floor(Nmax/2) : (Nmax - 1)); n1++) {
      Rcpp::NumericMatrix dbivar_00_1     = dbivar_00_list[n1 - 1];
      Rcpp::NumericMatrix dbivar_01_1     = dbivar_01_list[n1 - 1];
      Rcpp::NumericMatrix dbivar_10_1     = dbivar_10_list[n1 - 1];
      Rcpp::NumericMatrix dbivar_11_1     = dbivar_11_list[n1 - 1];
      Rcpp::NumericMatrix dbivar_max1_1 = dbivar_max1_list[n1 - 1];
      Rcpp::NumericMatrix dbivar_max2_1 = dbivar_max2_list[n1 - 1];
      for (int aR1 = 0; aR1 < n1; aR1++) {
        for (int aT1 = n1 - 1; aT1 > 0; aT1--) {
          NumericMatrix reject_11 = dbivar_11_1(Range(aR1 + 1, n1),
                                                Range(0, aT1 - 1));
          A_11[0] = 1 - sum(reject_11);
          if (A_11[0] > beta) {
            break;
          }
          NumericMatrix reject_00 = dbivar_00_1(Range(aR1 + 1, n1),
                                                Range(0, aT1 - 1));
          A_00[0] = 1 - sum(reject_00);
          NumericMatrix reject_01 = dbivar_01_1(Range(aR1 + 1, n1),
                                                Range(0, aT1 - 1));
          A_01[0] = 1 - sum(reject_01);
          NumericMatrix reject_10 = dbivar_10_1(Range(aR1 + 1, n1),
                                                Range(0, aT1 - 1));
          A_10[0] = 1 - sum(reject_10);
          NumericMatrix reject_max1 = dbivar_max1_1(Range(aR1 + 1, n1),
                                                    Range(0, aT1 - 1));
          A_max1[0] = 1 - sum(reject_max1);
          NumericMatrix reject_max2 = dbivar_max2_1(Range(aR1 + 1, n1),
                                                    Range(0, aT1 - 1));
          A_max2[0] = 1 - sum(reject_max2);
          for (int n2 = (equal_n == 1 ? n1 : max(1, Nmin - n1));
               n2 <= (equal_n == 1 ? n1 : Nmax - n1); n2++) {
            Rcpp::NumericMatrix dbivar_00_2    = dbivar_00_list[n2 - 1];
            Rcpp::NumericMatrix dbivar_01_2    = dbivar_01_list[n2 - 1];
            Rcpp::NumericMatrix dbivar_10_2    = dbivar_10_list[n2 - 1];
            Rcpp::NumericMatrix dbivar_11_2    = dbivar_11_list[n2 - 1];
            Rcpp::NumericMatrix dbivar_max1_2 = dbivar_max1_list[n2 - 1];
            Rcpp::NumericMatrix dbivar_max2_2 = dbivar_max2_list[n2 - 1];
            Rcpp::NumericMatrix pdf_2_00(n1 + n2 + 1, n1 + n2 + 1);
            Rcpp::NumericMatrix pdf_2_01(n1 + n2 + 1, n1 + n2 + 1);
            Rcpp::NumericMatrix pdf_2_10(n1 + n2 + 1, n1 + n2 + 1);
            Rcpp::NumericMatrix pdf_2_11(n1 + n2 + 1, n1 + n2 + 1);
            Rcpp::NumericMatrix pdf_2_max1(n1 + n2 + 1, n1 + n2 + 1);
            Rcpp::NumericMatrix pdf_2_max2(n1 + n2 + 1, n1 + n2 + 1);
            for (int sR1 = aR1 + 1; sR1 <= n1; sR1++) {
              for (int sT1 = 0; sT1 <= aT1 - 1; sT1++) {
                for (int sR2 = 0; sR2 <= n2; sR2++) {
                  for (int sT2 = 0; sT2 <= n2; sT2++) {
                    pdf_2_00(sR1 + sR2, sT1 + sT2) = pdf_2_00(sR1 + sR2, sT1 + sT2) + dbivar_00_1(sR1, sT1)*dbivar_00_2(sR2, sT2);
                    pdf_2_01(sR1 + sR2, sT1 + sT2) = pdf_2_01(sR1 + sR2, sT1 + sT2) + dbivar_01_1(sR1, sT1)*dbivar_01_2(sR2, sT2);
                    pdf_2_10(sR1 + sR2, sT1 + sT2) = pdf_2_10(sR1 + sR2, sT1 + sT2) + dbivar_10_1(sR1, sT1)*dbivar_10_2(sR2, sT2);
                    pdf_2_11(sR1 + sR2, sT1 + sT2) = pdf_2_11(sR1 + sR2, sT1 + sT2) + dbivar_11_1(sR1, sT1)*dbivar_11_2(sR2, sT2);
                    pdf_2_max1(sR1 + sR2, sT1 + sT2) = pdf_2_max1(sR1 + sR2, sT1 + sT2) + dbivar_max1_1(sR1, sT1)*dbivar_max1_2(sR2, sT2);
                    pdf_2_max2(sR1 + sR2, sT1 + sT2) = pdf_2_max2(sR1 + sR2, sT1 + sT2) + dbivar_max2_1(sR1, sT1)*dbivar_max2_2(sR2, sT2);
                  }
                }
              }
            }
            for (int aR2 = aR1 + 1; aR2 < n1 + n2; aR2++) {
              for (int aT2 = aT1 + n2 - 1; aT2 > 0; aT2--) {
                interrupt++;
                if (interrupt % 10000 == 0) {
                  Rcpp::checkUserInterrupt();
                  if (summary == 1) {
                    Rcpp::Rcout << "Over " << interrupt  << " designs assessed..." << std::endl;
                  }
                }
                NumericMatrix reject_11 = pdf_2_11(Range(aR2 + 1, n1 + n2),
                                                   Range(0, aT2 - 1));
                A_11[1] = sum(pdf_2_11) - sum(reject_11);
                if (A_11[0] + A_11[1] > beta) {
                  break;
                }
                NumericMatrix reject_00 = pdf_2_00(Range(aR2 + 1, n1 + n2),
                                                   Range(0, aT2 - 1));
                A_00[1] = sum(pdf_2_00) - sum(reject_00);
                if (1 - A_00[0] - A_00[1] <= alphaL) {
                  NumericMatrix reject_max1 = pdf_2_max1(Range(aR2 + 1, n1 + n2),
                                                         Range(0, aT2 - 1));
                  A_max1[1] = sum(pdf_2_max1) - sum(reject_max1);
                  NumericMatrix reject_max2 = pdf_2_max2(Range(aR2 + 1, n1 + n2),
                                                         Range(0, aT2 - 1));
                  A_max2[1] = sum(pdf_2_max2) - sum(reject_max2);
                  if (1 - min(sum(A_max1), sum(A_max2)) <= alphaG) {
                    NumericMatrix reject_01 = pdf_2_01(Range(aR2 + 1, n1 + n2),
                                                       Range(0, aT2 - 1));
                    A_01[1] = sum(pdf_2_01) - sum(reject_01);
                    NumericMatrix reject_10 = pdf_2_10(Range(aR2 + 1, n1 + n2),
                                                       Range(0, aT2 - 1));
                    A_10[1] = sum(pdf_2_10) - sum(reject_10);
                    double PET1_00 = A_00[0];
                    double PET1_01 = A_01[0];
                    double PET1_10 = A_10[0];
                    double PET1_11 = A_11[0];
                    double PET1_max1 = A_max1[0];
                    double PET1_max2 = A_max2[0];
                    double ESS_00 = n1 + (1 - PET1_00)*n2;
                    double ESS_01 = n1 + (1 - PET1_01)*n2;
                    double ESS_10 = n1 + (1 - PET1_10)*n2;
                    double ESS_11 = n1 + (1 - PET1_11)*n2;
                    double ESS_max1 = n1 + (1 - PET1_max1)*n2;
                    double ESS_max2 = n1 + (1 - PET1_max2)*n2;
                    double Med_00 = (PET1_00 < 0.5 ?
                                       n1 + n2 : (PET1_00 > 0.5 ?
                                       n1 : n1 + 0.5*n2));
                    double Med_01 = (PET1_01 < 0.5 ?
                                       n1 + n2 : (PET1_01 > 0.5 ?
                                       n1 : n1 + 0.5*n2));
                    double Med_10 = (PET1_10 < 0.5 ?
                                       n1 + n2 : (PET1_10 > 0.5 ?
                                       n1 : n1 + 0.5*n2));
                    double Med_11 = (PET1_11 < 0.5 ?
                                           n1 + n2 : (PET1_11 > 0.5 ?
                                           n1 : n1 + 0.5*n2));
                    double Med_max1 = (PET1_max1 < 0.5 ?
                                           n1 + n2 : (PET1_max1 > 0.5 ?
                                           n1 : n1 + 0.5*n2));
                    double Med_max2 = (PET1_max2 < 0.5 ?
                                         n1 + n2 : (PET1_max2 > 0.5 ?
                                         n1 : n1 + 0.5*n2));
                    double Var_00 = PET1_00*pow(n1, 2) +
                      (1 - PET1_00)*pow(n1 + n2, 2) -
                      pow(ESS_00, 2);
                    double Var_01 = PET1_01*pow(n1, 2) +
                      (1 - PET1_01)*pow(n1 + n2, 2) -
                      pow(ESS_01, 2);
                    double Var_10 = PET1_10*pow(n1, 2) +
                      (1 - PET1_10)*pow(n1 + n2, 2) -
                      pow(ESS_10, 2);
                    double Var_11 = PET1_11*pow(n1, 2) +
                      (1 - PET1_11)*pow(n1 + n2, 2) -
                      pow(ESS_11, 2);
                    double Var_max1 = PET1_max1*pow(n1, 2) +
                      (1 - PET1_max1)*pow(n1 + n2, 2) -
                      pow(ESS_max1, 2);
                    double Var_max2 = PET1_max2*pow(n1, 2) +
                      (1 - PET1_max2)*pow(n1 + n2, 2) -
                      pow(ESS_max2, 2);
                    feasible_designs_1(counter, _) = NumericVector::create(n1, n2, aR1,
                                       aR2, aT1, aT2, 1 - sum(A_00), 1 - sum(A_01),
                                       1 - sum(A_10), 1 - sum(A_11), 1 - sum(A_max1),
                                       1 - sum(A_max2), PET1_00,
                                       PET1_01, PET1_10, PET1_11, PET1_max1, PET1_max2);
                    feasible_designs_2(counter, _) = NumericVector::create(ESS_00, ESS_01, ESS_10, ESS_11, ESS_max1, ESS_max2,
                                       Med_00, Med_01, Med_10, Med_11, Med_max1, Med_max2,
                                       Var_00, Var_01, Var_10, Var_11, Var_max1, Var_max2);
                    counter++;
                  }
                }
              }
            }
          }
        }
      }
    }
    NumericMatrix feasible_designs = Rcpp::cbind(feasible_designs_1, feasible_designs_2);
    NumericMatrix output = feasible_designs(Range(0, 0 + (counter > 0 ?
                                                            counter - 1 : 0)),
                                                            Range(0, 35));
    return output;
  }
}

